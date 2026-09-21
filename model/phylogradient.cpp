/*
 * phylogradient.cpp
 *
 * See phylogradient.h and docs/analytical-gradients-design.md.
 */

#include "phylogradient.h"
#include "modelfactory.h"
#include "modelsubst.h"
#include "rateheterogeneity.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

/* ---------------------------------------------------------------------- */
/* RAII guard                                                              */
/* ---------------------------------------------------------------------- */

PhyloGradient::OutsideGuard::OutsideGuard(PhyloNeighbor *n, double *buf, UBYTE *scale)
    : nei(n), saved_lh(n->partial_lh), saved_scale(n->scale_num), saved_computed(n->partial_lh_computed) {
    // NOTE(design 5.2): attach BEFORE computeTraversalInfo runs; with partial_lh
    // set, reorientPartialLh() returns early and no inside slot is taken over.
    ASSERT(saved_lh == nullptr && "reverse neighbour already owns a partial");
    nei->partial_lh = buf;
    nei->scale_num = scale;
    nei->partial_lh_computed &= ~1;
}

PhyloGradient::OutsideGuard::~OutsideGuard() noexcept {
    nei->partial_lh = saved_lh;
    nei->scale_num = saved_scale;
    nei->partial_lh_computed = saved_computed;
}

/* ---------------------------------------------------------------------- */
/* construction / buffers                                                  */
/* ---------------------------------------------------------------------- */

PhyloGradient::PhyloGradient(PhyloTree *t) : tree(t) {
}

PhyloGradient::~PhyloGradient() {
    for (auto p : pool_lh) aligned_free(p);
    for (auto p : pool_scale) aligned_free(p);
}

void PhyloGradient::ensurePool(int depth) {
    while ((int)pool_lh.size() <= depth) {
        size_t lh_size = tree->getPartialLhSize();
        size_t sc_size = tree->getScaleNumSize();
        double *lh = aligned_alloc<double>(lh_size);
        UBYTE *sc = aligned_alloc<UBYTE>(sc_size);
        memset(sc, 0, sc_size * sizeof(UBYTE));
        pool_lh.push_back(lh);
        pool_scale.push_back(sc);
        buffer_bytes += lh_size * sizeof(double) + sc_size * sizeof(UBYTE);
    }
}

int PhyloGradient::treeDepth(PhyloNode *node, PhyloNode *dad) const {
    int d = 0;
    FOR_NEIGHBOR_DECLARE(node, dad, it) {
        PhyloNode *child = (PhyloNode*)(*it)->node;
        if (!child->isLeaf())
            d = max(d, 1 + treeDepth(child, node));
    }
    return d;
}

void PhyloGradient::cacheDimensions() {
    ModelSubst *model = tree->getModel();
    RateHeterogeneity *site_rate = tree->getRate();
    ModelFactory *mf = tree->getModelFactory();
    nstates = tree->aln->num_states;
    ncat = site_rate->getNRate();
    nmix = model->getNMixtures();
    fused = mf->fused_mix_rate;
    ncat_mix = fused ? ncat : ncat * nmix;
    denom = fused ? 1 : ncat;
    block = nstates * ncat_mix;
    tip_block = nstates * nmix;
    vsize = tree->vector_size;
    orig_nptn = tree->aln->size();
    nptn = get_safe_upper_limit(orig_nptn);
    safe = tree->safe_numeric;

    class_rate.assign(ncat_mix, 0.0);
    class_weight.assign(ncat_mix, 0.0);
    class_eval_offset.assign(ncat_mix, 0);
    class_mix.assign(ncat_mix, 0);
    size_t eval_stride = get_safe_upper_limit(nstates);   // padded per component (ModelMixture::initMem)
    for (size_t c = 0; c < ncat_mix; c++) {
        size_t mycat = c % ncat;
        size_t m = c / denom;
        class_rate[c] = site_rate->getRate(mycat);
        class_weight[c] = site_rate->getProp(mycat) * model->getMixtureWeight(m);
        class_eval_offset[c] = m * eval_stride;
        class_mix[c] = m;
    }
    val.assign(ncat_mix * nstates, 0.0);
    lval.assign(ncat_mix * nstates, 0.0);
}

/* ---------------------------------------------------------------------- */
/* main entry                                                              */
/* ---------------------------------------------------------------------- */

PhyloGradient::Side PhyloGradient::sideOf(PhyloNeighbor *nei_to_side, PhyloNode *side_node) const {
    Side s;
    if (side_node->isLeaf()) {
        s.leaf_id = side_node->id;
    } else {
        s.lh = nei_to_side->partial_lh;
        s.scale = nei_to_side->scale_num;
        ASSERT(s.lh && "inside partial not resident");
    }
    return s;
}

bool PhyloGradient::compute(Result &res) {
    res = Result();
    // forward pass: every inside partial is now oriented toward current_it
    tree->clearAllPartialLH();
    res.logl = tree->computeLikelihood();

    PhyloNeighbor *saved_it = tree->current_it;
    PhyloNeighbor *saved_back = tree->current_it_back;
    ASSERT(saved_it && saved_back);

    cacheDimensions();
    res.dlogl_dt.assign(max(tree->branchNum, 1), 0.0);

    PhyloNode *v = (PhyloNode*)saved_it->node;         // dad_branch side
    PhyloNode *u = (PhyloNode*)saved_back->node;       // dad side

    // root edge (u,v): both directions are resident
    {
        Side side_v = sideOf(saved_it, v);
        Side side_u = sideOf(saved_back, u);
        accumulateEdge(saved_it, side_v, side_u, res);
    }

    // buffer pool sized from the actual recursion depth (a caterpillar needs N-2)
    int depth = max(u->isLeaf() ? 0 : treeDepth(u, v), v->isLeaf() ? 0 : treeDepth(v, u));
    ensurePool(depth);

    if (!u->isLeaf()) visit(u, v, 0, res);
    if (!v->isLeaf()) visit(v, u, 0, res);

    // restore tree state the pass may have disturbed
    tree->current_it = saved_it;
    tree->current_it_back = saved_back;
    tree->theta_computed = false;   // NOTE(design 5.2): kernel sets it per edge; stale theta_all must not be reused

    for (double g : res.dlogl_dt)
        if (!std::isfinite(g)) res.valid = false;
    if (!std::isfinite(res.logl)) res.valid = false;
    return res.valid;
}

void PhyloGradient::visit(PhyloNode *node, PhyloNode *dad, int depth, Result &res) {
    FOR_NEIGHBOR_DECLARE(node, dad, it) {
        PhyloNeighbor *nei = (PhyloNeighbor*)(*it);          // node -> child (inside partial of child)
        PhyloNode *child = (PhyloNode*)nei->node;
        PhyloNeighbor *rev = (PhyloNeighbor*)child->findNeighbor(node);   // child -> node (outside partial)

        OutsideGuard guard(rev, poolLh(depth), poolScale(depth));
        // The kernel computes rev (node's side away from child) from node's other
        // neighbours: the parent-side reverse neighbour attached one level up (or
        // the resident root-edge partial) and the resident sibling partials, and
        // returns the tree lnL on this edge.
        double edge_logl = tree->computeLikelihoodBranch(rev, child, true);
        res.max_edge_logl_diff = max(res.max_edge_logl_diff, fabs(edge_logl - res.logl));

        Side side_child = sideOf(nei, child);
        Side side_parent;
        side_parent.lh = rev->partial_lh;
        side_parent.scale = rev->scale_num;
        accumulateEdge(nei, side_child, side_parent, res);

        if (!child->isLeaf())
            visit(child, node, depth + 1, res);
        // guard restores rev on scope exit
    }
}

/* ---------------------------------------------------------------------- */
/* per-edge accumulation                                                   */
/* ---------------------------------------------------------------------- */

void PhyloGradient::accumulateEdge(PhyloNeighbor *nei, const Side &child, const Side &parent, Result &res) {
    ModelSubst *model = tree->getModel();
    const double *eval = model->getEigenvalues();
    const double *tip_lh = tree->tip_partial_lh;
    const double *ptn_freq = tree->ptn_freq;
    const double *ptn_invar = tree->ptn_invar;
    const size_t V = vsize;

    // per-class exp(lambda * rate * t) and lambda * exp(...)
    for (size_t c = 0; c < ncat_mix; c++) {
        double tau = class_rate[c] * nei->getLength(c % ncat);
        const double *ev = eval + class_eval_offset[c];
        double *vc = &val[c * nstates], *lc = &lval[c * nstates];
        for (size_t i = 0; i < nstates; i++) {
            vc[i] = exp(ev[i] * tau);
            lc[i] = ev[i] * vc[i];
        }
    }

    // pattern chunks: contiguous, aligned to the SIMD block, one partial sum each
    int nchunks = max(1, tree->num_packets);
    size_t nblocks = (orig_nptn + V - 1) / V;
    vector<double> chunk_sum(nchunks, 0.0);
    vector<char> chunk_bad(nchunks, 0);

    const Alignment *aln = tree->aln;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(tree->num_threads)
#endif
    for (int ch = 0; ch < nchunks; ch++) {
        size_t b_lo = (nblocks * ch) / nchunks, b_hi = (nblocks * (ch + 1)) / nchunks;
        double sum = 0.0;
        vector<double> a_c(ncat_mix), b_c(ncat_mix);
        vector<unsigned> sc(ncat_mix);
        for (size_t ptn = b_lo * V; ptn < min(b_hi * V, orig_nptn); ptn++) {
            size_t base = (ptn / V) * V * block, lane = ptn % V;
            int st_child = child.leaf_id >= 0 ? (*aln)[ptn][child.leaf_id] : -1;
            int st_parent = parent.leaf_id >= 0 ? (*aln)[ptn][parent.leaf_id] : -1;
            unsigned min_scale = 0;
            for (size_t c = 0; c < ncat_mix; c++) {
                size_t m = class_mix[c];
                const double *xc, *xp;
                size_t stride_c = V, stride_p = V;
                if (st_child >= 0) { xc = tip_lh + st_child * tip_block + m * nstates; stride_c = 1; }
                else xc = child.lh + base + c * nstates * V + lane;
                if (st_parent >= 0) { xp = tip_lh + st_parent * tip_block + m * nstates; stride_p = 1; }
                else xp = parent.lh + base + c * nstates * V + lane;
                const double *vc = &val[c * nstates], *lc = &lval[c * nstates];
                double a = 0.0, b = 0.0;
                for (size_t i = 0; i < nstates; i++) {
                    double prod = xc[i * stride_c] * xp[i * stride_p];
                    a += prod * vc[i];
                    b += prod * lc[i];
                }
                a_c[c] = a; b_c[c] = b;
                unsigned s = 0;
                if (safe) {
                    if (child.scale) s += child.scale[ptn * ncat_mix + c];
                    if (parent.scale) s += parent.scale[ptn * ncat_mix + c];
                }
                sc[c] = s;
                if (c == 0 || s < min_scale) min_scale = s;
            }
            // class-relative scaling as in the kernel (design 5.3); the common
            // per-pattern factor cancels in dL/L
            double L = 0.0, dL = 0.0;
            for (size_t c = 0; c < ncat_mix; c++) {
                double f = 1.0;
                if (safe) {
                    if (sc[c] == min_scale + 1) f = SCALING_THRESHOLD;
                    else if (sc[c] > min_scale + 1) f = 0.0;
                }
                double w = class_weight[c] * f;
                L += w * a_c[c];
                dL += w * class_rate[c] * b_c[c];
            }
            L = fabs(L) + ptn_invar[ptn];
            if (L <= 0.0 || !std::isfinite(L)) { chunk_bad[ch] = 1; continue; }
            sum += ptn_freq[ptn] * dL / L;
        }
        chunk_sum[ch] = sum;
    }

    double total = 0.0;
    for (int ch = 0; ch < nchunks; ch++) {
        total += chunk_sum[ch];
        if (chunk_bad[ch]) res.valid = false;
    }
    int id = nei->id;
    if (id < 0 || id >= (int)res.dlogl_dt.size()) {
        res.dlogl_dt.resize(max((size_t)id + 1, res.dlogl_dt.size()), 0.0);
    }
    res.dlogl_dt[id] = total;
    res.num_edges++;
}
