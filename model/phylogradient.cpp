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
#include <stdexcept>
#include <cstdlib>
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
    // test hook (design 12): throw at the k-th edge of the outside pass to
    // prove that the RAII guard leaves the tree reusable
    if (const char *e = getenv("AG_TEST_FAULT_EDGE")) fault_edge_ = atoi(e);
}

PhyloGradient::~PhyloGradient() {
    for (auto p : pool_lh) aligned_free(p);
    for (auto p : pool_scale) aligned_free(p);
}

void PhyloGradient::setNeedQ(const std::vector<bool> &need) {
    need_q = need;
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
    p_invar = site_rate->getPInvar();
    // padded per-component strides of the shared eigen buffers (ModelMixture::initMem)
    eval_stride = get_safe_upper_limit(nstates);
    evec_stride = eval_stride * nstates;

    class_rate.assign(ncat_mix, 0.0);
    class_weight.assign(ncat_mix, 0.0);
    class_prop.assign(ncat_mix, 0.0);
    class_mix.assign(ncat_mix, 0);
    class_cat.assign(ncat_mix, 0);
    for (size_t c = 0; c < ncat_mix; c++) {
        size_t mycat = c % ncat;
        size_t m = c / denom;
        class_rate[c] = site_rate->getRate(mycat);
        class_prop[c] = site_rate->getProp(mycat);
        class_weight[c] = class_prop[c] * model->getMixtureWeight(m);
        class_mix[c] = m;
        class_cat[c] = mycat;
    }
    if (need_q.size() != nmix) need_q.assign(nmix, true);
    val.assign(ncat_mix * nstates, 0.0);
    lval.assign(ncat_mix * nstates, 0.0);
    nchunks = max(1, tree->num_packets);
    bool any_q = false;
    for (bool b : need_q) any_q |= b;
    size_t ss = nstates * nstates;
    xker.assign(any_q ? ncat_mix * ss : 0, 0.0);
    A.assign(any_q ? (size_t)nchunks * ncat_mix * ss : 0, 0.0);
    cacheInvariantStates();
}

// Mirror of PhyloTree::computePtnInvar: which states an invariant pattern's
// term p_inv * sum_x pi_bar_x runs over (DNA IUPAC and protein B/Z/J ambiguity).
void PhyloGradient::cacheInvariantStates() {
    inv_mask.assign(orig_nptn, 0);
    inv_gaponly.assign(orig_nptn, 0);
    if (p_invar <= 0.0) return;
    const Alignment *aln = tree->aln;
    const int unknown = aln->STATE_UNKNOWN;
    const int ambi_aa[] = { 4 + 8, 32 + 64, 512 + 1024 };   // B = N|D, Z = Q|E, J = I|L
    for (size_t ptn = 0; ptn < orig_nptn; ptn++) {
        int cstate = (unsigned char)aln->at(ptn).const_char;
        if (cstate > unknown) continue;                 // variable pattern
        if (cstate == unknown) { inv_gaponly[ptn] = 1; continue; }
        uint64_t mask = 0;
        if (cstate < (int)nstates) {
            mask = (uint64_t)1 << cstate;
        } else if (aln->seq_type == SEQ_DNA) {
            int astate = cstate - (int)nstates + 1;
            for (size_t x = 0; x < nstates; x++)
                if (astate & (1 << x)) mask |= (uint64_t)1 << x;
        } else if (aln->seq_type == SEQ_PROTEIN) {
            int astate = cstate - (int)nstates;
            if (astate >= 0 && astate <= 2)
                for (int x = 0; x < 11; x++)
                    if (ambi_aa[astate] & (1 << x)) mask |= (uint64_t)1 << x;
        }
        inv_mask[ptn] = mask;
    }
}

/* ---------------------------------------------------------------------- */
/* X kernel: (e^{a tau} - e^{b tau}) / (a - b), tau e^{a tau} on the diagonal */
/* ---------------------------------------------------------------------- */

double PhyloGradient::xKernel(double lam_j, double lam_k, double tau) {
    double d = lam_j - lam_k;
    // expm1(tau d)/d is accurate for ANY non-zero d, so the diagonal formula is
    // only needed when the gap is exactly zero (or would underflow); a wider
    // window here would cost relative accuracy tau*|d|/2 (design doc 5.3)
    if (fabs(d) < 1e-200)
        return tau * exp(lam_j * tau);
    if (fabs(d * tau) < 0.1)
        return exp(lam_k * tau) * expm1(tau * d) / d;    // algebraically exact, no cancellation
    return (exp(lam_j * tau) - exp(lam_k * tau)) / d;
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
    size_t ss = nstates * nstates;
    res.dlogl_dt.assign(max(tree->branchNum, 1), 0.0);
    res.dlogl_dQ.assign(nmix, vector<double>());
    res.G.assign(nmix, vector<double>());
    for (size_t m = 0; m < nmix; m++)
        if (need_q[m]) { res.G[m].assign(ss, 0.0); res.dlogl_dQ[m].assign(ss, 0.0); }
    res.dlogl_dclass.assign(ncat_mix, 0.0);
    res.dlogl_drate.assign(ncat, 0.0);
    res.root_term.assign(nmix, vector<double>(nstates, 0.0));
    res.inv_state_sum.assign(nstates, 0.0);
    res.inv_total = 0.0;

    PhyloNode *v = (PhyloNode*)saved_it->node;         // dad_branch side
    PhyloNode *u = (PhyloNode*)saved_back->node;       // dad side
    // NOTE(design 5.3): every edge is accumulated with the outside partial on
    // the side containing u, so dlogL/dQ is the derivative of the likelihood
    // rooted at u (this matters only for non-reversible perturbations of Q)
    res.root_dad = u;
    res.root_node = v;

    // root edge (u,v): both directions are resident; u is the "outside" (row) side
    {
        Side side_v = sideOf(saved_it, v);
        Side side_u = sideOf(saved_back, u);
        accumulateEdge(saved_it, side_v, side_u, true, res);
    }

    // buffer pool sized from the actual recursion depth (a caterpillar needs N-2)
    int depth = max(u->isLeaf() ? 0 : treeDepth(u, v), v->isLeaf() ? 0 : treeDepth(v, u));
    ensurePool(depth);

    try {
        if (!u->isLeaf()) visit(u, v, 0, res);
        if (!v->isLeaf()) visit(v, u, 0, res);
    } catch (...) {
        // the guards have restored every neighbour on the way out; restore the
        // tree state too, then let the caller decide (design 12)
        tree->current_it = saved_it;
        tree->current_it_back = saved_back;
        tree->theta_computed = false;
        throw;
    }

    finishQ(res);

    // restore tree state the pass may have disturbed
    tree->current_it = saved_it;
    tree->current_it_back = saved_back;
    tree->theta_computed = false;   // NOTE(design 5.2): kernel sets it per edge; stale theta_all must not be reused

    for (double g : res.dlogl_dt) if (!std::isfinite(g)) res.valid = false;
    for (auto &d : res.dlogl_dQ) for (double g : d) if (!std::isfinite(g)) res.valid = false;
    if (!std::isfinite(res.logl)) res.valid = false;
    return res.valid;
}

void PhyloGradient::visit(PhyloNode *node, PhyloNode *dad, int depth, Result &res) {
    FOR_NEIGHBOR_DECLARE(node, dad, it) {
        PhyloNeighbor *nei = (PhyloNeighbor*)(*it);          // node -> child (inside partial of child)
        PhyloNode *child = (PhyloNode*)nei->node;
        PhyloNeighbor *rev = (PhyloNeighbor*)child->findNeighbor(node);   // child -> node (outside partial)

        OutsideGuard guard(rev, poolLh(depth), poolScale(depth));
        if (fault_edge_ >= 0 && res.num_edges == fault_edge_) {
            fault_edge_ = -1;   // one-shot: later gradients must be analytic again
            throw std::runtime_error("AG_TEST_FAULT_EDGE: injected fault inside the outside pass");
        }
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
        accumulateEdge(nei, side_child, side_parent, false, res);

        if (!child->isLeaf())
            visit(child, node, depth + 1, res);
        // guard restores rev on scope exit
    }
}

/* ---------------------------------------------------------------------- */
/* per-edge accumulation                                                   */
/* ---------------------------------------------------------------------- */

void PhyloGradient::accumulateEdge(PhyloNeighbor *nei, const Side &child, const Side &parent, bool root_edge, Result &res) {
    ModelSubst *model = tree->getModel();
    const double *eval = model->getEigenvalues();
    const double *evec = model->getEigenvectors();
    const double *tip_lh = tree->tip_partial_lh;
    const double *ptn_freq = tree->ptn_freq;
    const double *ptn_invar = tree->ptn_invar;
    const size_t V = vsize, S = nstates, SS = S * S;
    const bool any_q = !A.empty();

    // per-class exp(lambda * rate * t), lambda * exp(...), and X kernels
    for (size_t c = 0; c < ncat_mix; c++) {
        double tau = class_rate[c] * nei->getLength(class_cat[c]);
        const double *ev = eval + class_mix[c] * eval_stride;
        double *vc = &val[c * S], *lc = &lval[c * S];
        for (size_t i = 0; i < S; i++) {
            vc[i] = exp(ev[i] * tau);
            lc[i] = ev[i] * vc[i];
        }
        if (any_q && need_q[class_mix[c]]) {
            double *X = &xker[c * SS];
            for (size_t j = 0; j < S; j++)
                for (size_t k = 0; k < S; k++)
                    X[j * S + k] = xKernel(ev[j], ev[k], tau);
        }
    }

    size_t nblocks = (orig_nptn + V - 1) / V;
    vector<double> dt_chunk(nchunks, 0.0);
    vector<double> drate_chunk((size_t)nchunks * ncat, 0.0);
    vector<double> dclass_chunk(root_edge ? (size_t)nchunks * ncat_mix : 0, 0.0);
    vector<double> inv_state_chunk(root_edge ? (size_t)nchunks * S : 0, 0.0);
    vector<double> inv_total_chunk(root_edge ? nchunks : 0, 0.0);
    vector<double> root_chunk(root_edge ? (size_t)nchunks * nmix * S : 0, 0.0);
    vector<char> chunk_bad(nchunks, 0);
    if (any_q) fill(A.begin(), A.end(), 0.0);
    const Alignment *aln = tree->aln;
    const double t_edge = nei->length;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(tree->num_threads)
#endif
    for (int ch = 0; ch < nchunks; ch++) {
        size_t b_lo = (nblocks * ch) / nchunks, b_hi = (nblocks * (ch + 1)) / nchunks;
        double sum_dt = 0.0;
        double *drate = &drate_chunk[(size_t)ch * ncat];
        double *dclass = root_edge ? &dclass_chunk[(size_t)ch * ncat_mix] : nullptr;
        double *inv_state = root_edge ? &inv_state_chunk[(size_t)ch * S] : nullptr;
        double *rootc = root_edge ? &root_chunk[(size_t)ch * nmix * S] : nullptr;
        double *Ach = any_q ? &A[(size_t)ch * ncat_mix * SS] : nullptr;
        vector<double> a_c(ncat_mix), b_c(ncat_mix), wf(ncat_mix);
        vector<const double*> xc_ptr(ncat_mix), xp_ptr(ncat_mix);
        vector<size_t> sc_ptr(ncat_mix), sp_ptr(ncat_mix);
        vector<unsigned> sc(ncat_mix);
        vector<double> o_state(S), q_state(S);
        for (size_t ptn = b_lo * V; ptn < min(b_hi * V, orig_nptn); ptn++) {
            size_t base = (ptn / V) * V * block, lane = ptn % V;
            int st_child = child.leaf_id >= 0 ? (*aln)[ptn][child.leaf_id] : -1;
            int st_parent = parent.leaf_id >= 0 ? (*aln)[ptn][parent.leaf_id] : -1;
            unsigned min_scale = 0;
            for (size_t c = 0; c < ncat_mix; c++) {
                size_t m = class_mix[c];
                const double *xc, *xp;
                size_t stride_c = V, stride_p = V;
                if (st_child >= 0) { xc = tip_lh + st_child * tip_block + m * S; stride_c = 1; }
                else xc = child.lh + base + c * S * V + lane;
                if (st_parent >= 0) { xp = tip_lh + st_parent * tip_block + m * S; stride_p = 1; }
                else xp = parent.lh + base + c * S * V + lane;
                xc_ptr[c] = xc; xp_ptr[c] = xp; sc_ptr[c] = stride_c; sp_ptr[c] = stride_p;
                const double *vc = &val[c * S], *lc = &lval[c * S];
                double a = 0.0, b = 0.0;
                for (size_t i = 0; i < S; i++) {
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
                wf[c] = class_weight[c] * f;
                L += wf[c] * a_c[c];
                dL += wf[c] * class_rate[c] * b_c[c];
            }
            L = fabs(L) + ptn_invar[ptn];
            if (L <= 0.0 || !std::isfinite(L)) { chunk_bad[ch] = 1; continue; }
            const double g = ptn_freq[ptn] / L;
            sum_dt += g * dL;
            for (size_t c = 0; c < ncat_mix; c++)
                drate[class_cat[c]] += g * wf[c] * t_edge * b_c[c];

            if (root_edge) {
                // class sums for weights/proportions; invariant-site sums
                for (size_t c = 0; c < ncat_mix; c++)
                    dclass[c] += g * (wf[c] / max(class_weight[c], 1e-300)) * a_c[c];
                if (p_invar > 0.0 && (inv_mask[ptn] || inv_gaponly[ptn])) {
                    inv_total_chunk[ch] += g * ptn_invar[ptn] / p_invar;
                    for (size_t x = 0; x < S; x++)
                        if (inv_mask[ptn] & ((uint64_t)1 << x)) inv_state[x] += g * p_invar;
                }
                // root-frequency term: (U o~)_k (U (val o p~))_k per class, state space
                for (size_t c = 0; c < ncat_mix; c++) {
                    if (wf[c] == 0.0) continue;
                    size_t m = class_mix[c];
                    const double *U = evec + m * evec_stride;
                    const double *xc = xc_ptr[c], *xp = xp_ptr[c];
                    const double *vc = &val[c * S];
                    for (size_t x = 0; x < S; x++) {
                        double o = 0.0, q = 0.0;
                        const double *Ux = U + x * S;
                        for (size_t i = 0; i < S; i++) {
                            o += Ux[i] * xp[i * sp_ptr[c]];
                            q += Ux[i] * vc[i] * xc[i * sc_ptr[c]];
                        }
                        rootc[m * S + x] += g * wf[c] * o * q;
                    }
                }
            }
            // rank-one accumulation for dlogL/dQ (row = outside/parent side, column = child side)
            if (any_q) {
                for (size_t c = 0; c < ncat_mix; c++) {
                    if (!need_q[class_mix[c]] || wf[c] == 0.0) continue;
                    double coef = g * wf[c];
                    const double *xc = xc_ptr[c], *xp = xp_ptr[c];
                    double *Acc = Ach + c * SS;
                    for (size_t j = 0; j < S; j++) {
                        double oj = coef * xp[j * sp_ptr[c]];
                        if (oj == 0.0) continue;
                        double *row = Acc + j * S;
                        for (size_t k = 0; k < S; k++)
                            row[k] += oj * xc[k * sc_ptr[c]];
                    }
                }
            }
        }
        dt_chunk[ch] = sum_dt;
    }

    // fixed-order reductions
    double total = 0.0;
    for (int ch = 0; ch < nchunks; ch++) {
        total += dt_chunk[ch];
        if (chunk_bad[ch]) res.valid = false;
        for (size_t r = 0; r < ncat; r++) res.dlogl_drate[r] += drate_chunk[(size_t)ch * ncat + r];
        if (root_edge) {
            for (size_t c = 0; c < ncat_mix; c++) res.dlogl_dclass[c] += dclass_chunk[(size_t)ch * ncat_mix + c];
            for (size_t x = 0; x < S; x++) res.inv_state_sum[x] += inv_state_chunk[(size_t)ch * S + x];
            res.inv_total += inv_total_chunk[ch];
            for (size_t m = 0; m < nmix; m++)
                for (size_t x = 0; x < S; x++)
                    res.root_term[m][x] += root_chunk[((size_t)ch * nmix + m) * S + x];
        }
    }
    if (any_q) {
        vector<double> Ac(SS);
        for (size_t c = 0; c < ncat_mix; c++) {
            size_t m = class_mix[c];
            if (!need_q[m]) continue;
            fill(Ac.begin(), Ac.end(), 0.0);
            for (int ch = 0; ch < nchunks; ch++) {
                const double *src = &A[((size_t)ch * ncat_mix + c) * SS];
                for (size_t e = 0; e < SS; e++) Ac[e] += src[e];
            }
            const double *X = &xker[c * SS];
            double *G = res.G[m].data();
            for (size_t e = 0; e < SS; e++) G[e] += X[e] * Ac[e];
        }
    }

    int id = nei->id;
    if (id < 0 || id >= (int)res.dlogl_dt.size())
        res.dlogl_dt.resize(max((size_t)id + 1, res.dlogl_dt.size()), 0.0);
    res.dlogl_dt[id] = total;
    res.num_edges++;
}

// D_m = U^-T G_m U^T with U = eigenvectors (row = state, column = eigen index)
// and U^-1 stored transposed as inv_eigenvectors_transposed (row = state).
void PhyloGradient::finishQ(Result &res) {
    ModelSubst *model = tree->getModel();
    const double *evec = model->getEigenvectors();
    const double *inv_t = model->getInverseEigenvectorsTransposed();
    const size_t S = nstates;
    vector<double> T(S * S);
    for (size_t m = 0; m < nmix; m++) {
        if (!need_q[m]) continue;
        const double *U = evec + m * evec_stride;        // U[b][k]
        const double *W = inv_t + m * evec_stride;       // W[a][j] = (U^-1)[j][a]
        const double *G = res.G[m].data();               // G[j][k]
        double *D = res.dlogl_dQ[m].data();
        // T[j][b] = sum_k G[j][k] U[b][k]
        for (size_t j = 0; j < S; j++)
            for (size_t b = 0; b < S; b++) {
                double s = 0.0;
                for (size_t k = 0; k < S; k++) s += G[j * S + k] * U[b * S + k];
                T[j * S + b] = s;
            }
        // D[a][b] = sum_j W[a][j] T[j][b]
        for (size_t a = 0; a < S; a++)
            for (size_t b = 0; b < S; b++) {
                double s = 0.0;
                for (size_t j = 0; j < S; j++) s += W[a * S + j] * T[j * S + b];
                D[a * S + b] = s;
            }
    }
}
