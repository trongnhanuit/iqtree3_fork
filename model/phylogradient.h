/*
 * phylogradient.h
 *
 * Exact reverse-mode gradients of the tree log-likelihood for reversible
 * models and mixtures. Used by --analytical-gradients (model/gradientoptimizer).
 *
 * Design: docs/analytical-gradients-design.md, sections 5 and 6.
 *
 * Idea: after tree->computeLikelihood() every internal node holds one
 * "inside" partial oriented toward the current root branch. The engine walks
 * the tree from that branch and, for every edge (n -> c), computes the
 * "outside" partial of c (the likelihood of everything on n's side) by
 * attaching a private buffer to the reverse neighbour c->findNeighbor(n) and
 * calling the tree's own likelihood kernel on it. The kernel then computes
 * exactly that partial with its own scaling rules and SIMD code, and returns
 * the tree log-likelihood evaluated on that edge as a free self-check.
 *
 * Invariant (load-bearing): the buffer is attached BEFORE the traversal runs.
 * PhyloTree::reorientPartialLh() returns early when partial_lh is set, so the
 * kernel never steals an inside slot (MemSlotVector::takeover moves pointers
 * unconditionally otherwise). The RAII guard restores the neighbour on every
 * exit path, including exceptions.
 *
 * Per edge e and class c = (mixture component m, rate category r), with
 * tau = rate_r * t_e, outside partial o~ (eigen coordinates, without pi
 * because U^-1 = U^T Pi) and inside partial p~:
 *   L_ptn      = sum_c w_c f_c sum_i o~_i e^{lambda_i tau} p~_i + ptn_invar
 *   dL/dt_e    = sum_c w_c f_c rate_r sum_i o~_i lambda_i e^{lambda_i tau} p~_i
 *   dL/dQ_m    = U^-T [ sum_e sum_ptn (freq/L) w f X(Lambda_m, tau) o (o~ p~^T) ] U^T
 * with X_jk = (e^{lambda_j tau} - e^{lambda_k tau}) / (lambda_j - lambda_k)
 * (diagonal tau e^{lambda_j tau}; three numerical regimes) and f_c the
 * kernel's per-class safe-numeric rescaling. Class sums, the root-frequency
 * term and invariant-site sums are taken at the root edge.
 */

#ifndef PHYLOGRADIENT_H
#define PHYLOGRADIENT_H

#include <vector>
#include <cstddef>
#include <cstdint>
#include "tree/phylotree.h"

class PhyloGradient {
public:
    struct Result {
        /** tree log-likelihood from the forward pass */
        double logl = 0.0;
        /** dlogL/dt per branch, indexed by Neighbor::id */
        std::vector<double> dlogl_dt;
        /** per mixture component m: dlogL/dQ_m, S*S row-major, all entries independent */
        std::vector<std::vector<double>> dlogl_dQ;
        /** per mixture component m: accumulated G_m in the eigen basis (for tests) */
        std::vector<std::vector<double>> G;
        /** per class c: dlogL/d(omega_c), omega_c = prop_r * w_m, taken at the root edge */
        std::vector<double> dlogl_dclass;
        /** per rate category r: dlogL/d(rate_r) summed over edges and components */
        std::vector<double> dlogl_drate;
        /** per component m, per state k: root-frequency term sum_ptn (freq/L) sum_r w f (U o~)_k (U e^{L tau} p~)_k */
        std::vector<std::vector<double>> root_term;
        /** per state x: sum over constant patterns containing x of (freq/L) * p_inv */
        std::vector<double> inv_state_sum;
        /** sum over patterns with an invariant term of (freq/L) * ptn_invar / p_inv (0 if p_inv == 0) */
        double inv_total = 0.0;
        /** endpoints of the root branch: dlogL/dQ is the derivative of the likelihood rooted at root_dad */
        PhyloNode *root_dad = nullptr, *root_node = nullptr;
        /** number of edges visited (== number of branches) */
        int num_edges = 0;
        /** self-check: max |lnL(edge) - lnL(root)| over all edges the kernel evaluated */
        double max_edge_logl_diff = 0.0;
        /** false if any non-finite value was produced */
        bool valid = true;
    };

    explicit PhyloGradient(PhyloTree *tree);
    ~PhyloGradient();
    PhyloGradient(const PhyloGradient&) = delete;
    PhyloGradient& operator=(const PhyloGradient&) = delete;

    /**
     * Choose which mixture components need dlogL/dQ (the S^2 rank-one
     * accumulation is skipped for the others, e.g. fixed C10 profiles).
     * Default: all components.
     */
    void setNeedQ(const std::vector<bool> &need);

    /**
     * Run the forward pass, the outside pass and the per-edge accumulation.
     * Leaves the tree's inside partials untouched, restores current_it and
     * clears theta_computed. Overwrites _pattern_lh / _pattern_lh_cat, so
     * callers needing posteriors must recompute them afterwards.
     * @return res.valid
     */
    bool compute(Result &res);

    /** bytes currently held in the outside-partial pool */
    size_t bufferBytes() const { return buffer_bytes; }

    /** the divided-difference kernel X(Lambda, tau) (public for --ag-selftest) */
    static double xKernel(double lam_j, double lam_k, double tau);

private:
    /** one side of an edge: an eigen-space partial (with scale counts) or a leaf */
    struct Side {
        const double *lh = nullptr;   // partial in eigen coordinates, kernel layout
        const UBYTE *scale = nullptr; // scale counts (per pattern, or per pattern x class when safe)
        int leaf_id = -1;             // >= 0: use the tip vectors of this leaf instead
    };

    /** RAII attach/restore of a private buffer on a reverse neighbour */
    class OutsideGuard {
    public:
        OutsideGuard(PhyloNeighbor *nei, double *buf, UBYTE *scale);
        ~OutsideGuard() noexcept;
        OutsideGuard(const OutsideGuard&) = delete;
        OutsideGuard& operator=(const OutsideGuard&) = delete;
    private:
        PhyloNeighbor *nei;
        double *saved_lh;
        UBYTE *saved_scale;
        int saved_computed;
    };

    void cacheDimensions();
    void cacheInvariantStates();
    int treeDepth(PhyloNode *node, PhyloNode *dad) const;
    void ensurePool(int depth);
    double *poolLh(int depth) { return pool_lh[depth]; }
    UBYTE *poolScale(int depth) { return pool_scale[depth]; }

    /** visit all edges below node (away from dad); `depth` indexes the pool */
    void visit(PhyloNode *node, PhyloNode *dad, int depth, Result &res);

    /** accumulate one edge given both sides and the neighbour carrying its length/id */
    void accumulateEdge(PhyloNeighbor *nei, const Side &child, const Side &parent, bool root_edge, Result &res);

    /** finish: G_m -> dlogL/dQ_m in the state basis */
    void finishQ(Result &res);

    Side sideOf(PhyloNeighbor *nei_to_side, PhyloNode *side_node) const;

    PhyloTree *tree;
    // cached dimensions (valid during compute())
    size_t nstates = 0, ncat = 0, nmix = 0, ncat_mix = 0, denom = 1, block = 0, tip_block = 0;
    size_t vsize = 1, orig_nptn = 0, nptn = 0, eval_stride = 0, evec_stride = 0;
    bool safe = false, fused = false;
    double p_invar = 0.0;
    std::vector<double> class_rate, class_weight, class_prop;   // per class c
    std::vector<size_t> class_mix, class_cat;                    // per class c
    std::vector<bool> need_q;                                    // per component m
    std::vector<double> val, lval;                               // per edge scratch: exp(lambda tau), lambda*exp(lambda tau)
    std::vector<double> xker;                                    // per edge scratch: X per class, S*S each
    std::vector<double> A;                                       // per edge scratch: per chunk, per class, S*S
    std::vector<uint64_t> inv_mask;                              // per pattern: states of an invariant pattern (0 = none)
    std::vector<char> inv_gaponly;                               // per pattern: gap-only invariant pattern
    int nchunks = 1;

    int fault_edge_ = -1;                                        // AG_TEST_FAULT_EDGE (tests only)
    std::vector<double*> pool_lh;
    std::vector<UBYTE*> pool_scale;
    size_t buffer_bytes = 0;
};

#endif /* PHYLOGRADIENT_H */
