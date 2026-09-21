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
 * Stage 1 scope: outside pass, per-branch dlogL/dt, self-check. Later stages
 * add dlogL/dQ (Frechet adjoint), weights, rates and the chain rule.
 */

#ifndef PHYLOGRADIENT_H
#define PHYLOGRADIENT_H

#include <vector>
#include <cstddef>
#include "tree/phylotree.h"

class PhyloGradient {
public:
    struct Result {
        /** tree log-likelihood from the forward pass */
        double logl = 0.0;
        /** dlogL/dt per branch, indexed by Neighbor::id */
        std::vector<double> dlogl_dt;
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
     * Run the forward pass, the outside pass and the per-edge accumulation.
     * Leaves the tree's inside partials untouched, restores current_it and
     * clears theta_computed. Overwrites _pattern_lh / _pattern_lh_cat, so
     * callers needing posteriors must recompute them afterwards.
     * @return res.valid
     */
    bool compute(Result &res);

    /** bytes currently held in the outside-partial pool */
    size_t bufferBytes() const { return buffer_bytes; }

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
    int treeDepth(PhyloNode *node, PhyloNode *dad) const;
    void ensurePool(int depth);
    double *poolLh(int depth) { return pool_lh[depth]; }
    UBYTE *poolScale(int depth) { return pool_scale[depth]; }

    /** visit all edges below node (away from dad); `depth` indexes the pool */
    void visit(PhyloNode *node, PhyloNode *dad, int depth, Result &res);

    /** accumulate one edge given both sides and the neighbour carrying its length/id */
    void accumulateEdge(PhyloNeighbor *nei, const Side &child, const Side &parent, Result &res);

    Side sideOf(PhyloNeighbor *nei_to_side, PhyloNode *side_node) const;

    PhyloTree *tree;
    // cached dimensions (valid during compute())
    size_t nstates = 0, ncat = 0, nmix = 0, ncat_mix = 0, denom = 1, block = 0, tip_block = 0;
    size_t vsize = 1, orig_nptn = 0, nptn = 0;
    bool safe = false, fused = false;
    std::vector<double> class_rate, class_weight;   // per class c
    std::vector<size_t> class_eval_offset;            // per class c: offset into eigenvalues
    std::vector<size_t> class_mix;                    // per class c: mixture component m
    std::vector<double> val, lval;                    // per edge scratch: exp(lambda tau), lambda*exp(lambda tau)

    std::vector<double*> pool_lh;
    std::vector<UBYTE*> pool_scale;
    size_t buffer_bytes = 0;
};

#endif /* PHYLOGRADIENT_H */
