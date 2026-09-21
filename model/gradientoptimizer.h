/*
 * gradientoptimizer.h
 *
 * Analytical-gradient model parameter optimisation for IQ-TREE 3:
 * per-axis EM warm start followed by a joint BFGS polish driven by exact
 * reverse-mode gradients (tree/phylogradient.h). Engaged only when
 * --analytical-gradients is set and supports() returns true; otherwise
 * ModelFactory::optimizeParameters runs its unchanged alternating loop.
 *
 * Design: docs/analytical-gradients-design.md
 */

#ifndef GRADIENTOPTIMIZER_H
#define GRADIENTOPTIMIZER_H

#include <string>

class ModelFactory;
class PhyloTree;

class GradientOptimizer {
public:
    /**
     * Decide whether the analytic pipeline handles this model/tree configuration.
     * Has no side effects. When it returns false, `why` holds a one-line reason
     * that the caller prints once as a NOTE before falling back to the default
     * optimiser.
     */
    static bool supports(ModelFactory *factory, std::string &why);

    /**
     * Debug tool behind --ag-gradient-check-only / --ag-dump-gradient: compute
     * the analytic gradient at the current parameter values, compare it with
     * independent references, write <prefix>.gradcheck.tsv (and
     * <prefix>.aggrad.tsv for the external oracle) and print a summary.
     * Stage 1 covers branch lengths, compared against the tree's Newton
     * derivative (computeLikelihoodDerv) and central differences.
     * @return 0 if every row passed the tolerance, 1 otherwise
     */
    static int gradientCheckOnly(ModelFactory *factory, PhyloTree *tree);
};

#endif /* GRADIENTOPTIMIZER_H */
