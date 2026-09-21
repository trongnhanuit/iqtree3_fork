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
class ModelParamMap;

class GradientOptimizer {
public:
    /**
     * --ag-selftest: the X kernel's three regimes against a long-double
     * reference (and the naive formula shown to fail for tiny gaps), and the
     * pack/unpack round trip of the parameter map. Prints SELFTEST lines.
     * @return 0 if all pass, 1 otherwise
     */
    static int selfTest(ModelParamMap *map, PhyloTree *tree);

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
     * Branch lengths are compared against the tree's Newton derivative
     * (computeLikelihoodDerv) and central differences; every model parameter
     * of the ModelParamMap vector against central differences in theta space,
     * with the closed-form Q chain rule cross-checked against a
     * finite-difference-of-Q-builder chain rule and two zero-cost identities.
     * @return 0 if every row passed the tolerance, 1 otherwise
     */
    static int gradientCheckOnly(ModelFactory *factory, PhyloTree *tree);
};

#endif /* GRADIENTOPTIMIZER_H */
