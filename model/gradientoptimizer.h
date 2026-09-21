/*
 * gradientoptimizer.h
 *
 * Analytical-gradient model parameter optimisation for IQ-TREE 3:
 * per-axis EM warm start followed by a joint BFGS polish driven by exact
 * reverse-mode gradients. Engaged only when --analytical-gradients is set
 * and supports() returns true; otherwise ModelFactory::optimizeParameters
 * runs its unchanged alternating loop.
 *
 * Design: docs/analytical-gradients-design.md
 */

#ifndef GRADIENTOPTIMIZER_H
#define GRADIENTOPTIMIZER_H

#include <string>

class ModelFactory;

class GradientOptimizer {
public:
    /**
     * Decide whether the analytic pipeline handles this model/tree configuration.
     * Has no side effects. When it returns false, `why` holds a one-line reason
     * that the caller prints once as a NOTE before falling back to the default
     * optimiser.
     */
    static bool supports(ModelFactory *factory, std::string &why);
};

#endif /* GRADIENTOPTIMIZER_H */
