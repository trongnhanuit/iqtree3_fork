/*
 * gradientoptimizer.cpp
 *
 * See gradientoptimizer.h and docs/analytical-gradients-design.md.
 */

#include "gradientoptimizer.h"
#include "modelfactory.h"

bool GradientOptimizer::supports(ModelFactory *factory, std::string &why) {
    (void)factory;
    // Stage 0 scaffolding: the pipeline is not available yet, so every model
    // falls back to the default optimiser. Later commits replace this with the
    // real capability check (design doc, section 4).
    why = "analytical-gradient pipeline not yet available";
    return false;
}
