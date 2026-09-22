/*
 * gradientoptimizer.h
 *
 * Analytical-gradient model parameter optimisation for IQ-TREE 3: a joint
 * BFGS polish of all substitution-model parameters driven by exact
 * reverse-mode gradients (model/phylogradient.h, model/modelparammap.h),
 * alternated with the existing Newton branch-length optimisation. Engaged
 * only when --analytical-gradients is set and supports() returns true;
 * otherwise ModelFactory::optimizeParameters runs its unchanged alternating
 * loop.
 *
 * Design: docs/analytical-gradients-design.md
 */

#ifndef GRADIENTOPTIMIZER_H
#define GRADIENTOPTIMIZER_H

#include <string>
#include <vector>
#include <memory>
#include "utils/optimization.h"
#include "utils/tools.h"

class ModelFactory;
class PhyloTree;
class ModelParamMap;
class PhyloGradient;

class GradientOptimizer : public Optimization {
public:
    explicit GradientOptimizer(ModelFactory *factory);
    virtual ~GradientOptimizer();

    /**
     * Decide whether the analytic pipeline handles this model/tree configuration.
     * Has no lasting side effects. When it returns false, `why` holds a one-line
     * reason that the caller prints once as a NOTE before falling back to the
     * default optimiser.
     */
    static bool supports(ModelFactory *factory, std::string &why);

    /**
     * Replacement for the alternating loop of ModelFactory::optimizeParameters
     * (design 3.2): per round a branch-length step according to fixed_len,
     * then one joint analytic-gradient BFGS polish of all model parameters;
     * rounds continue while a round improves the log-likelihood by more than
     * logl_epsilon. Ends with the same terminal branch optimisation as the
     * default loop, then restores the best state seen if the final
     * log-likelihood is below both the entry score and the best round.
     * @param entry_logl the log-likelihood at entry (the prologue's cur_lh)
     * @return the final log-likelihood
     */
    double optimize(int fixed_len, bool write_info, double logl_epsilon, double gradient_epsilon, double entry_logl);

    /** number of rounds the last optimize() ran (for the "took N rounds" message) */
    int rounds() const { return rounds_; }

    /**
     * --ag-selftest: the X kernel's three regimes against a long-double
     * reference (and the naive formula shown to fail for tiny gaps), and the
     * pack/unpack round trip of the parameter map. Prints SELFTEST lines.
     * @return 0 if all pass, 1 otherwise
     */
    static int selfTest(ModelParamMap *map, PhyloTree *tree);

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

    // ---- Optimization interface (minimised function is -logL over theta) ----
    virtual int getNDim();
    virtual double targetFunk(double x[]);
    virtual double derivativeFunk(double x[], double dfx[]);
    /** never randomise at a bound: the bounds are numerical fences only (design 6.1) */
    virtual bool restartParameters(double guess[], int ndim, double lower[], double upper[], bool bound_check[], int iteration) { return false; }
    /** improvement-based stopping rule (design 6.4) */
    virtual bool stopEarly(int iter, double f_prev, double f_new);

private:
    struct BestState { std::vector<double> theta; DoubleVector brlen; double logl = -1e300; bool valid = false; };

    void snapshot(BestState &st, double logl);
    void restore(const BestState &st);
    /** joint BFGS (or L-BFGS-B) over theta from the current model state; returns logL */
    double polish(double gradient_epsilon);
    /** identical +FO profiles are a fixed point of every gradient method (design 6.4); make them distinct */
    void breakSymmetry(bool write_info);
    /** --ag-start cold: profiles jittered around the empirical frequencies, equal weights (design 11) */
    void coldStart(bool write_info);
    /** --ag-multistart: score jittered candidates, refine the best few with EM, keep the best (design 11) */
    void multiStart(bool write_info, double gradient_epsilon);
    /** per component ||U U^-1 - I||_inf < 1e-8 (design 5.5) */
    bool eigenResidualOk() const;
    /** append gradient-check rows for the current theta/gradient to <prefix>.gradcheck.tsv */
    void gradientCheckStep(const std::vector<double> &theta, const std::vector<double> &g);
    void say(const std::string &line) const;
    /** --ag-abort-after <phase>: write the checkpoint (never inside an OpenMP region) and exit, for resume tests */
    void abortAfter(const std::string &phase);

    // ---- EM axes (design 10): each is one M-step, accepted only if the
    //      production log-likelihood does not decrease; returns the new logL ----
    /** W: mixture weights by ModelMixture::optimizeWeights */
    double emWeights(double cur_lh);
    /** R: free rates/proportions by RateFree::optimizeWithEM, alpha/p_inv by the rate model's own optimiser */
    double emRates(double cur_lh, double gradient_epsilon);
    /** F: profiles from posterior-weighted site compositions */
    double emProfiles(double cur_lh);
    /** restore `before` if `after_lh` is below `before_lh`; returns the accepted logL */
    double acceptOrRevert(const BestState &before, double before_lh, double after_lh, const char *axis);
    /** put the model back on the parametrised manifold (floors, mean rate 1, ptn_invar) and return logL */
    double renormalise();

    ModelFactory *factory_;
    PhyloTree *tree_;
    std::unique_ptr<ModelParamMap> map_;
    std::unique_ptr<PhyloGradient> engine_;
    int rounds_ = 0;
    double logl_epsilon_ = 0.001;
    // stopping rule state for one polish
    int stop_iter_ = 0;
    double delta_max_ = 0.0;
    // counters (--ag-stats)
    long n_lh_ = 0, n_grad_ = 0, n_fd_fallback_ = 0, n_bfgs_iter_ = 0;
    long n_em_w_ = 0, n_em_r_ = 0, n_em_f_ = 0, n_em_revert_ = 0;
    bool em_enabled_ = false;
    long n_multistart_ = 0;
    std::string em_axes_;
    bool warned_fallback_ = false;
    bool gradcheck_started_ = false;
    std::vector<double> theta_, grad_;
};

#endif /* GRADIENTOPTIMIZER_H */
