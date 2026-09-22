/*
 * modelparammap.h
 *
 * Unconstrained parameter vector for the analytical-gradient optimiser and
 * the chain rule from the engine's natural gradients (dlogL/dQ_m, class
 * sums, rate sums, root term, invariant-site sums) to that vector.
 *
 * Design: docs/analytical-gradients-design.md, section 6.
 *
 * Blocks, in order (each block may be absent):
 *   S  exchangeabilities: log of every free rate variable of a component, in
 *      the component's own variable layout (DNA rate groups, GTR20 entries;
 *      the last entry stays pinned at 1 as in IQ-TREE). One shared block when
 *      exchangeabilities are linked across the mixture (`+Fk` / --link-exchange
 *      with a GTR-type matrix), else one block per estimated component.
 *   F  per component with estimated frequencies: log ratios relative to the
 *      pinned state (the largest frequency at construction); states with a
 *      frequency at or below ZERO_FREQ are excluded (the kernel drops them).
 *   W  mixture weights: K-1 log ratios relative to the last class.
 *   R  free rates: K-1 log proportion ratios (last pinned) and K-1 rate
 *      parameters rho with r_c = e^{rho_c} / ((1-p_inv) sum_c s_c e^{rho_c}),
 *      rho_{K-1} = 0, so the mean rate stays 1 and the rate/branch-length
 *      scale redundancy is removed; gamma: log(alpha); invariant sites:
 *      logit(p_inv).
 *
 * unpack() writes through public setters only (the component's own variable
 * layout for rates, setStateFrequency, setMixtureWeight, setRate/setProp,
 * setGammaShape, setPInvar), then decomposeRateMatrix(), computePtnInvar()
 * and clearAllPartialLH(). The linked-GTR num_params toggling inside
 * ModelMixture::setVariables is never used.
 */

#ifndef MODELPARAMMAP_H
#define MODELPARAMMAP_H

#include <string>
#include <vector>
#include "phylogradient.h"

class ModelFactory;
class PhyloTree;
class ModelSubst;
class ModelMarkov;
class RateHeterogeneity;

class ModelParamMap {
public:
    ModelParamMap(ModelFactory *factory, PhyloTree *tree);

    /** empty when the model is representable; otherwise a one-line reason */
    const std::string &unsupported() const { return unsupported_; }

    int ndim() const { return (int)names_.size(); }
    const std::string &name(int i) const { return names_[i]; }
    const std::vector<double> &lower() const { return lower_; }
    const std::vector<double> &upper() const { return upper_; }

    /** which mixture components need dlogL/dQ from the engine */
    const std::vector<bool> &needQ() const { return need_q_; }
    /** theta indices of the +FO profile entries (F blocks) */
    const std::vector<int> &profileParams() const { return profile_params_; }
    /** number of components with estimated profiles */
    int numProfileClasses() const { return (int)fblocks_.size(); }

    /** current model -> theta */
    void pack(std::vector<double> &theta) const;
    /** theta -> model (then decomposeRateMatrix, computePtnInvar, clearAllPartialLH) */
    void unpack(const std::vector<double> &theta);

    /**
     * Chain rule: engine result -> dlogL/dtheta. Also fills the natural
     * gradients (for the external oracle), see NaturalGradient.
     */
    void gradient(const PhyloGradient::Result &res, std::vector<double> &g);

    /**
     * Natural-parameter gradients. Conventions (what an oracle must
     * differentiate): dR per rate entry with the kernel's normalisation of Q;
     * dpi per raw frequency entry with pi renormalised to sum 1 inside Q and
     * at the root; dw per raw mixture weight without renormalisation;
     * drate/dprop per raw category rate/proportion; dpinv_direct is the
     * derivative through ptn_invar only (rates and proportions held fixed),
     * dpinv the total including the coupled rates/proportions.
     */
    struct NaturalGradient {
        std::vector<std::vector<double>> dR;      // per component: per rate entry (getNumRateEntries)
        std::vector<std::vector<double>> dpi;     // per component: per state
        std::vector<double> dw;                    // per component
        std::vector<double> drate, dprop;          // per rate category
        double dpinv_direct = 0.0, dpinv = 0.0, dalpha = 0.0;
    };
    const NaturalGradient &natural() const { return nat_; }

    /**
     * Same chain rule but with the derivative of the normalised rate matrix
     * taken by central differences of the Q-builder (no likelihood
     * evaluations); the closed form in gradient() must agree with it.
     */
    void gradientFDQ(const PhyloGradient::Result &res, std::vector<double> &g);

    /** the normalised rate matrix of component m as the kernel sees it (S*S, row-major) */
    void buildQ(int m, std::vector<double> &Q) const;

    /** sum_m <dlogL/dQ_m, Q_m>; equals sum_e t_e dlogL/dt_e (zero-cost identity) */
    double qIdentity(const PhyloGradient::Result &res) const;

private:
    struct SBlock { std::vector<int> comps; int nvar = 0; int offset = 0;
                    std::vector<std::vector<int>> entries; /* per variable: rate entries it drives */ };
    struct FBlock { int comp = 0; int pinned = 0; int offset = 0; std::vector<int> states; /* free states */ };

    ModelMarkov *component(int m) const;
    void readRateVariables(ModelMarkov *comp, std::vector<double> &v) const;
    void writeRateVariables(ModelMarkov *comp, const std::vector<double> &v) const;
    void discoverGroups(ModelMarkov *comp, SBlock &blk) const;
    void buildQFrom(ModelMarkov *comp, const double *entries, const double *freq_raw, std::vector<double> &Q) const;
    void naturalFromQ(const PhyloGradient::Result &res, bool fd_q);
    void naturalToTheta(std::vector<double> &g) const;
    void rateModelJacobian();
    void addParam(const std::string &name, double lo, double hi);

    ModelFactory *factory_;
    PhyloTree *tree_;
    ModelSubst *model_;
    RateHeterogeneity *rate_;
    int nstates_ = 0, nmix_ = 0, ncat_ = 0;
    bool fused_ = false;

    std::vector<SBlock> sblocks_;
    std::vector<FBlock> fblocks_;
    int w_offset_ = -1, r_prop_offset_ = -1, r_rate_offset_ = -1, alpha_offset_ = -1, pinv_offset_ = -1;
    bool rates_free_ = false, props_free_ = false, alpha_free_ = false, pinv_free_ = false;
    std::vector<std::string> names_;
    std::vector<double> lower_, upper_;
    std::vector<bool> need_q_;
    std::string unsupported_;
    // entries held at a floor by the last unpack() (design 6): their theta
    // gradient is reported as zero, the direction being flat there
    std::vector<char> clamped_;
    std::vector<int> profile_params_;
    double freq_floor_ = 1e-4, weight_floor_ = 1e-3;

    NaturalGradient nat_;
    std::vector<double> drate_dalpha_;   // d rate_c / d alpha, by central differences of the gamma quantiles
};

#endif /* MODELPARAMMAP_H */
