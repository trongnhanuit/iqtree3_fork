/*
 * modelparammap.cpp
 *
 * See modelparammap.h and docs/analytical-gradients-design.md, section 6.
 */

#include "modelparammap.h"
#include "modelfactory.h"
#include "modelmarkov.h"
#include "modelmixture.h"
#include "rateheterogeneity.h"
#include "utils/tools.h"
#include <cmath>
#include <algorithm>

using namespace std;

namespace {
const double LOG_RATE_LO = log(MIN_RATE), LOG_RATE_HI = log(MAX_RATE);
const double LOG_RATIO_LO = -12.0, LOG_RATIO_HI = 12.0;
const double LOG_FREE_RATE_LO = log(1e-3), LOG_FREE_RATE_HI = log(1e3);
const double LOGIT_PINV_LO = -14.0, LOGIT_PINV_HI = 7.0;

inline double logistic(double x) { return 1.0 / (1.0 + exp(-x)); }
inline double logit(double p) { return log(p / (1.0 - p)); }
}

/* ---------------------------------------------------------------------- */
/* construction                                                            */
/* ---------------------------------------------------------------------- */

ModelParamMap::ModelParamMap(ModelFactory *factory, PhyloTree *tree)
    : factory_(factory), tree_(tree), model_(tree->getModel()), rate_(tree->getRate()) {
    nstates_ = model_->num_states;
    nmix_ = model_->getNMixtures();
    ncat_ = rate_->getNRate();
    fused_ = factory->fused_mix_rate;
    need_q_.assign(nmix_, false);

    if (fused_) { unsupported_ = "fused mixture-rate model (*R)"; return; }
    if (model_->isSiteSpecificModel()) { unsupported_ = "site-specific model"; return; }
    if (rate_->isSiteSpecificRate()) { unsupported_ = "site-specific rates"; return; }
    if (rate_->isHeterotachy()) { unsupported_ = "heterotachy model"; return; }
    for (int m = 0; m < nmix_; m++) {
        ModelMarkov *comp = component(m);
        if (!comp) { unsupported_ = "model class without a rate matrix"; return; }
        if (!comp->isReversible()) { unsupported_ = "non-reversible model"; return; }
        if (!comp->half_matrix) { unsupported_ = "full (non-symmetric) rate matrix layout"; return; }
        if (comp->ignore_state_freq) { unsupported_ = "rate matrix ignoring state frequencies"; return; }
        if (comp->getNumRateEntries() != nstates_ * (nstates_ - 1) / 2) { unsupported_ = "unexpected rate entry count"; return; }
    }

    // ---- S blocks --------------------------------------------------------
    bool linked = Params::getInstance().optimize_linked_gtr && model_->isMixture() && nmix_ > 1
                  && component(0)->getNDim() > 0 && component(0)->getNParams() > 0;
    if (linked) {
        // NOTE(design 4): linking is keyed on num_params > 0 of the first class;
        // `+Fk` sets optimize_linked_gtr even for LG-type components (0 params).
        SBlock blk;
        ModelMarkov *c0 = component(0);
        for (int m = 0; m < nmix_; m++) {
            ModelMarkov *c = component(m);
            if (c->getNParams() != c0->getNParams()) { unsupported_ = "linked exchangeabilities across components with different rate layouts"; return; }
            blk.comps.push_back(m);
        }
        blk.nvar = c0->getNParams();
        discoverGroups(c0, blk);
        blk.offset = ndim();
        for (int i = 0; i < blk.nvar; i++)
            addParam("S[*][" + convertIntToString(i) + "]", LOG_RATE_LO, LOG_RATE_HI);
        sblocks_.push_back(blk);
        for (int m = 0; m < nmix_; m++) need_q_[m] = true;
    } else {
        for (int m = 0; m < nmix_; m++) {
            ModelMarkov *c = component(m);
            if (c->getNDim() == 0 || c->getNParams() <= 0) continue;
            SBlock blk;
            blk.comps.push_back(m);
            blk.nvar = c->getNParams();
            discoverGroups(c, blk);
            blk.offset = ndim();
            for (int i = 0; i < blk.nvar; i++)
                addParam("S[" + convertIntToString(m) + "][" + convertIntToString(i) + "]", LOG_RATE_LO, LOG_RATE_HI);
            sblocks_.push_back(blk);
            need_q_[m] = true;
        }
    }
    // a component whose free dimensions are neither rates nor plain frequencies
    // (e.g. the special DNA frequency parametrisations) is not representable
    for (int m = 0; m < nmix_; m++) {
        ModelMarkov *c = component(m);
        if (c->getNDim() == 0) continue;
        int expect = max(0, c->getNParams()) + (c->getFreqType() == FREQ_ESTIMATE ? nstates_ - 1 : 0);
        if (c->getNDim() != expect) { unsupported_ = "unsupported frequency parametrisation"; return; }
    }

    // ---- F blocks --------------------------------------------------------
    for (int m = 0; m < nmix_; m++) {
        ModelMarkov *c = component(m);
        if (c->getNDim() == 0 || c->getFreqType() != FREQ_ESTIMATE) continue;
        FBlock blk;
        blk.comp = m;
        const double *pi = c->state_freq;
        blk.pinned = (int)(max_element(pi, pi + nstates_) - pi);
        blk.offset = ndim();
        for (int k = 0; k < nstates_; k++) {
            if (k == blk.pinned || pi[k] <= ZERO_FREQ) continue;
            blk.states.push_back(k);
            addParam("F[" + convertIntToString(m) + "][" + tree_->aln->convertStateBackStr(k) + "]", LOG_RATIO_LO, LOG_RATIO_HI);
        }
        fblocks_.push_back(blk);
        need_q_[m] = true;
    }

    // ---- W block ---------------------------------------------------------
    if (model_->isMixture() && nmix_ > 1) {
        ModelMixture *mix = dynamic_cast<ModelMixture*>(model_);
        if (mix && !mix->fix_prop && !mix->fixed_parameters) {
            w_offset_ = ndim();
            for (int m = 0; m < nmix_ - 1; m++)
                addParam("W[" + convertIntToString(m) + "]", LOG_RATIO_LO, LOG_RATIO_HI);
        }
    }

    // ---- R block ---------------------------------------------------------
    pinv_free_ = !rate_->isFixPInvar();
    alpha_free_ = rate_->isGammaRate() && !rate_->isFixGammaShape();
    if (rate_->isFreeRate()) {
        int d = rate_->getNDim() - (pinv_free_ ? 1 : 0);
        rates_free_ = d >= ncat_ - 1;
        props_free_ = d >= 2 * ncat_ - 2;
        if (props_free_) {
            r_prop_offset_ = ndim();
            for (int c = 0; c < ncat_ - 1; c++)
                addParam("R.prop[" + convertIntToString(c) + "]", LOG_RATIO_LO, LOG_RATIO_HI);
        }
        if (rates_free_) {
            r_rate_offset_ = ndim();
            for (int c = 0; c < ncat_ - 1; c++)
                addParam("R.rate[" + convertIntToString(c) + "]", LOG_FREE_RATE_LO, LOG_FREE_RATE_HI);
        }
    }
    if (alpha_free_) {
        alpha_offset_ = ndim();
        addParam("alpha", log(max(Params::getInstance().min_gamma_shape, 1e-6)), log(MAX_GAMMA_SHAPE));
    }
    if (pinv_free_) {
        pinv_offset_ = ndim();
        addParam("pinv", LOGIT_PINV_LO, LOGIT_PINV_HI);
    }
}

void ModelParamMap::addParam(const string &name, double lo, double hi) {
    names_.push_back(name);
    lower_.push_back(lo);
    upper_.push_back(hi);
}

ModelMarkov *ModelParamMap::component(int m) const {
    if (model_->isMixture()) {
        ModelMixture *mix = dynamic_cast<ModelMixture*>(model_);
        return mix ? mix->at(m) : nullptr;
    }
    return dynamic_cast<ModelMarkov*>(model_);
}

// The component's own 1-indexed variable vector (rates first, then estimated
// frequencies); only the first nvar entries are used here.
void ModelParamMap::readRateVariables(ModelMarkov *comp, vector<double> &v) const {
    v.assign(comp->getNDim() + 1, 0.0);
    comp->setVariables(v.data());
}

void ModelParamMap::writeRateVariables(ModelMarkov *comp, const vector<double> &v) const {
    vector<double> tmp(v);
    comp->getVariables(tmp.data());
}

// Which rate entries each variable drives (DNA rate groups such as HKY's
// "010010", GTR20's one-to-one layout): perturb one variable at a time and
// see which entries change. The mapping is a plain copy, so this is exact.
void ModelParamMap::discoverGroups(ModelMarkov *comp, SBlock &blk) const {
    int nrate = comp->getNumRateEntries();
    vector<double> base, cur(nrate), base_rates(nrate);
    readRateVariables(comp, base);
    comp->getRateMatrix(base_rates.data());
    blk.entries.assign(blk.nvar, vector<int>());
    for (int i = 0; i < blk.nvar; i++) {
        vector<double> v(base);
        v[i + 1] = base[i + 1] * 2.0 + 1.0;
        writeRateVariables(comp, v);
        comp->getRateMatrix(cur.data());
        for (int e = 0; e < nrate; e++)
            if (cur[e] != base_rates[e]) blk.entries[i].push_back(e);
        writeRateVariables(comp, base);
    }
    comp->setRateMatrix(base_rates.data());
}

/* ---------------------------------------------------------------------- */
/* pack / unpack                                                           */
/* ---------------------------------------------------------------------- */

void ModelParamMap::pack(vector<double> &theta) const {
    theta.assign(ndim(), 0.0);
    for (auto &blk : sblocks_) {
        vector<double> v;
        readRateVariables(component(blk.comps[0]), v);
        for (int i = 0; i < blk.nvar; i++) theta[blk.offset + i] = log(max(v[i + 1], 1e-300));
    }
    for (auto &blk : fblocks_) {
        const double *pi = component(blk.comp)->state_freq;
        for (size_t j = 0; j < blk.states.size(); j++)
            theta[blk.offset + j] = log(pi[blk.states[j]] / pi[blk.pinned]);
    }
    if (w_offset_ >= 0) {
        double wl = model_->getMixtureWeight(nmix_ - 1);
        for (int m = 0; m < nmix_ - 1; m++) theta[w_offset_ + m] = log(model_->getMixtureWeight(m) / wl);
    }
    if (props_free_) {
        double pl = rate_->getProp(ncat_ - 1);
        for (int c = 0; c < ncat_ - 1; c++) theta[r_prop_offset_ + c] = log(rate_->getProp(c) / pl);
    }
    if (rates_free_) {
        double rl = rate_->getRate(ncat_ - 1);
        for (int c = 0; c < ncat_ - 1; c++) theta[r_rate_offset_ + c] = log(rate_->getRate(c) / rl);
    }
    if (alpha_free_) theta[alpha_offset_] = log(rate_->getGammaShape());
    if (pinv_free_) theta[pinv_offset_] = logit(min(max(rate_->getPInvar(), 1e-12), 1.0 - 1e-12));
}

void ModelParamMap::unpack(const vector<double> &theta) {
    ASSERT((int)theta.size() == ndim());
    // S: rate variables in the component's own layout; linked blocks copy the
    // first component's entries to the others
    for (auto &blk : sblocks_) {
        ModelMarkov *c0 = component(blk.comps[0]);
        vector<double> v;
        readRateVariables(c0, v);
        for (int i = 0; i < blk.nvar; i++) v[i + 1] = exp(theta[blk.offset + i]);
        writeRateVariables(c0, v);
        if (blk.comps.size() > 1) {
            vector<double> entries(c0->getNumRateEntries());
            c0->getRateMatrix(entries.data());
            for (size_t j = 1; j < blk.comps.size(); j++)
                component(blk.comps[j])->setRateMatrix(entries.data());
        }
    }
    // F: raw entries relative to the pinned state, then normalised to sum 1;
    // excluded (ZERO_FREQ) states keep their current ratio to the pinned state
    for (auto &blk : fblocks_) {
        ModelMarkov *c = component(blk.comp);
        vector<double> pi(c->state_freq, c->state_freq + nstates_);
        double pinned = pi[blk.pinned];
        for (int k = 0; k < nstates_; k++) pi[k] /= pinned;
        for (size_t j = 0; j < blk.states.size(); j++) pi[blk.states[j]] = exp(theta[blk.offset + j]);
        double sum = 0.0;
        for (int k = 0; k < nstates_; k++) sum += pi[k];
        for (int k = 0; k < nstates_; k++) pi[k] /= sum;
        c->setStateFrequency(pi.data());
    }
    // W
    if (w_offset_ >= 0) {
        double sum = 1.0;
        for (int m = 0; m < nmix_ - 1; m++) sum += exp(theta[w_offset_ + m]);
        for (int m = 0; m < nmix_ - 1; m++) model_->setMixtureWeight(m, exp(theta[w_offset_ + m]) / sum);
        model_->setMixtureWeight(nmix_ - 1, 1.0 / sum);
    }
    // R: p_inv first (the free-rate proportions carry the factor 1-p_inv)
    double p = rate_->getPInvar();
    if (pinv_free_) p = logistic(theta[pinv_offset_]);
    if (alpha_free_) rate_->setGammaShape(exp(theta[alpha_offset_]));   // keeps the current rate scale
    if (pinv_free_) rate_->setPInvar(p);                                 // Gamma+I: rescales rates to 1/(1-p)
    if (rate_->isFreeRate() && (rates_free_ || props_free_)) {
        vector<double> s(ncat_), r(ncat_);
        double sum = 0.0;
        for (int c = 0; c < ncat_; c++) {
            s[c] = props_free_ ? (c < ncat_ - 1 ? exp(theta[r_prop_offset_ + c]) : 1.0) : rate_->getProp(c);
            sum += s[c];
        }
        for (int c = 0; c < ncat_; c++) s[c] /= sum;
        double Z = 0.0;
        for (int c = 0; c < ncat_; c++) {
            r[c] = rates_free_ ? (c < ncat_ - 1 ? exp(theta[r_rate_offset_ + c]) : 1.0) : rate_->getRate(c);
            Z += s[c] * r[c];
        }
        for (int c = 0; c < ncat_; c++) {
            r[c] /= (1.0 - p) * Z;      // mean rate sum_c prop_c r_c = 1
            rate_->setProp(c, (1.0 - p) * s[c]);
            rate_->setRate(c, r[c]);
        }
    }
    model_->decomposeRateMatrix();
    if (rate_->getPInvar() > 0.0 || pinv_free_) tree_->computePtnInvar();
    tree_->clearAllPartialLH();
}

/* ---------------------------------------------------------------------- */
/* rate matrix as the kernel sees it                                       */
/* ---------------------------------------------------------------------- */

// Mirror of ModelMarkov::decomposeRateMatrix (Eigen3 path): states with a
// frequency <= ZERO_FREQ are dropped, pi is renormalised over the rest,
// Q_ab = R_ab pi_b, rows sum to zero, scale = total_num_subst / (pi^T R pi).
void ModelParamMap::buildQFrom(ModelMarkov *comp, const double *entries, const double *freq_raw, vector<double> &Q) const {
    const int S = nstates_;
    Q.assign((size_t)S * S, 0.0);
    vector<double> pi(S, 0.0);
    double sum = 0.0;
    for (int a = 0; a < S; a++) if (freq_raw[a] > ZERO_FREQ) sum += freq_raw[a];
    for (int a = 0; a < S; a++) if (freq_raw[a] > ZERO_FREQ) pi[a] = freq_raw[a] / sum;
    int k = 0;
    for (int a = 0; a < S; a++)
        for (int b = a + 1; b < S; b++, k++) {
            if (pi[a] == 0.0 || pi[b] == 0.0) continue;
            Q[a * S + b] = entries[k] * pi[b];
            Q[b * S + a] = entries[k] * pi[a];
        }
    double Z = 0.0;
    for (int a = 0; a < S; a++) {
        double row = 0.0;
        for (int b = 0; b < S; b++) if (b != a) row += Q[a * S + b];
        Q[a * S + a] = -row;
        Z += row * pi[a];
    }
    if (comp->normalize_matrix && Z > 0.0) {
        double scale = comp->total_num_subst / Z;
        for (double &q : Q) q *= scale;
    }
}

void ModelParamMap::buildQ(int m, vector<double> &Q) const {
    ModelMarkov *comp = component(m);
    vector<double> entries(comp->getNumRateEntries());
    comp->getRateMatrix(entries.data());
    buildQFrom(comp, entries.data(), comp->state_freq, Q);
}

double ModelParamMap::qIdentity(const PhyloGradient::Result &res) const {
    double total = 0.0;
    vector<double> Q;
    for (int m = 0; m < nmix_; m++) {
        if (res.dlogl_dQ[m].empty()) continue;
        buildQ(m, Q);
        for (size_t e = 0; e < Q.size(); e++) total += res.dlogl_dQ[m][e] * Q[e];
    }
    return total;
}

/* ---------------------------------------------------------------------- */
/* chain rule                                                              */
/* ---------------------------------------------------------------------- */

void ModelParamMap::naturalFromQ(const PhyloGradient::Result &res, bool fd_q) {
    const int S = nstates_;
    nat_.dR.assign(nmix_, vector<double>());
    nat_.dpi.assign(nmix_, vector<double>(S, 0.0));
    nat_.dw.assign(nmix_, 0.0);
    nat_.drate.assign(ncat_, 0.0);
    nat_.dprop.assign(ncat_, 0.0);

    for (int m = 0; m < nmix_; m++) {
        ModelMarkov *comp = component(m);
        int nrate = comp->getNumRateEntries();
        nat_.dR[m].assign(nrate, 0.0);
        vector<double> entries(nrate);
        comp->getRateMatrix(entries.data());
        const double *raw = comp->state_freq;
        double raw_sum = 0.0;
        for (int a = 0; a < S; a++) raw_sum += raw[a];
        // explicit dependence on pi at the root and in the invariant-site term
        vector<double> dpi_hat(S, 0.0);
        double w_m = model_->getMixtureWeight(m);
        for (int k = 0; k < S; k++)
            dpi_hat[k] = res.root_term[m][k] + w_m * res.inv_state_sum[k];

        if (!res.dlogl_dQ[m].empty()) {
            const double *D = res.dlogl_dQ[m].data();
            if (fd_q) {
                vector<double> Qp, Qm;
                for (int e = 0; e < nrate; e++) {
                    double h = 1e-5 * max(1.0, fabs(entries[e]));
                    vector<double> ep(entries), em(entries);
                    ep[e] += h; em[e] -= h;
                    buildQFrom(comp, ep.data(), raw, Qp);
                    buildQFrom(comp, em.data(), raw, Qm);
                    double s = 0.0;
                    for (int i = 0; i < S * S; i++) s += D[i] * (Qp[i] - Qm[i]);
                    nat_.dR[m][e] = s / (2 * h);
                }
                // raw frequency entries: the builder renormalises, so this is
                // already the raw-entry derivative of the Q part
                vector<double> dpi_raw_q(S, 0.0);
                for (int k = 0; k < S; k++) {
                    if (raw[k] <= ZERO_FREQ) continue;
                    double h = 1e-6 * max(raw_sum, raw[k]);
                    vector<double> fp(raw, raw + S), fm(raw, raw + S);
                    fp[k] += h; fm[k] -= h;
                    buildQFrom(comp, entries.data(), fp.data(), Qp);
                    buildQFrom(comp, entries.data(), fm.data(), Qm);
                    double s = 0.0;
                    for (int i = 0; i < S * S; i++) s += D[i] * (Qp[i] - Qm[i]);
                    dpi_raw_q[k] = s / (2 * h);
                }
                // project the explicit part through pi_hat = pi / sum(pi)
                double dot = 0.0;
                for (int j = 0; j < S; j++) dot += (raw[j] / raw_sum) * dpi_hat[j];
                for (int k = 0; k < S; k++)
                    nat_.dpi[m][k] = dpi_raw_q[k] + (dpi_hat[k] - dot) / raw_sum;
            } else {
                // closed form (design doc 6.2): Q = nu R Pi, nu = T / (pi^T R pi)
                vector<double> pi(S, 0.0), R((size_t)S * S, 0.0), Q;
                double act_sum = 0.0;
                for (int a = 0; a < S; a++) if (raw[a] > ZERO_FREQ) act_sum += raw[a];
                for (int a = 0; a < S; a++) if (raw[a] > ZERO_FREQ) pi[a] = raw[a] / act_sum;
                int k = 0;
                for (int a = 0; a < S; a++)
                    for (int b = a + 1; b < S; b++, k++)
                        if (pi[a] > 0.0 && pi[b] > 0.0) R[a * S + b] = R[b * S + a] = entries[k];
                buildQFrom(comp, entries.data(), raw, Q);
                vector<double> Rpi(S, 0.0);
                double Z = 0.0;
                for (int a = 0; a < S; a++) {
                    for (int b = 0; b < S; b++) Rpi[a] += R[a * S + b] * pi[b];
                    Z += pi[a] * Rpi[a];
                }
                bool norm = comp->normalize_matrix && Z > 0.0;
                double nu = norm ? comp->total_num_subst / Z : 1.0;
                double DQ = 0.0;
                for (int i = 0; i < S * S; i++) DQ += D[i] * Q[i];
                k = 0;
                for (int a = 0; a < S; a++)
                    for (int b = a + 1; b < S; b++, k++) {
                        if (pi[a] == 0.0 || pi[b] == 0.0) continue;
                        double g = nu * (D[a * S + b] * pi[b] + D[b * S + a] * pi[a]
                                         - D[a * S + a] * pi[b] - D[b * S + b] * pi[a]);
                        if (norm) g -= 2.0 * pi[a] * pi[b] / Z * DQ;
                        nat_.dR[m][k] = g;
                    }
                vector<double> dhat(S, 0.0);   // d/d pi_hat_k of the Q part
                for (int kk = 0; kk < S; kk++) {
                    if (pi[kk] == 0.0) continue;
                    double g = 0.0;
                    for (int a = 0; a < S; a++)
                        if (a != kk) g += R[a * S + kk] * (D[a * S + kk] - D[a * S + a]);
                    g *= nu;
                    if (norm) g -= 2.0 * Rpi[kk] / Z * DQ;
                    dhat[kk] = g;
                }
                // total d/d pi_hat, then the raw-entry projection
                for (int kk = 0; kk < S; kk++) dhat[kk] += dpi_hat[kk];
                double dot = 0.0;
                for (int j = 0; j < S; j++) dot += (raw[j] / raw_sum) * dhat[j];
                for (int kk = 0; kk < S; kk++) nat_.dpi[m][kk] = (dhat[kk] - dot) / raw_sum;
            }
        } else {
            double dot = 0.0;
            for (int j = 0; j < S; j++) dot += (raw[j] / raw_sum) * dpi_hat[j];
            for (int k = 0; k < S; k++) nat_.dpi[m][k] = (dpi_hat[k] - dot) / raw_sum;
        }
        // mixture weight (raw entry): classes of this component plus the invariant term
        double g_w = 0.0;
        for (int r = 0; r < ncat_; r++) g_w += rate_->getProp(r) * res.dlogl_dclass[m * ncat_ + r];
        for (int k = 0; k < S; k++) g_w += (raw[k] / raw_sum) * res.inv_state_sum[k];
        nat_.dw[m] = g_w;
    }
    for (int r = 0; r < ncat_; r++) {
        nat_.drate[r] = res.dlogl_drate[r];
        double g = 0.0;
        for (int m = 0; m < nmix_; m++) g += model_->getMixtureWeight(m) * res.dlogl_dclass[m * ncat_ + r];
        nat_.dprop[r] = g;
    }
    // invariant sites: direct term plus the coupled rates (x 1/(1-p)) and proportions (x (1-p))
    double p = rate_->getPInvar();
    nat_.dpinv_direct = res.inv_total;
    nat_.dpinv = res.inv_total;
    if (p > 0.0 || pinv_free_) {
        for (int c = 0; c < ncat_; c++)
            nat_.dpinv += nat_.drate[c] * rate_->getRate(c) / (1.0 - p) - nat_.dprop[c] * rate_->getProp(c) / (1.0 - p);
    }
    nat_.dalpha = 0.0;
    if (alpha_free_) {
        rateModelJacobian();
        for (int c = 0; c < ncat_; c++) nat_.dalpha += nat_.drate[c] * drate_dalpha_[c];
    }
}

// d rate_c / d alpha by Richardson-extrapolated central differences of the
// gamma quantile computation. RateGamma::computeRates rescales relative to
// the rates currently stored, so they are restored before EVERY evaluation.
void ModelParamMap::rateModelJacobian() {
    drate_dalpha_.assign(ncat_, 0.0);
    double alpha = rate_->getGammaShape();
    vector<double> r0(ncat_);
    for (int c = 0; c < ncat_; c++) r0[c] = rate_->getRate(c);
    auto restore = [&]() { for (int c = 0; c < ncat_; c++) rate_->setRate(c, r0[c]); };
    auto eval = [&](double a, vector<double> &out) {
        restore();
        rate_->setGammaShape(a);
        out.resize(ncat_);
        for (int c = 0; c < ncat_; c++) out[c] = rate_->getRate(c);
    };
    double h = 1e-4 * max(1.0, alpha);
    vector<double> rp, rm, rp2, rm2;
    eval(alpha + h, rp); eval(alpha - h, rm);
    eval(alpha + h / 2, rp2); eval(alpha - h / 2, rm2);
    for (int c = 0; c < ncat_; c++) {
        double d1 = (rp[c] - rm[c]) / (2 * h), d2 = (rp2[c] - rm2[c]) / h;
        drate_dalpha_[c] = (4 * d2 - d1) / 3;
    }
    restore();
    rate_->setGammaShape(alpha);
    restore();
}

void ModelParamMap::naturalToTheta(vector<double> &g) const {
    g.assign(ndim(), 0.0);
    for (auto &blk : sblocks_) {
        vector<double> v;
        readRateVariables(component(blk.comps[0]), v);
        for (int i = 0; i < blk.nvar; i++) {
            double s = 0.0;
            for (int m : blk.comps)
                for (int e : blk.entries[i]) s += nat_.dR[m][e];
            g[blk.offset + i] = s * v[i + 1];               // theta = log(variable)
        }
    }
    for (auto &blk : fblocks_) {
        const double *raw = component(blk.comp)->state_freq;
        for (size_t j = 0; j < blk.states.size(); j++) {
            int k = blk.states[j];
            g[blk.offset + j] = nat_.dpi[blk.comp][k] * raw[k];  // raw entry = e^theta * pinned
        }
    }
    if (w_offset_ >= 0) {
        double dot = 0.0;
        for (int m = 0; m < nmix_; m++) dot += model_->getMixtureWeight(m) * nat_.dw[m];
        for (int m = 0; m < nmix_ - 1; m++) {
            double w = model_->getMixtureWeight(m);
            g[w_offset_ + m] = w * (nat_.dw[m] - dot);
        }
    }
    if (rates_free_ || props_free_) {
        double p = rate_->getPInvar();
        vector<double> s(ncat_), r(ncat_);
        for (int c = 0; c < ncat_; c++) { s[c] = rate_->getProp(c) / (1.0 - p); r[c] = rate_->getRate(c); }
        double sum_gr = 0.0, sum_gp_s = 0.0;
        for (int c = 0; c < ncat_; c++) { sum_gr += nat_.drate[c] * r[c]; sum_gp_s += nat_.dprop[c] * s[c]; }
        if (rates_free_)
            for (int k = 0; k < ncat_ - 1; k++)
                g[r_rate_offset_ + k] = nat_.drate[k] * r[k] - (1.0 - p) * s[k] * r[k] * sum_gr;
        if (props_free_)
            for (int k = 0; k < ncat_ - 1; k++) {
                double gp = (1.0 - p) * s[k] * (nat_.dprop[k] - sum_gp_s);
                double gr = -s[k] * ((1.0 - p) * r[k] - 1.0) * sum_gr;
                g[r_prop_offset_ + k] = gp + gr;
            }
    }
    if (alpha_free_) g[alpha_offset_] = nat_.dalpha * rate_->getGammaShape();
    if (pinv_free_) { double p = rate_->getPInvar(); g[pinv_offset_] = nat_.dpinv * p * (1.0 - p); }
}

void ModelParamMap::gradient(const PhyloGradient::Result &res, vector<double> &g) {
    naturalFromQ(res, false);
    naturalToTheta(g);
}

void ModelParamMap::gradientFDQ(const PhyloGradient::Result &res, vector<double> &g) {
    naturalFromQ(res, true);
    naturalToTheta(g);
}
