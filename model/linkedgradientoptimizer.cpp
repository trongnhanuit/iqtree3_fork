/*
 * linkedgradientoptimizer.cpp
 *
 * See linkedgradientoptimizer.h and docs/analytical-gradients-design.md, section 15.
 */

#include "linkedgradientoptimizer.h"
#include "gradientoptimizer.h"
#include "modelparammap.h"
#include "phylogradient.h"
#include "modelfactory.h"
#include "modelmarkov.h"
#include "modelmixture.h"
#include "tree/phylosupertree.h"
#include "utils/timeutil.h"
#include "utils/tools.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>

using namespace std;

namespace {
ModelMarkov *component(ModelSubst *model, int m) {
    if (model->isMixture()) {
        ModelMixture *mix = dynamic_cast<ModelMixture*>(model);
        return mix ? mix->at(m) : nullptr;
    }
    return dynamic_cast<ModelMarkov*>(model);
}
}

LinkedGradientOptimizer::LinkedGradientOptimizer(PhyloSuperTree *stree, const string &model_name)
    : stree_(stree), name_(model_name) {
    dfpmin_scaled = true;   // Shanno-Phua scaling and the curvature-sign fix (design 9)
}

LinkedGradientOptimizer::~LinkedGradientOptimizer() {
}

int LinkedGradientOptimizer::getNDim() {
    return maps_.empty() ? 0 : maps_[0]->ndim();
}

/* ---------------------------------------------------------------------- */
/* entry checks                                                            */
/* ---------------------------------------------------------------------- */

string LinkedGradientOptimizer::prepare() {
    trees_.clear();
    maps_.clear();
    if (stree_->part_order.empty()) stree_->computePartitionOrder();
    vector<int> index_of(stree_->size(), -1);
    for (int i = 0; i < (int)stree_->size(); i++)
        if (stree_->at(i)->getModel()->getName() == name_) {
            index_of[i] = (int)trees_.size();
            trees_.push_back(stree_->at(i));
        }
    if (trees_.empty()) return "no partition uses the model";
    order_.clear();
    for (int i : stree_->part_order)
        if (index_of[i] >= 0) order_.push_back(index_of[i]);
    checkCopiesIdentical();
    for (PhyloTree *t : trees_) {
        string why;
        if (!GradientOptimizer::supports(t->getModelFactory(), why)) return t->aln->name + ": " + why;
    }
    // maps are rebuilt for every step: a map pins each profile's largest state when built
    maps_.resize(trees_.size());
    for (size_t p = 0; p < trees_.size(); p++) {
        maps_[p].reset(new ModelParamMap(trees_[p]->getModelFactory(), trees_[p], true));
        const ModelParamMap &m = *maps_[p], &m0 = *maps_[0];
        if (!m.unsupported().empty()) return trees_[p]->aln->name + ": " + m.unsupported();
        if (m.ndim() != m0.ndim() || m.lower() != m0.lower() || m.upper() != m0.upper())
            return "different parameter layouts across partitions";
        for (int i = 0; i < m.ndim(); i++)
            if (m.name(i) != m0.name(i)) return "different parameter layouts across partitions";
    }
    // states at zero frequency drop out of the map but not out of the model
    if (maps_[0]->ndim() != trees_[0]->getModel()->getNDim())
        return "parameter layout differs from the model's (zero-frequency states)";
    lh_part_.assign(trees_.size(), 0.0);
    g_part_.assign(trees_.size(), vector<double>());
    return "";
}

void LinkedGradientOptimizer::checkCopiesIdentical() const {
    ModelSubst *a = trees_[0]->getModel();
    const int S = a->num_states;
    for (size_t p = 1; p < trees_.size(); p++) {
        ModelSubst *b = trees_[p]->getModel();
        bool same = a->getNMixtures() == b->getNMixtures();
        for (int k = 0; same && k < a->getNMixtures(); k++) {
            ModelMarkov *ca = component(a, k), *cb = component(b, k);
            if (!ca || !cb) { same = false; break; }
            vector<double> ra(ca->getNumRateEntries()), rb(cb->getNumRateEntries());
            ca->getRateMatrix(ra.data());
            cb->getRateMatrix(rb.data());
            same = ra == rb && memcmp(ca->state_freq, cb->state_freq, sizeof(double) * S) == 0
                   && ca->total_num_subst == cb->total_num_subst
                   && a->getMixtureWeight(k) == b->getMixtureWeight(k);
        }
        if (!same)
            outError("Linked model " + name_ + ": partitions " + trees_[0]->aln->name + " and " + trees_[p]->aln->name +
                     " hold different values. With --analytical-gradients, start the profiles with --start-points "
                     "warm|cold on GTR20+Fk, or give them explicitly, e.g. --model-joint \"GTR20+FMIX{C10pi1,...}\" "
                     "--link-exchange-rates -mfopt");
    }
}

/* ---------------------------------------------------------------------- */
/* Optimization interface                                                  */
/* ---------------------------------------------------------------------- */

double LinkedGradientOptimizer::sumLikelihood() {
    const int np = (int)trees_.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if(stree_->num_threads > 1)
#endif
    for (int j = 0; j < np; j++) {
        int p = order_[j];
        lh_part_[p] = trees_[p]->computeLikelihood();
    }
    double sum = 0.0;   // fixed order, so the result does not depend on the thread count
    for (int p = 0; p < np; p++) sum += lh_part_[p];
    return sum;
}

double LinkedGradientOptimizer::targetFunk(double x[]) {
    double t = getRealTime();
    const int n = getNDim(), np = (int)trees_.size();
    theta_.assign(x + 1, x + 1 + n);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if(stree_->num_threads > 1)
#endif
    for (int j = 0; j < np; j++) {
        int p = order_[j];
        maps_[p]->unpack(theta_);
        lh_part_[p] = trees_[p]->computeLikelihood();
    }
    double sum = 0.0;
    for (int p = 0; p < np; p++) sum += lh_part_[p];
    n_lh_++;
    t_lh_ += getRealTime() - t;
    if (!std::isfinite(sum)) return 1e30;   // a rejected line-search step, never a crash
    return -sum;
}

double LinkedGradientOptimizer::derivativeFunk(double x[], double dfx[]) {
    double t = getRealTime();
    const int n = getNDim(), np = (int)trees_.size();
    theta_.assign(x + 1, x + 1 + n);
    vector<char> ok(np, 0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if(stree_->num_threads > 1)
#endif
    for (int j = 0; j < np; j++) {
        int p = order_[j];
        maps_[p]->unpack(theta_);
        // one engine per call: its outside-partial pool lives only while this partition runs
        PhyloGradient engine(trees_[p]);
        engine.setNeedQ(maps_[p]->needQ());
        PhyloGradient::Result res;
        try {
            ok[p] = engine.compute(res);
        } catch (std::exception &e) {
            ok[p] = 0;   // the engine's guard has restored the tree; this call falls back below
        }
        if (ok[p]) {
            lh_part_[p] = res.logl;
            maps_[p]->gradient(res, g_part_[p]);   // raw: the floor rule is applied to the sum
        }
    }
    bool all_ok = GradientOptimizer::eigenResidualOk(trees_[0]);   // the copies share one Q
    for (int p = 0; p < np && all_ok; p++) all_ok = ok[p];
    n_grad_++;
    if (!all_ok) {
        // joint finite differences for this call, as GradientOptimizer does for one tree
        n_fd_++;
        if (!warned_fd_) {
            GradientOptimizer::say("NOTE: analytic gradient of the linked step unavailable at this point; "
                                   "using finite differences for this step");
            warned_fd_ = true;
        }
        double f = Optimization::derivativeFunk(x, dfx);
        targetFunk(x);   // every copy back at x
        t_grad_ += getRealTime() - t;
        return f;
    }
    vector<double> g(n, 0.0);
    double sum = 0.0;
    for (int p = 0; p < np; p++) {
        for (int i = 0; i < n; i++) g[i] += g_part_[p][i];
        sum += lh_part_[p];
    }
    // flat below a floor (design 6): applied once to the sum; the copies hold the same
    // state, so partition 0's flags are everyone's
    const vector<char> &clamped = maps_[0]->clamped();
    for (int i = 0; i < n; i++) {
        if (clamped[i] && g[i] < 0.0) g[i] = 0.0;
        dfx[i + 1] = -g[i];
    }
    t_grad_ += getRealTime() - t;
    return -sum;
}

bool LinkedGradientOptimizer::stopEarly(int iter, double f_prev, double f_new) {
    // NOTE(design 9): stop once a gain is below logl_epsilon and below 1% of the largest gain
    double delta = f_prev - f_new;
    if (iter == 1) { delta_max_ = delta; stop_iter_ = 1; return false; }
    delta_max_ = max(delta_max_, delta);
    stop_iter_ = iter;
    return delta < min(0.01 * delta_max_, logl_epsilon_);
}

/* ---------------------------------------------------------------------- */
/* joint gradient check                                                    */
/* ---------------------------------------------------------------------- */

int LinkedGradientOptimizer::gradientCheck() {
    const int n = getNDim();
    const double tol = Params::getInstance().ag_gradient_check_tol;
    vector<double> theta;
    maps_[0]->pack(theta);
    vector<double> x(n + 1), dfx(n + 1);
    for (int i = 0; i < n; i++) x[i + 1] = theta[i];
    double f0 = derivativeFunk(x.data(), dfx.data());
    double gmax = 0.0;
    for (int i = 1; i <= n; i++) gmax = max(gmax, fabs(dfx[i]));
    auto fd = [&](int i, double h) {
        double xi = x[i];
        x[i] = xi + h; double fp = targetFunk(x.data());
        x[i] = xi - h; double fm = targetFunk(x.data());
        x[i] = xi;
        return -(fp - fm) / (2 * h);   // targetFunk is -lnL
    };
    ofstream out(Params::getInstance().out_prefix + string(".gradcheck.tsv"));
    out << "idx\tname\tx\tanalytic\tfd_h\tfd_h2\trichardson\trel_err\tstatus" << endl;
    int n_fail = 0, n_noise = 0, worst = 0;
    double max_rel = 0.0;
    for (int i = 1; i <= n; i++) {
        double h = 1e-4 * max(1.0, fabs(x[i]));
        double a = -dfx[i], d1 = fd(i, h), d2 = fd(i, h / 2), r = (4 * d2 - d1) / 3;
        double rel = fabs(a - r) / max(max(fabs(a), fabs(r)), 1e-6 * gmax);
        // a disagreement no larger than the FD step's own instability is rounding noise
        string status = rel <= tol ? "PASS" : (fabs(a - r) <= 4 * fabs(d1 - d2) ? "PASS-FDNOISE" : "FAIL");
        if (status == "FAIL") n_fail++;
        if (status == "PASS-FDNOISE") n_noise++;
        if (rel > max_rel && status != "PASS-FDNOISE") { max_rel = rel; worst = i - 1; }
        out << i - 1 << "\t" << maps_[0]->name(i - 1) << "\t" << x[i] << "\t" << a << "\t" << d1 << "\t" << d2
            << "\t" << r << "\t" << rel << "\t" << status << endl;
    }
    targetFunk(x.data());   // every copy back at the checked point
    ostringstream os;
    os << "GRADCHECK joint n=" << n << " partitions=" << trees_.size() << " logl=" << fixed << setprecision(4) << -f0
       << " max_rel=" << scientific << setprecision(2) << max_rel << " worst=" << (n ? maps_[0]->name(worst) : "")
       << " n_fail=" << n_fail << " n_fdnoise=" << n_noise;
    GradientOptimizer::say(os.str());
    return n_fail;
}

/* ---------------------------------------------------------------------- */
/* model file                                                              */
/* ---------------------------------------------------------------------- */

void LinkedGradientOptimizer::writeModelFile(const string &file) const {
    ModelSubst *model = trees_[0]->getModel();
    const int K = model->getNMixtures(), S = model->num_states;
    ModelMarkov *c0 = component(model, 0);
    if (!c0 || !c0->isReversible() || c0->getNumRateEntries() != S * (S - 1) / 2) return;
    vector<double> rates(c0->getNumRateEntries()), other(rates.size());
    c0->getRateMatrix(rates.data());
    for (int k = 1; k < K; k++) {   // one matrix only for linked exchangeabilities
        component(model, k)->getRateMatrix(other.data());
        if (other != rates) return;
    }
    if (rates.back() > 0.0) { double pin = rates.back(); for (double &r : rates) r /= pin; }
    vector<vector<double>> pi(K, vector<double>(S));
    vector<double> w(K), mean(S, 0.0);
    double wsum = 0.0;
    for (int k = 0; k < K; k++) {
        const double *f = component(model, k)->state_freq;
        double s = 0.0;
        for (int a = 0; a < S; a++) s += f[a];
        for (int a = 0; a < S; a++) pi[k][a] = f[a] / s;   // training keeps them unnormalised
        w[k] = model->getMixtureWeight(k);
        wsum += w[k];
    }
    for (int k = 0; k < K; k++) {
        w[k] /= wsum;
        for (int a = 0; a < S; a++) mean[a] += w[k] * pi[k][a];
    }
    vector<int> by_weight(K);
    for (int k = 0; k < K; k++) by_weight[k] = k;
    stable_sort(by_weight.begin(), by_weight.end(), [&](int a, int b) { return w[a] > w[b]; });
    // the protein-matrix lookup upper-cases the name, so the matrix name must be upper case
    const string name = "MM" + convertIntToString(K), qname = K > 1 ? name + "Q" : name;
    ofstream out(file);
    out << setprecision(17);
    out << "#nexus" << endl
        << "[ mixture_maker: " << (K > 1 ? convertIntToString(K) + " profiles with one" : "one")
        << " exchangeability matrix, estimated jointly from " << trees_.size() << " alignments. " << endl
        << "  The matrix belongs with its profiles; it is not a standalone matrix. Exchangeabilities" << endl
        << "  at the lower bound (1e-4 of the last entry) are not estimable from the data." << endl
        << "  Use: -mdef " << file << " -m " << name << "+G4 ]" << endl
        << "begin models;" << endl
        << "model " << qname << " =" << endl;
    // PAML lower triangle (row a, columns b < a); the rates are stored as the upper triangle
    for (int a = 1; a < S; a++) {
        for (int b = 0; b < a; b++) out << " " << rates[b * S - b * (b + 1) / 2 + (a - b - 1)];
        out << endl;
    }
    for (int a = 0; a < S; a++) out << " " << (K > 1 ? mean[a] : pi[0][a]);
    out << ";" << endl;
    if (K > 1) {
        for (int j = 0; j < K; j++) {
            out << "frequency " << name << "pi" << j + 1 << " =";
            for (int a = 0; a < S; a++) out << " " << pi[by_weight[j]][a];
            out << ";" << endl;
        }
        out << "model " << name << " = " << qname << "+FMIX{";
        for (int j = 0; j < K; j++) {
            int k = by_weight[j];
            out << (j ? "," : "") << name << "pi" << j + 1 << ":" << component(model, k)->total_num_subst << ":" << w[k];
        }
        out << "};" << endl;
    }
    out << "end;" << endl;
}

/* ---------------------------------------------------------------------- */
/* one linked step                                                         */
/* ---------------------------------------------------------------------- */

void LinkedGradientOptimizer::snapshot(RawState &st) const {
    ModelSubst *model = trees_[0]->getModel();
    const int K = model->getNMixtures(), S = model->num_states;
    st.rates.assign(K, vector<double>());
    st.freqs.assign(K, vector<double>());
    st.weights.assign(K, 1.0);
    for (int k = 0; k < K; k++) {
        ModelMarkov *c = component(model, k);
        st.rates[k].resize(c->getNumRateEntries());
        c->getRateMatrix(st.rates[k].data());
        st.freqs[k].assign(c->state_freq, c->state_freq + S);
        st.weights[k] = model->getMixtureWeight(k);
    }
}

void LinkedGradientOptimizer::restore(const RawState &st) {
    for (PhyloTree *t : trees_) {
        ModelSubst *model = t->getModel();
        for (int k = 0; k < model->getNMixtures(); k++) {
            ModelMarkov *c = component(model, k);
            vector<double> rates(st.rates[k]), freqs(st.freqs[k]);   // the setters take non-const pointers
            c->setRateMatrix(rates.data());
            c->setStateFrequency(freqs.data());
            model->setMixtureWeight(k, st.weights[k]);
        }
        model->decomposeRateMatrix();
        if (t->getRate()->getPInvar() > 0.0) t->computePtnInvar();
        t->clearAllPartialLH();
    }
}

double LinkedGradientOptimizer::step(bool write_info, double logl_epsilon) {
    double t0 = getRealTime();
    double t_step_b = last_end_ > 0.0 ? t0 - last_end_ : 0.0;   // partitions' own step since the last linked step
    const int n = getNDim();
    logl_epsilon_ = logl_epsilon;
    vector<double> theta, x(n + 1), lower(n + 1), upper(n + 1);
    maps_[0]->pack(theta);
    for (int i = 0; i < n; i++) {
        x[i + 1] = theta[i];
        lower[i + 1] = maps_[0]->lower()[i];
        upper[i + 1] = maps_[0]->upper()[i];
    }
    RawState entry_state;
    snapshot(entry_state);
    double entry = -targetFunk(x.data());
    long lh0 = n_lh_, g0 = n_grad_, fd0 = n_fd_;
    double tl0 = t_lh_, tg0 = t_grad_;
    delta_max_ = 0.0;
    stop_iter_ = 0;
    // dfpmin's own relative gradient test would stop at |g| ~ 1e-4 |logL|; with gtol
    // negligible the search ends through stopEarly, the step-size test or the cap
    unique_ptr<bool[]> bc(new bool[n + 1]());   // no bound checks: restartParameters() never fires
    minimizeMultiDimen(x.data(), n, lower.data(), upper.data(), bc.get(), 1e-10);
    double final_lh = -targetFunk(x.data());   // the returned point, in every copy
    bool restored = false;
    if (final_lh < entry - 1e-6) {
        restore(entry_state);
        final_lh = sumLikelihood();
        restored = true;
    }
    if (write_info || Params::getInstance().ag_stats) {
        int up_hits = 0, floor_hits = 0;
        for (int i = 1; i <= n; i++) if (x[i] >= upper[i] - 1e-9) up_hits++;
        for (char c : maps_[0]->clamped()) floor_hits += c;
        if (maps_[0]->lastWeightClamped()) floor_hits++;
        long nl = n_lh_ - lh0, ng = n_grad_ - g0;
        ostringstream os;
        os << "AG joint: step " << steps_ + 1 << ", " << trees_.size() << " partitions, " << n << " parameters: "
           << stop_iter_ << " iterations, logL " << fixed << setprecision(3) << entry << " -> " << final_lh
           << (restored ? " (entry restored)" : "") << ", partitions' step " << setprecision(1) << t_step_b << " s"
           << ", upper bounds " << up_hits << ", floors " << floor_hits << ", fd fallbacks " << n_fd_ - fd0
           << ", lnL " << nl << " x " << setprecision(2) << (nl ? 1000.0 * (t_lh_ - tl0) / nl : 0.0) << " ms"
           << ", gradient " << ng << " x " << (ng ? 1000.0 * (t_grad_ - tg0) / ng : 0.0) << " ms";
        GradientOptimizer::say(os.str());
    }
    steps_++;
    writeModelFile(Params::getInstance().out_prefix + string(".mm.nex"));
    last_end_ = getRealTime();
    return final_lh;
}
