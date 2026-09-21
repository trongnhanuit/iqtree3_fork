/*
 * gradientoptimizer.cpp
 *
 * See gradientoptimizer.h and docs/analytical-gradients-design.md.
 */

#include "gradientoptimizer.h"
#include "modelfactory.h"
#include "modelmarkov.h"
#include "modelmixture.h"
#include "tree/phylotree.h"
#include "phylogradient.h"
#include "modelparammap.h"
#include "utils/tools.h"
#include "utils/MPIHelper.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <random>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

/* ---------------------------------------------------------------------- */
/* capability check                                                        */
/* ---------------------------------------------------------------------- */

bool GradientOptimizer::supports(ModelFactory *factory, std::string &why) {
    ModelSubst *model = factory->model;
    RateHeterogeneity *rate = factory->site_rate;
    PhyloTree *tree = rate->getTree();
    // NOTE(design 4): every exclusion of design section 2 with its reason
    if (tree->isSuperTree() || tree->isTreeMix()) { why = "partition or tree-mixture container"; return false; }
    if (tree->isMixlen() || rate->isHeterotachy()) { why = "heterotachy (+H) model"; return false; }
    if (factory->fused_mix_rate) { why = "fused mixture-rate model (*G/*R)"; return false; }
    if (tree->aln->seq_type == SEQ_CODON) { why = "codon model"; return false; }
    if (model->isPolymorphismAware()) { why = "polymorphism-aware (PoMo) model"; return false; }
    if (model->containDNAerror()) { why = "DNA error model"; return false; }
    if (!model->useRevKernel()) { why = "non-reversible model or kernel"; return false; }
    if (model->isSiteSpecificModel()) { why = "site-specific model"; return false; }
    if (rate->isSiteSpecificRate()) { why = "site-specific rates"; return false; }
    if (!factory->unobserved_ptns.empty()) { why = "ascertainment-bias correction (+ASC)"; return false; }
    if (Params::getInstance().lh_mem_save == LM_MEM_SAVE) { why = "memory-saving mode (-mem)"; return false; }
    if (MPIHelper::getInstance().getNumProcesses() > 1) { why = "MPI with more than one process"; return false; }
    if (model->getNDim() + rate->getNDim() == 0) { why = "no free model parameters"; return false; }
    ModelParamMap map(factory, tree);
    if (!map.unsupported().empty()) { why = map.unsupported(); return false; }
    if (map.ndim() == 0) { why = "no free model parameters"; return false; }
    return true;
}

/* ---------------------------------------------------------------------- */
/* helpers                                                                 */
/* ---------------------------------------------------------------------- */

namespace {

// leaf names on the `node` side of the branch (node, dad), sorted, joined by '|'
void collectLeafNames(Node *node, Node *dad, vector<string> &names) {
    if (node->isLeaf()) { names.push_back(node->name); return; }
    FOR_NEIGHBOR_DECLARE(node, dad, it)
        collectLeafNames((*it)->node, node, names);
}

string sideLabel(Node *node, Node *dad) {
    vector<string> names;
    collectLeafNames(node, dad, names);
    sort(names.begin(), names.end());
    string s;
    for (size_t i = 0; i < names.size(); i++) { if (i) s += '|'; s += names[i]; }
    return s;
}

double relErr(double a, double n, double gmax) {
    double d = max(max(fabs(a), fabs(n)), 1e-6 * gmax);
    return d > 0 ? fabs(a - n) / d : 0.0;
}

// long-double reference for the divided-difference kernel
long double xRef(long double a, long double b, long double tau) {
    long double d = a - b;
    if (d == 0.0L) return tau * expl(a * tau);
    if (fabsl(d * tau) < 0.1L) return expl(b * tau) * expm1l(tau * d) / d;
    return (expl(a * tau) - expl(b * tau)) / d;
}

struct CheckRow {
    string idx, name, label;
    double x, analytic, newton, fd_h, fd_h2, rich, rel;
    bool has_newton, pass, noisy;
};

struct CheckStats { int n_fail = 0, n_noisy = 0; double max_rel = 0.0; string worst; };

// A row whose two central-difference estimates disagree with each other by
// more than they disagree with the analytic value is limited by rounding in
// the likelihood (tiny gradients after a 20x20 eigendecomposition); it is
// counted as noisy, not as a failure.
void judgeRow(CheckRow &r, double tol, double gmax, CheckStats &st) {
    r.rel = relErr(r.analytic, r.rich, gmax);
    double fd_noise = fabs(r.fd_h - r.fd_h2);
    r.noisy = false;
    r.pass = std::isfinite(r.analytic) && r.rel <= tol;
    if (!r.pass && std::isfinite(r.analytic) && fabs(r.analytic - r.rich) <= 4.0 * fd_noise) { r.pass = true; r.noisy = true; st.n_noisy++; }
    if (!r.pass) st.n_fail++;
    if (r.rel > st.max_rel && !r.noisy) { st.max_rel = r.rel; st.worst = r.name; }
}

// central differences of the tree log-likelihood in theta space, through unpack()
void paramRows(ModelParamMap &map, PhyloTree *tree, const vector<double> &theta, const vector<double> &g,
               double tol, double gmax, vector<CheckRow> &rows, CheckStats &st) {
    for (int i = 0; i < map.ndim(); i++) {
        CheckRow r;
        r.idx = convertIntToString(i);
        r.name = map.name(i);
        r.x = theta[i];
        r.analytic = g[i];
        r.newton = NAN;
        r.has_newton = false;
        double h = 1e-4 * max(1.0, fabs(r.x));
        double lh[4];
        double steps[4] = { h, -h, h / 2, -h / 2 };
        for (int k = 0; k < 4; k++) {
            vector<double> th(theta);
            th[i] += steps[k];
            map.unpack(th);
            lh[k] = tree->computeLikelihood();
        }
        r.fd_h = (lh[0] - lh[1]) / (2 * h);
        r.fd_h2 = (lh[2] - lh[3]) / h;
        r.rich = (4 * r.fd_h2 - r.fd_h) / 3;
        judgeRow(r, tol, gmax, st);
        rows.push_back(r);
    }
    map.unpack(theta);
}

void writeRows(ostream &out, const vector<CheckRow> &rows, long iter, const string &level) {
    out << setprecision(10);
    for (auto &r : rows) {
        out << iter << '\t' << level << '\t' << r.idx << '\t' << r.name << '\t' << r.x << '\t'
            << r.analytic << '\t' << r.rich << '\t' << r.fd_h << '\t' << r.fd_h2 << '\t';
        if (r.has_newton) out << r.newton; else out << "NA";
        out << '\t' << fabs(r.analytic - r.rich) << '\t' << r.rel << '\t' << (r.noisy ? "PASS-FDNOISE" : r.pass ? "PASS" : "FAIL") << '\n';
    }
}

const char *GRADCHECK_HEADER = "iter\tlevel\tidx\tname\tx\tanalytic\tnumeric\tfd_h\tfd_h2\tnewton\tabs_err\trel_err\tstatus\n";

} // namespace

/* ---------------------------------------------------------------------- */
/* construction, Optimization interface                                    */
/* ---------------------------------------------------------------------- */

GradientOptimizer::GradientOptimizer(ModelFactory *factory)
    : factory_(factory), tree_(factory->site_rate->getTree()) {
    map_.reset(new ModelParamMap(factory, tree_));
    engine_.reset(new PhyloGradient(tree_));
    engine_->setNeedQ(map_->needQ());
}

GradientOptimizer::~GradientOptimizer() {
}

void GradientOptimizer::say(const string &line) const {
#ifdef _OPENMP
#pragma omp critical(ag_output)
#endif
    {
        cout << line << endl;
    }
}

int GradientOptimizer::getNDim() {
    return map_->ndim();
}

double GradientOptimizer::targetFunk(double x[]) {
    int n = map_->ndim();
    theta_.assign(x + 1, x + 1 + n);
    map_->unpack(theta_);
    n_lh_++;
    double lh = tree_->computeLikelihood();
    if (!std::isfinite(lh)) return 1e30;    // a rejected line-search step, never a crash (design 6.3)
    return -lh;
}

double GradientOptimizer::derivativeFunk(double x[], double dfx[]) {
    int n = map_->ndim();
    theta_.assign(x + 1, x + 1 + n);
    map_->unpack(theta_);
    PhyloGradient::Result res;
    bool ok = engine_->compute(res) && eigenResidualOk();
    n_grad_++;
    n_lh_++;   // the engine's forward pass is one likelihood evaluation
    if (!ok) {
        // NOTE(design 5.5): an invalid gradient (non-finite value or an
        // eigendecomposition that fails the residual check) falls back to the
        // base class's finite differences for this step, warning once.
        n_fd_fallback_++;
        if (!warned_fallback_) {
            say("NOTE: analytic gradient unavailable at this point (non-finite value or eigen residual); using finite differences for this step");
            warned_fallback_ = true;
        }
        n_lh_ += n + 1;
        return Optimization::derivativeFunk(x, dfx);
    }
    map_->gradient(res, grad_);
    for (int i = 0; i < n; i++) dfx[i + 1] = -grad_[i];
    if (Params::getInstance().ag_gradient_check &&
        (n_grad_ % max(1, Params::getInstance().ag_gradient_check_every)) == 0)
        gradientCheckStep(theta_, grad_);
    return -res.logl;
}

bool GradientOptimizer::stopEarly(int iter, double f_prev, double f_new) {
    // NOTE(design 6.4): stop when a step gains less than 1% of the largest
    // step so far AND less than the absolute floor (the floor is mandatory:
    // a slow ridge must not be abandoned while steps still exceed logl_epsilon)
    double delta = f_prev - f_new;
    if (iter == 1) { delta_max_ = delta; stop_iter_ = 1; return false; }
    delta_max_ = max(delta_max_, delta);
    stop_iter_ = iter;
    return delta < min(0.01 * delta_max_, logl_epsilon_);
}

bool GradientOptimizer::eigenResidualOk() const {
    ModelSubst *model = tree_->getModel();
    const size_t S = tree_->aln->num_states, stride = get_safe_upper_limit(S) * S;
    const double *U = model->getEigenvectors(), *Uinv = model->getInverseEigenvectors();
    for (int m = 0; m < model->getNMixtures(); m++) {
        const double *u = U + m * stride, *v = Uinv + m * stride;
        for (size_t a = 0; a < S; a++)
            for (size_t b = 0; b < S; b++) {
                double s = 0.0;
                for (size_t k = 0; k < S; k++) s += u[a * S + k] * v[k * S + b];
                if (fabs(s - (a == b ? 1.0 : 0.0)) > 1e-8) return false;
            }
    }
    return true;
}

/* ---------------------------------------------------------------------- */
/* state, symmetry breaking, polish                                        */
/* ---------------------------------------------------------------------- */

void GradientOptimizer::snapshot(BestState &st, double logl) {
    map_->pack(st.theta);
    st.brlen.clear();
    tree_->saveBranchLengths(st.brlen);
    st.logl = logl;
    st.valid = true;
}

void GradientOptimizer::restore(const BestState &st) {
    map_->unpack(st.theta);
    DoubleVector brlen(st.brlen);
    tree_->restoreBranchLengths(brlen);
    tree_->clearAllPartialLH();
}

void GradientOptimizer::breakSymmetry(bool write_info) {
    ModelSubst *model = tree_->getModel();
    if (!model->isMixture()) return;
    ModelMixture *mix = dynamic_cast<ModelMixture*>(model);
    if (!mix) return;
    const int S = tree_->aln->num_states;
    vector<ModelMarkov*> est;
    for (size_t m = 0; m < mix->size(); m++)
        if (mix->at(m)->getFreqType() == FREQ_ESTIMATE && mix->at(m)->getNDim() > 0) est.push_back(mix->at(m));
    if (est.size() < 2) return;
    for (size_t m = 1; m < est.size(); m++)
        for (int k = 0; k < S; k++)
            if (fabs(est[m]->state_freq[k] - est[0]->state_freq[k]) > 1e-12) return;   // already distinct
    // NOTE(design 6.4): identical classes are a fixed point of every gradient
    // method (the gradient is identical for all of them), so the `+Fk` default
    // start must be perturbed. Protein: the first k profiles of the smallest
    // C-series with >= k classes; otherwise a light log-normal jitter from a
    // private, seeded generator (no global RNG state touched).
    string how;
    bool done = false;
    if (S == 20 && est.size() <= 60) {
        int kk = 10;
        while (kk < (int)est.size()) kk += 10;
        ModelsBlock *models_block = readModelsDefinition(Params::getInstance());
        bool all = true;
        vector<string> descr;
        for (size_t m = 0; m < est.size() && all; m++) {
            NxsModel *fm = models_block->findModel("C" + convertIntToString(kk) + "pi" + convertIntToString((int)m + 1));
            if (!fm || !(fm->flag & NM_FREQ)) all = false; else descr.push_back(fm->description);
        }
        delete models_block;
        if (all) {
            for (size_t m = 0; m < est.size(); m++) est[m]->readStateFreq(descr[m]);
            how = "C" + convertIntToString(kk) + " profiles";
            done = true;
        }
    }
    if (!done) {
        std::mt19937_64 rng((unsigned long long)Params::getInstance().ran_seed + 7919ULL);
        std::normal_distribution<double> normal(0.0, 1.0);
        const double A = 0.1;
        for (size_t m = 0; m < est.size(); m++) {
            double sum = 0.0;
            for (int k = 0; k < S; k++) { est[m]->state_freq[k] *= exp(A * normal(rng)); sum += est[m]->state_freq[k]; }
            for (int k = 0; k < S; k++) est[m]->state_freq[k] /= sum;
        }
        how = "log-normal jitter";
    }
    model->decomposeRateMatrix();
    if (tree_->getRate()->getPInvar() > 0.0) tree_->computePtnInvar();
    tree_->clearAllPartialLH();
    if (write_info || verbose_mode >= VB_MED)
        say("AG: " + convertIntToString((int)est.size()) + " identical starting profiles; symmetry broken with " + how);
}

double GradientOptimizer::polish(double gradient_epsilon) {
    int n = map_->ndim();
    vector<double> x(n + 1), lower(n + 1), upper(n + 1);
    bool *bound_check = new bool[n + 1];
    for (int i = 0; i <= n; i++) bound_check[i] = false;
    vector<double> theta;
    map_->pack(theta);
    for (int i = 0; i < n; i++) { x[i + 1] = theta[i]; lower[i + 1] = map_->lower()[i]; upper[i + 1] = map_->upper()[i]; }
    delta_max_ = 0.0;
    stop_iter_ = 0;
    // NOTE(design 9): dfpmin's own test stops when |g_i| max(|x_i|,1) / |f| < gtol,
    // i.e. at |g| ~ 0.5 per parameter for |logL| ~ 5000 with the default 1e-4;
    // on a slow ridge that leaves several log-likelihood units on the table
    // (the default path creeps along the same ridge over dozens of rounds).
    // Gradients are cheap here, so the test is made negligible and the search
    // ends through stopEarly() (gain below logl_epsilon), the step-size test or
    // the iteration cap instead.
    double gtol = 1e-10;
    (void)gradient_epsilon;
    if (Params::getInstance().ag_optalg == "LBFGSB") {
        L_BFGS_B(n, x.data() + 1, lower.data() + 1, upper.data() + 1, gtol, 200, 20);
    } else {
        // bound_check all false: restartParameters() never fires anyway
        minimizeMultiDimen(x.data(), n, lower.data(), upper.data(), bound_check, gtol);
    }
    n_bfgs_iter_ += stop_iter_;
    delete [] bound_check;
    theta.assign(x.begin() + 1, x.end());
    map_->unpack(theta);
    n_lh_++;
    return tree_->computeLikelihood();
}

/* ---------------------------------------------------------------------- */
/* the replacement for the alternating loop                                */
/* ---------------------------------------------------------------------- */

double GradientOptimizer::optimize(int fixed_len, bool write_info, double logl_epsilon, double gradient_epsilon, double entry_logl) {
    Params &params = Params::getInstance();
    logl_epsilon_ = logl_epsilon;
    rounds_ = 0;
    const double t0 = getRealTime();

    BestState entry, best;
    snapshot(entry, entry_logl);
    best = entry;

    breakSymmetry(write_info);
    double cur_lh = entry_logl;
    {
        vector<double> th;
        map_->pack(th);
        bool moved = false;
        for (size_t i = 0; i < th.size(); i++) if (th[i] != entry.theta[i]) moved = true;
        if (moved) { n_lh_++; cur_lh = tree_->computeLikelihood(); }
    }

    // observable contract of the default loop (design 3.2): branch step per
    // fixed_len, a parameter step, VB_MED progress lines, the iteration cap,
    // and the terminal branch optimisation
    int max_rounds = max(1, params.num_param_iterations - 2);
    for (int k = 1; k <= max_rounds; k++) {
        double new_lh;
        if (fixed_len == BRLEN_OPTIMIZE)
            new_lh = tree_->optimizeAllBranches(min(k + 1, 3), logl_epsilon);
        else if (fixed_len == BRLEN_SCALE) {
            double scaling = 1.0;
            new_lh = tree_->optimizeTreeLengthScaling(MIN_BRLEN_SCALE, scaling, MAX_BRLEN_SCALE, gradient_epsilon);
        } else
            new_lh = cur_lh;
        (void)new_lh;

        new_lh = polish(gradient_epsilon);
        rounds_ = k;

        if (verbose_mode >= VB_MED) {
            tree_->getModel()->writeInfo(cout);
            tree_->getRate()->writeInfo(cout);
            if (fixed_len == BRLEN_SCALE)
                cout << "Scaled tree length: " << tree_->treeLength() << endl;
        }
        if (new_lh > best.logl) snapshot(best, new_lh);
        if (new_lh > cur_lh + logl_epsilon) {
            cur_lh = new_lh;
            if (write_info) {
                ostringstream os;
                os << (k + 1) << ". Current log-likelihood: " << cur_lh;
                if (verbose_mode >= VB_MED) os << " (after " << (getRealTime() - t0) << " wall-clock sec)";
                say(os.str());
            }
        } else {
            cur_lh = max(cur_lh, new_lh);
            break;
        }
    }
    if (fixed_len == BRLEN_OPTIMIZE)
        cur_lh = tree_->optimizeAllBranches(100, logl_epsilon);
    else if (fixed_len == BRLEN_SCALE) {
        double scaling = 1.0;
        cur_lh = tree_->optimizeTreeLengthScaling(MIN_BRLEN_SCALE, scaling, MAX_BRLEN_SCALE, gradient_epsilon);
    }
    if (cur_lh > best.logl) snapshot(best, cur_lh);

    // NOTE(design 6.4): whole-call guard. IQTree::optimizeModelParameters
    // aborts on a regression > 1 below the entry score; per-round acceptance
    // does not bound the composition, so restore the best state seen.
    if (cur_lh < max(entry.logl, best.logl) - 1e-6) {
        const BestState &back = best.logl >= entry.logl ? best : entry;
        restore(back);
        n_lh_++;
        cur_lh = tree_->computeLikelihood();
        ostringstream os;
        os << "NOTE: analytic optimiser ended below its best state (" << setprecision(10) << back.logl << "); restored it";
        say(os.str());
    }
    if (params.ag_stats) {
        ostringstream os;
        os << "AG stats: rounds=" << rounds_ << " bfgs_iterations=" << n_bfgs_iter_ << " likelihood_evaluations=" << n_lh_ << " gradient_evaluations=" << n_grad_
           << " fd_fallbacks=" << n_fd_fallback_ << " parameters=" << map_->ndim()
           << " time=" << setprecision(3) << (getRealTime() - t0) << " sec";
        say(os.str());
    }
    tree_->setCurScore(cur_lh);
    return cur_lh;
}

/* ---------------------------------------------------------------------- */
/* gradient check during optimisation                                      */
/* ---------------------------------------------------------------------- */

void GradientOptimizer::gradientCheckStep(const vector<double> &theta, const vector<double> &g) {
    Params &params = Params::getInstance();
    string prefix = params.out_prefix ? string(params.out_prefix) : string("iqtree");
    double gmax = 0.0;
    for (double v : g) gmax = max(gmax, fabs(v));
    vector<CheckRow> rows;
    CheckStats st;
    paramRows(*map_, tree_, theta, g, params.ag_gradient_check_tol, gmax, rows, st);
    string level = "polish:" + convertIntToString(rounds_ + 1);
    string fname = prefix + ".gradcheck.tsv";
#ifdef _OPENMP
#pragma omp critical(ag_gradcheck_file)
#endif
    {
        bool fresh = !gradcheck_started_;
        gradcheck_started_ = true;
        ofstream out(fname.c_str(), fresh ? ios::out : ios::app);
        if (fresh) out << GRADCHECK_HEADER;
        writeRows(out, rows, n_grad_, level);
    }
    ostringstream os;
    os << scientific << setprecision(2) << "GRADCHECK iter=" << n_grad_ << " level=" << level << " n=" << rows.size()
       << " max_rel=" << st.max_rel << " worst=" << st.worst << " n_fail=" << st.n_fail << " n_fdnoise=" << st.n_noisy;
    say(os.str());
    if (st.n_fail > 0 && params.ag_gradient_check_strict)
        outError("--ag-gradient-check-strict: analytic and numerical gradients disagree (see " + fname + ")");
}

/* ---------------------------------------------------------------------- */
/* self-test and check-only tools                                          */
/* ---------------------------------------------------------------------- */

int GradientOptimizer::selfTest(ModelParamMap *map, PhyloTree *tree) {
    int fails = 0;
    // 1. X kernel: three regimes against a long-double reference; the naive
    //    double formula must be shown to lose accuracy for tiny eigenvalue gaps
    double worst = 0.0, naive_worst = 0.0;
    const double taus[] = { 0.01, 0.5, 3.0 };
    const double lam = -1.3;
    for (double tau : taus) {
        const double gaps[] = { 0.0, 1e-12, 1e-6, 0.05 / tau, 3.0 };
        for (double gap : gaps) {
            double a = lam, b = lam - gap;
            double ref = (double)xRef(a, b, tau);
            double got = PhyloGradient::xKernel(a, b, tau);
            worst = max(worst, fabs(got - ref) / fabs(ref));
            if (gap == 1e-12) {
                double naive = (exp(a * tau) - exp(b * tau)) / (a - b);
                naive_worst = max(naive_worst, fabs(naive - ref) / fabs(ref));
            }
        }
    }
    bool xk_ok = worst <= 1e-13 && naive_worst > 1e-8;
    if (!xk_ok) fails++;
    ios::fmtflags saved_flags = cout.flags();
    cout << scientific << setprecision(3);
    cout << "SELFTEST xkernel max_rel=" << worst
         << " naive_at_1e-12=" << naive_worst << " -> " << (xk_ok ? "PASS" : "FAIL") << endl;
    // 2. pack/unpack round trip
    bool rt_ok = true;
    double rt_worst = 0.0, dl = 0.0;
    if (map && map->unsupported().empty()) {
        tree->clearAllPartialLH();
        double l0 = tree->computeLikelihood();
        vector<double> t0, t1;
        map->pack(t0);
        map->unpack(t0);
        map->pack(t1);
        for (int i = 0; i < map->ndim(); i++) rt_worst = max(rt_worst, fabs(t1[i] - t0[i]));
        double l1 = tree->computeLikelihood();
        dl = l1 - l0;
        rt_ok = rt_worst <= 1e-10;
    }
    if (!rt_ok) fails++;
    cout << "SELFTEST roundtrip ndim=" << (map ? map->ndim() : 0) << " max_abs=" << rt_worst
         << " logl_change=" << dl << " -> " << (rt_ok ? "PASS" : "FAIL") << endl;
    cout.flags(saved_flags);
    return fails ? 1 : 0;
}

int GradientOptimizer::gradientCheckOnly(ModelFactory *factory, PhyloTree *tree) {
    Params &params = Params::getInstance();
    const double tol = params.ag_gradient_check_tol;
    string prefix = params.out_prefix ? string(params.out_prefix) : string("iqtree");
    ModelSubst *model = tree->getModel();
    RateHeterogeneity *rate = tree->getRate();
    const int nst = tree->aln->num_states;
    const int nmix = model->getNMixtures();
    const int ncat = rate->getNRate();

    // 0. parameter map (model parameters); branch lengths are checked regardless
    ModelParamMap map(factory, tree);
    const bool have_map = map.unsupported().empty();
    if (!have_map)
        cout << "AG: model parameters not representable (" << map.unsupported() << "); checking branch lengths only" << endl;
    int rc_self = 0;
    if (params.ag_selftest) rc_self = selfTest(have_map ? &map : nullptr, tree);
    vector<double> theta;
    if (have_map) {
        map.pack(theta);
        map.unpack(theta);   // put the model exactly on the parametrised manifold (mean rate 1, sum pi 1)
    }

    // 1. analytic gradient
    PhyloGradient engine(tree);
    if (have_map) engine.setNeedQ(map.needQ());
    PhyloGradient::Result res;
    bool ok = engine.compute(res);
    cout << "AG: analytic gradient at initial point: logl=" << setprecision(10) << res.logl
         << ", edges=" << res.num_edges
         << ", max|lnL(edge)-lnL(root)|=" << scientific << setprecision(2) << res.max_edge_logl_diff
         << ", buffers=" << fixed << setprecision(3) << (engine.bufferBytes() / (1024.0 * 1024.0)) << " MB"
         << (ok ? "" : ", NON-FINITE VALUES") << endl;

    vector<double> g, g_fdq;
    double qchain = 0.0;
    if (have_map && map.ndim() > 0) {
        map.gradient(res, g);
        map.gradientFDQ(res, g_fdq);
        double gmax = 0.0;
        for (double v : g) gmax = max(gmax, fabs(v));
        for (int i = 0; i < map.ndim(); i++) qchain = max(qchain, relErr(g[i], g_fdq[i], gmax));
    }

    // 2. references per branch: Newton derivative and central differences of the tree lnL
    BranchVector branches;
    tree->getBranches(branches);
    double gmax = 0.0;
    for (double v : res.dlogl_dt) gmax = max(gmax, fabs(v));
    for (double v : g) gmax = max(gmax, fabs(v));

    vector<CheckRow> rows;
    CheckStats st;
    double sum_t_dt = 0.0;

    for (auto &br : branches) {
        Node *n1 = br.first, *n2 = br.second;
        PhyloNeighbor *nei = (PhyloNeighbor*)n1->findNeighbor(n2);
        PhyloNeighbor *rev = (PhyloNeighbor*)n2->findNeighbor(n1);
        CheckRow r;
        int id = nei->id;
        r.idx = convertIntToString(id);
        r.name = "t[" + r.idx + "]";
        r.label = sideLabel(n2, n1);
        r.x = nei->length;
        r.analytic = (id >= 0 && id < (int)res.dlogl_dt.size()) ? res.dlogl_dt[id] : NAN;
        sum_t_dt += r.x * r.analytic;

        // Newton reference: the kernel caches its edge buffer (theta_all) between
        // calls, so the flag must be reset for every branch or stale data is reused.
        double df = 0.0, ddf = 0.0;
        tree->theta_computed = false;
        tree->computeLikelihoodDerv(nei, (PhyloNode*)n1, &df, &ddf);
        tree->theta_computed = false;
        r.newton = df;
        r.has_newton = true;

        // central differences (two step sizes) of the full log-likelihood; the
        // step must stay inside (0, t) so a branch at the minimum length is not
        // differenced one-sidedly (that halves the estimate)
        double h = min(1e-4 * max(1.0, r.x), 0.5 * r.x);
        double lh[4];
        double steps[4] = { h, -h, h / 2, -h / 2 };
        for (int k = 0; k < 4; k++) {
            nei->length = rev->length = r.x + steps[k];
            tree->clearAllPartialLH();
            lh[k] = tree->computeLikelihood();
        }
        nei->length = rev->length = r.x;
        r.fd_h = (lh[0] - lh[1]) / (2 * h);
        r.fd_h2 = (lh[2] - lh[3]) / h;
        r.rich = (4 * r.fd_h2 - r.fd_h) / 3;

        // "numeric" = Richardson-extrapolated central difference (the user-facing
        // comparison); the Newton derivative is reported as an extra column.
        judgeRow(r, tol, gmax, st);
        rows.push_back(r);
    }

    // 3. model parameters: central differences in theta space through unpack()
    if (have_map) paramRows(map, tree, theta, g, tol, gmax, rows, st);
    // leave the tree in a consistent state
    tree->clearAllPartialLH();
    tree->computeLikelihood();

    // 4. zero-cost identities (design 5)
    //    sum_m <D_m, Q_m> = sum_e t_e dlogL/dt_e  (scaling Q == scaling every branch)
    //    sum_m sum_k pi_mk root_term[m][k] = sum_c omega_c dlogL/domega_c
    double q_id = have_map ? map.qIdentity(res) : 0.0;
    bool any_q = false;
    for (auto &d : res.dlogl_dQ) any_q |= !d.empty();
    double q_rel = any_q ? fabs(q_id - sum_t_dt) / max(max(fabs(q_id), fabs(sum_t_dt)), 1e-6 * max(gmax, 1.0)) : 0.0;
    double root_lhs = 0.0, root_rhs = 0.0;
    {
        vector<double> pi(nst);
        for (int m = 0; m < nmix; m++) {
            model->getStateFrequency(pi.data(), m);
            for (int k = 0; k < nst; k++) root_lhs += pi[k] * res.root_term[m][k];
            for (int r = 0; r < ncat; r++)
                root_rhs += rate->getProp(r) * model->getMixtureWeight(m) * res.dlogl_dclass[m * ncat + r];
        }
    }
    double root_rel = fabs(root_lhs - root_rhs) / max(max(fabs(root_lhs), fabs(root_rhs)), 1e-6 * max(gmax, 1.0));
    bool id_ok = q_rel <= 1e-8 && root_rel <= 1e-8 && qchain <= 1e-5;

    // 5. outputs
    {
        ofstream out((prefix + ".gradcheck.tsv").c_str());
        out << GRADCHECK_HEADER;
        writeRows(out, rows, 0, "init");
    }
    if (params.ag_dump_gradient) {
        ofstream out((prefix + ".aggrad.tsv").c_str());
        out << setprecision(17);
        out << "#logl\t" << res.logl << '\n';
        out << "#tree\t";
        tree->printTree(out, WT_BR_LEN);
        out << '\n';
        // model parameters, enough for an external oracle to rebuild Q per component
        out << "#seqtype\t" << (tree->aln->seq_type == SEQ_DNA ? "DNA" : tree->aln->seq_type == SEQ_PROTEIN ? "AA" : "OTHER") << '\n';
        out << "#nstates\t" << nst << '\n';
        out << "#nmix\t" << nmix << '\n';
        for (int m = 0; m < nmix; m++) {
            ModelSubst *comp = model->isMixture() ? model->getMixtureClass(m) : model;
            double weight = model->isMixture() ? model->getMixtureWeight(m) : 1.0;
            vector<double> rates(comp->getNumRateEntries(), 0.0), freqs(nst, 0.0);
            comp->getRateMatrix(rates.data());
            comp->getStateFrequency(freqs.data());
            out << "#mix\t" << m << "\tweight\t" << weight << "\treversible\t" << (comp->isReversible() ? 1 : 0) << '\n';
            out << "#exchangeabilities\t" << m;
            for (double x : rates) out << '\t' << x;
            out << '\n';
            out << "#freqs\t" << m;
            for (double x : freqs) out << '\t' << x;
            out << '\n';
            if (!res.dlogl_dQ[m].empty()) {
                out << "#dlogl_dQ\t" << m;
                for (double x : res.dlogl_dQ[m]) out << '\t' << x;
                out << '\n';
            }
        }
        out << "#pinvar\t" << rate->getPInvar() << '\n';
        out << "#ncat\t" << ncat << '\n';
        out << "#cat_rates";
        for (int c = 0; c < ncat; c++) out << '\t' << rate->getRate(c);
        out << '\n';
        out << "#cat_props";
        for (int c = 0; c < ncat; c++) out << '\t' << rate->getProp(c);
        out << '\n';
        out << "#fused\t" << (factory->fused_mix_rate ? 1 : 0) << '\n';
        // dlogL/dQ is the derivative of the likelihood rooted at this side's node
        out << "#root_side\t" << sideLabel(res.root_dad, res.root_node) << '\n';
        if (have_map && map.ndim() > 0) {
            const ModelParamMap::NaturalGradient &nat = map.natural();
            for (int m = 0; m < nmix; m++) {
                if (!map.needQ()[m]) continue;   // fixed matrix and frequencies: no Q gradient computed
                out << "#nat_dR\t" << m;
                for (double x : nat.dR[m]) out << '\t' << x;
                out << '\n';
                out << "#nat_dpi\t" << m;
                for (double x : nat.dpi[m]) out << '\t' << x;
                out << '\n';
            }
            out << "#nat_dw";
            for (double x : nat.dw) out << '\t' << x;
            out << '\n';
            out << "#nat_drate";
            for (double x : nat.drate) out << '\t' << x;
            out << '\n';
            out << "#nat_dprop";
            for (double x : nat.dprop) out << '\t' << x;
            out << '\n';
            out << "#nat_dpinv_direct\t" << nat.dpinv_direct << '\n';
            out << "#nat_dpinv\t" << nat.dpinv << '\n';
            out << "#nat_dalpha\t" << nat.dalpha << '\n';
            out << "#theta";
            for (int i = 0; i < map.ndim(); i++) out << '\t' << map.name(i) << '=' << theta[i] << ':' << g[i];
            out << '\n';
        }
        out << "id\tlength\tdlogl_dt\tside_taxa\n";
        for (auto &r : rows)
            if (r.has_newton)
                out << r.idx << '\t' << setprecision(17) << r.x << '\t' << r.analytic << '\t' << r.label << '\n';
        cout << "AG: raw gradient written to " << prefix << ".aggrad.tsv" << endl;
    }
    bool edge_ok = res.max_edge_logl_diff <= 1e-6 * max(1.0, fabs(res.logl));
    ios::fmtflags saved_flags = cout.flags();
    cout << scientific << setprecision(2);
    cout << "GRADCHECK iter=0 n=" << rows.size() << " max_rel=" << st.max_rel
         << " worst=" << st.worst << " n_fail=" << st.n_fail << " n_fdnoise=" << st.n_noisy
         << " edge_lnl_check=" << (edge_ok ? "PASS" : "FAIL")
         << " qchain=" << qchain << " q_identity=" << q_rel << " root_identity=" << root_rel
         << " identities=" << (id_ok ? "PASS" : "FAIL")
         << " (" << prefix << ".gradcheck.tsv)" << endl;
    cout.flags(saved_flags);
    return (st.n_fail == 0 && ok && edge_ok && id_ok && rc_self == 0) ? 0 : 1;
}
