/*
 * gradientoptimizer.cpp
 *
 * See gradientoptimizer.h and docs/analytical-gradients-design.md.
 */

#include "gradientoptimizer.h"
#include "modelfactory.h"
#include "tree/phylotree.h"
#include "phylogradient.h"
#include "modelparammap.h"
#include "utils/tools.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>

using namespace std;

bool GradientOptimizer::supports(ModelFactory *factory, std::string &why) {
    (void)factory;
    // Stage 0-2 scaffolding: the pipeline is not available yet, so every model
    // falls back to the default optimiser. A later commit replaces this with the
    // real capability check (design doc, section 4).
    why = "analytical-gradient pipeline not yet available";
    return false;
}

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

} // namespace

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

    struct Row { string idx, name, label; double x, analytic, newton, fd_h, fd_h2, rich, rel; bool has_newton, pass, noisy; };
    vector<Row> rows;
    int n_fail = 0, n_noisy = 0; double max_rel = 0.0; string worst;
    // A row whose two central-difference estimates disagree with each other by
    // more than they disagree with the analytic value is limited by rounding
    // in the likelihood (tiny gradients after a 20x20 eigendecomposition); it
    // is counted as noisy, not as a failure.
    auto judge = [&](Row &r) {
        r.rel = relErr(r.analytic, r.rich, gmax);
        double fd_noise = fabs(r.fd_h - r.fd_h2);
        r.noisy = false;
        r.pass = std::isfinite(r.analytic) && r.rel <= tol;
        if (!r.pass && std::isfinite(r.analytic) && fabs(r.analytic - r.rich) <= 4.0 * fd_noise) { r.pass = true; r.noisy = true; n_noisy++; }
        if (!r.pass) n_fail++;
        if (r.rel > max_rel && !r.noisy) { max_rel = r.rel; worst = r.name; }
    };
    double sum_t_dt = 0.0;

    for (auto &br : branches) {
        Node *n1 = br.first, *n2 = br.second;
        PhyloNeighbor *nei = (PhyloNeighbor*)n1->findNeighbor(n2);
        PhyloNeighbor *rev = (PhyloNeighbor*)n2->findNeighbor(n1);
        Row r;
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
        judge(r);
        rows.push_back(r);
    }

    // 3. model parameters: central differences in theta space through unpack()
    for (int i = 0; have_map && i < map.ndim(); i++) {
        Row r;
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
        judge(r);
        rows.push_back(r);
    }
    // leave the tree in a consistent state
    if (have_map) map.unpack(theta);
    tree->clearAllPartialLH();
    tree->computeLikelihood();

    // 4. zero-cost identities (design doc 6.2)
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
        out << "iter\tlevel\tidx\tname\tx\tanalytic\tnumeric\tfd_h\tfd_h2\tnewton\tabs_err\trel_err\tstatus\n";
        out << setprecision(10);
        for (auto &r : rows) {
            out << 0 << "\tinit\t" << r.idx << '\t' << r.name << '\t' << r.x << '\t'
                << r.analytic << '\t' << r.rich << '\t' << r.fd_h << '\t' << r.fd_h2 << '\t';
            if (r.has_newton) out << r.newton; else out << "NA";
            out << '\t' << fabs(r.analytic - r.rich) << '\t' << r.rel << '\t' << (r.noisy ? "PASS-FDNOISE" : r.pass ? "PASS" : "FAIL") << '\n';
        }
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
    cout << "GRADCHECK iter=0 n=" << rows.size() << " max_rel=" << max_rel
         << " worst=" << worst << " n_fail=" << n_fail << " n_fdnoise=" << n_noisy
         << " edge_lnl_check=" << (edge_ok ? "PASS" : "FAIL")
         << " qchain=" << qchain << " q_identity=" << q_rel << " root_identity=" << root_rel
         << " identities=" << (id_ok ? "PASS" : "FAIL")
         << " (" << prefix << ".gradcheck.tsv)" << endl;
    cout.flags(saved_flags);
    return (n_fail == 0 && ok && edge_ok && id_ok && rc_self == 0) ? 0 : 1;
}
