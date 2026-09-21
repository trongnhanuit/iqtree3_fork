/*
 * gradientoptimizer.cpp
 *
 * See gradientoptimizer.h and docs/analytical-gradients-design.md.
 */

#include "gradientoptimizer.h"
#include "modelfactory.h"
#include "tree/phylotree.h"
#include "phylogradient.h"
#include "utils/tools.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>

using namespace std;

bool GradientOptimizer::supports(ModelFactory *factory, std::string &why) {
    (void)factory;
    // Stage 0/1 scaffolding: the pipeline is not available yet, so every model
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

} // namespace

int GradientOptimizer::gradientCheckOnly(ModelFactory *factory, PhyloTree *tree) {
    (void)factory;
    Params &params = Params::getInstance();
    const double tol = params.ag_gradient_check_tol;
    string prefix = params.out_prefix ? string(params.out_prefix) : string("iqtree");

    // 1. analytic gradient
    PhyloGradient engine(tree);
    PhyloGradient::Result res;
    bool ok = engine.compute(res);
    cout << "AG: analytic gradient at initial point: logl=" << setprecision(10) << res.logl
         << ", edges=" << res.num_edges
         << ", max|lnL(edge)-lnL(root)|=" << setprecision(3) << res.max_edge_logl_diff
         << ", buffers=" << (engine.bufferBytes() / (1024.0 * 1024.0)) << " MB"
         << (ok ? "" : ", NON-FINITE VALUES") << endl;

    // 2. references per branch: Newton derivative and central differences of the tree lnL
    BranchVector branches;
    tree->getBranches(branches);
    double gmax = 0.0;
    for (double g : res.dlogl_dt) gmax = max(gmax, fabs(g));

    struct Row { int id; string label; double t, analytic, newton, fd_h, fd_h2, rich, rel; bool pass; };
    vector<Row> rows;
    int n_fail = 0; double max_rel = 0.0; string worst;

    for (auto &br : branches) {
        Node *n1 = br.first, *n2 = br.second;
        PhyloNeighbor *nei = (PhyloNeighbor*)n1->findNeighbor(n2);
        PhyloNeighbor *rev = (PhyloNeighbor*)n2->findNeighbor(n1);
        Row r;
        r.id = nei->id;
        r.label = sideLabel(n2, n1);
        r.t = nei->length;
        r.analytic = (r.id >= 0 && r.id < (int)res.dlogl_dt.size()) ? res.dlogl_dt[r.id] : NAN;

        // Newton reference: the kernel caches its edge buffer (theta_all) between
        // calls, so the flag must be reset for every branch or stale data is reused.
        double df = 0.0, ddf = 0.0;
        tree->theta_computed = false;
        tree->computeLikelihoodDerv(nei, (PhyloNode*)n1, &df, &ddf);
        tree->theta_computed = false;
        r.newton = df;

        // central differences (two step sizes) of the full log-likelihood; the
        // step must stay inside (0, t) so a branch at the minimum length is not
        // differenced one-sidedly (that halves the estimate)
        double h = min(1e-4 * max(1.0, r.t), 0.5 * r.t);
        double lh[4];
        double steps[4] = { h, -h, h / 2, -h / 2 };
        for (int k = 0; k < 4; k++) {
            double tk = r.t + steps[k];
            nei->length = rev->length = tk;
            tree->clearAllPartialLH();
            lh[k] = tree->computeLikelihood();
        }
        nei->length = rev->length = r.t;
        r.fd_h = (lh[0] - lh[1]) / (2 * h);
        r.fd_h2 = (lh[2] - lh[3]) / h;
        r.rich = (4 * r.fd_h2 - r.fd_h) / 3;

        // "numeric" = Richardson-extrapolated central difference (the user-facing
        // comparison); the Newton derivative is reported as an extra column.
        r.rel = relErr(r.analytic, r.rich, gmax);
        r.pass = std::isfinite(r.analytic) && r.rel <= tol;
        if (!r.pass) n_fail++;
        if (r.rel > max_rel) { max_rel = r.rel; worst = "t[" + convertIntToString(r.id) + "]"; }
        rows.push_back(r);
    }
    // leave the tree in a consistent state
    tree->clearAllPartialLH();
    tree->computeLikelihood();

    // 3. outputs
    {
        ofstream out((prefix + ".gradcheck.tsv").c_str());
        out << "iter\tlevel\tidx\tname\tx\tanalytic\tnumeric\tfd_h\tfd_h2\tnewton\tabs_err\trel_err\tstatus\n";
        out << setprecision(10);
        for (auto &r : rows)
            out << 0 << "\tinit\t" << r.id << "\tt[" << r.id << "]\t" << r.t << '\t' << r.analytic << '\t'
                << r.rich << '\t' << r.fd_h << '\t' << r.fd_h2 << '\t' << r.newton << '\t'
                << fabs(r.analytic - r.rich) << '\t' << r.rel << '\t' << (r.pass ? "PASS" : "FAIL") << '\n';
    }
    if (params.ag_dump_gradient) {
        ofstream out((prefix + ".aggrad.tsv").c_str());
        out << setprecision(17);
        out << "#logl\t" << res.logl << '\n';
        out << "#tree\t";
        tree->printTree(out, WT_BR_LEN);
        out << '\n';
        // model parameters, enough for an external oracle to rebuild Q per component
        ModelSubst *model = tree->getModel();
        RateHeterogeneity *rate = tree->getRate();
        int nst = tree->aln->num_states;
        int nmix = model->getNMixtures();
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
        }
        out << "#pinvar\t" << rate->getPInvar() << '\n';
        out << "#ncat\t" << rate->getNRate() << '\n';
        out << "#cat_rates";
        for (int c = 0; c < rate->getNRate(); c++) out << '\t' << rate->getRate(c);
        out << '\n';
        out << "#cat_props";
        for (int c = 0; c < rate->getNRate(); c++) out << '\t' << rate->getProp(c);
        out << '\n';
        out << "#fused\t" << (factory->fused_mix_rate ? 1 : 0) << '\n';
        out << "id\tlength\tdlogl_dt\tside_taxa\n";
        for (auto &r : rows)
            out << r.id << '\t' << setprecision(17) << r.t << '\t' << r.analytic << '\t' << r.label << '\n';
        cout << "AG: raw gradient written to " << prefix << ".aggrad.tsv" << endl;
    }
    cout << "GRADCHECK iter=0 n=" << rows.size() << " max_rel=" << setprecision(3) << max_rel
         << " worst=" << worst << " n_fail=" << n_fail
         << " edge_lnl_check=" << (res.max_edge_logl_diff <= 1e-6 * max(1.0, fabs(res.logl)) ? "PASS" : "FAIL")
         << " (" << prefix << ".gradcheck.tsv)" << endl;
    bool edge_ok = res.max_edge_logl_diff <= 1e-6 * max(1.0, fabs(res.logl));
    return (n_fail == 0 && ok && edge_ok) ? 0 : 1;
}
