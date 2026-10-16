/*
 * agtrace.cpp
 * --ag-trace writer, see agtrace.h. Record types (field "type"):
 *   header  taxa, model, states (once per run segment; "resume" if appended)
 *   round   one GradientOptimizer round: call id, kind (start/em/polish), lnl, params
 *   call    end of a model-optimisation call: call id, iteration count, lnl, tree, params
 *   phase   "search" when the tree search starts
 *   iter    one tree-search iteration: it, lnl, tree (taxon ids, branch lengths)
 *   end     normal close
 */

#include "agtrace.h"
#include "tree/phylotree.h"
#include "utils/tools.h"
#include <cstdio>

namespace AgTrace {

bool on = false;
static FILE *out = nullptr;
static PhyloTree *traced = nullptr;
static int ncall = 0;

static void num(string &s, double x) {
    if (x != x || x > 1e300 || x < -1e300) { s += "null"; return; }
    char buf[32];
    snprintf(buf, sizeof(buf), "%.10g", x);
    s += buf;
}

static void str(string &s, const string &x) {
    s += '"';
    for (char c : x) {
        if (c == '"' || c == '\\') s += '\\';
        if (c != '\n' && c != '\r') s += c;
    }
    s += '"';
}

static void arr(string &s, const vector<double> &v) {
    s += '[';
    for (size_t i = 0; i < v.size(); i++) { if (i) s += ','; num(s, v[i]); }
    s += ']';
}

static void line(const string &s) {
    fputs(s.c_str(), out);
    fputc('\n', out);
    fflush(out);
}

// model parameters through plain getters; rates of a component equal to
// component 0's are written as null to keep mixture snapshots small
static string params(PhyloTree *tree) {
    ModelSubst *model = tree->getModel();
    RateHeterogeneity *rate = tree->getRate();
    int nmix = model->isMixture() ? model->getNMixtures() : 1;
    int nst = tree->aln->num_states;
    string s = "{\"tl\":";
    num(s, tree->treeLength());
    vector<double> w(nmix, 1.0), rates0;
    string rs = "[", fs = "[";
    for (int m = 0; m < nmix; m++) {
        ModelSubst *comp = model->isMixture() ? model->getMixtureClass(m) : model;
        if (model->isMixture()) w[m] = model->getMixtureWeight(m);
        vector<double> r(comp->getNumRateEntries(), 0.0), f(nst, 0.0);
        comp->getRateMatrix(r.data());
        comp->getStateFrequency(f.data());
        if (m) { rs += ','; fs += ','; }
        if (m == 0) rates0 = r;
        if (m > 0 && r == rates0) rs += "null"; else arr(rs, r);
        arr(fs, f);
    }
    s += ",\"w\":"; arr(s, w);
    s += ",\"rates\":" + rs + "],\"freqs\":" + fs + "]";
    s += ",\"alpha\":"; num(s, rate->getGammaShape());
    s += ",\"pinv\":"; num(s, rate->getPInvar());
    vector<double> cr, cp;
    for (int c = 0; c < rate->getNRate(); c++) { cr.push_back(rate->getRate(c)); cp.push_back(rate->getProp(c)); }
    s += ",\"cat_rates\":"; arr(s, cr);
    s += ",\"cat_props\":"; arr(s, cp);
    return s + "}";
}

void open(PhyloTree *tree) {
    Params &p = Params::getInstance();
    if (!p.ag_trace) return;
    if (tree->isSuperTree()) {
        outWarning("--ag-trace does not support partitioned models; no trace written");
        return;
    }
    string fn = string(p.out_prefix) + ".agtrace.jsonl";
    bool resume = !p.ignore_checkpoint && fileExists(fn);
    out = fopen(fn.c_str(), resume ? "ab" : "wb");
    if (!out) {
        outWarning("--ag-trace: cannot write " + fn);
        return;
    }
    traced = tree;
    ncall = 0;
    on = true;
    string s = "{\"type\":\"header\",\"resume\":";
    s += resume ? "true" : "false";
    s += ",\"model\":"; str(s, p.model_name);
    s += ",\"ag\":"; s += p.analytical_gradients ? "true" : "false";
    s += ",\"seq_type\":"; str(s, tree->aln->seq_type == SEQ_DNA ? "DNA" : tree->aln->seq_type == SEQ_PROTEIN ? "AA" : "OTHER");
    s += ",\"nstates\":" + convertIntToString(tree->aln->num_states);
    s += ",\"taxa\":[";
    for (int i = 0; i < tree->aln->getNSeq(); i++) { if (i) s += ','; str(s, tree->aln->getSeqName(i)); }
    line(s + "]}");
    cout << "NOTE: writing optimisation trace to " << fn << endl;
}

void close() {
    if (!on) return;
    line("{\"type\":\"end\"}");
    fclose(out);
    out = nullptr;
    traced = nullptr;
    on = false;
}

void phase(PhyloTree *tree, const char *name) {
    if (tree != traced) return;
    line(string("{\"type\":\"phase\",\"name\":\"") + name + "\"}");
}

void iteration(PhyloTree *tree, int it, double lnl, const string &tree_str) {
    if (tree != traced) return;
    string s = "{\"type\":\"iter\",\"it\":" + convertIntToString(it) + ",\"lnl\":";
    num(s, lnl);
    s += ",\"tree\":"; str(s, tree_str);
    line(s + "}");
}

void callEnd(PhyloTree *tree, int it, double lnl, const string &tree_str) {
    if (tree != traced) return;
    string s = "{\"type\":\"call\",\"call\":" + convertIntToString(ncall++) + ",\"it\":" + convertIntToString(it) + ",\"lnl\":";
    num(s, lnl);
    s += ",\"tree\":"; str(s, tree_str);
    s += ",\"params\":" + params(tree);
    line(s + "}");
}

void round(PhyloTree *tree, const char *kind, int round, double lnl) {
    if (tree != traced) return;
    string s = "{\"type\":\"round\",\"call\":" + convertIntToString(ncall) + ",\"kind\":\"" + kind
        + "\",\"round\":" + convertIntToString(round) + ",\"lnl\":";
    num(s, lnl);
    s += ",\"params\":" + params(tree);
    line(s + "}");
}

}
