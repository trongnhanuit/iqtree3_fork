#!/usr/bin/env python3
"""
Benchmark reports for the C10 re-estimation (layout written by the batch script):

  report.py --scores SCORES --runs RUNS --c10pf C10PF.nex --out DIR

SCORES/<S1|S2|S3>/<tb|hog>_<model>/scores.tsv (+ sitelh/), models C10, C10pf,
warm_sel, cold_sel (all settings), WAGHSSP (S1), C20pf (S3) and the other
finalists <arm>_<chain> (S1, TreeBase). RUNS/<warm|cold>/summary.tsv,
finalists/, halving.json.

DIR/report1.md  main result: the selected warm and cold models against C10pf
                (primary baseline) and built-in C10, three settings, two test
                sets, and the confirmatory Holm-corrected Wilcoxon tests (S1, TreeBase).
DIR/report2.md  extras: AIC/BIC per site against WAG_HSSP, KH counts against
                WAG_HSSP, the C10->C20 gap closed, the lucky-optimum check,
                profile matching against C10 and training diagnostics.
"""
import os
import sys
import json
import argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from evaluate import summarize
from common import read_mdef, match_profiles, AA
import numpy as np

SETS = {"S1": "S1 (primary): topology searched per model, weights fixed",
        "S2": "S2: topology searched per model, weights re-optimised per alignment",
        "S3": "S3: topology fixed to the supplied tree, weights fixed"}
TESTS = {"tb": "TreeBase-57", "hog": "HOGENOM-300"}
NAMES = {"warm_sel": "warm", "cold_sel": "cold", "C10pf": "C10pf", "C10": "built-in C10"}


def top3(runs, arm):
    rows = [l.rstrip("\n").split("\t") for l in list(open(os.path.join(runs, arm, "summary.tsv")))[1:]]
    rows = [r for r in rows if r[2] in ("converged", "not converged")]
    rows.sort(key=lambda r: -float(r[4]))
    return rows[:3], rows


def sc(scores, s, ts, m):
    return os.path.join(scores, s, "%s_%s" % (ts, m), "scores.tsv")


def row(label, r):
    return ("| %s | %d | %+.2f [%+.2f, %+.2f] | %d | %d / %d | %d / %d | %d / %d | %+.4f [%+.4f, %+.4f] | %.2g |"
            % (label, r["n"], r["dlnl_median"], r["dlnl_q1"], r["dlnl_q3"], r["dlnl_pos"], r["daic10_new"],
               r["daic10_base"], r["kh01"][0], r["kh01"][1], r["kh05"][0], r["kh05"][1], r["aic_site_mean"],
               r["aic_site_ci"][0], r["aic_site_ci"][1], r["wilcoxon_p"]))


HEAD = ("| comparison | n | dlnL median [IQR] | dlnL > 0 | dAIC>10 new / base | KH p<0.01 better / worse "
        "| KH p<0.05 better / worse | mean dAIC/site [95% CI] | Wilcoxon p |\n|---|---|---|---|---|---|---|---|---|")


def report1(a, sel):
    L = ["# Report 1: re-estimated C10 vs the published C10", "",
         "Selected models (best penalised training lnL per arm): warm = %s, cold = %s." % (sel["warm"], sel["cold"]),
         "Baselines: C10pf (the paper's C10: published profiles and weights, paper normalisation; primary) and "
         "IQ-TREE's built-in C10. Every model has the same number of free parameters within a setting, so "
         "dAIC = dBIC = 2 dlnL. KH: two-sided RELL test per alignment on site log-likelihoods.", ""]
    conf = {}
    for s in ("S1", "S2", "S3"):
        L.append("## " + SETS[s])
        for ts in ("tb", "hog"):
            L += ["", "### " + TESTS[ts], "", HEAD]
            for new in ("cold_sel", "warm_sel"):
                for base in ("C10pf", "C10"):
                    r = summarize(sc(a.scores, s, ts, new), sc(a.scores, s, ts, base))
                    L.append(row("%s vs %s" % (NAMES[new], NAMES[base]), r))
                    if s == "S1" and ts == "tb" and base == "C10pf":
                        conf[NAMES[new]] = r
        L.append("")
    L += ["## Confirmatory tests (S1, TreeBase-57, vs C10pf; one-sided Wilcoxon, Holm at familywise 0.01)", ""]
    order = sorted(conf, key=lambda k: conf[k]["wilcoxon_p"])
    holds = True
    for k, name in enumerate(order):
        thr = 0.01 / (len(order) - k)
        ok = holds and conf[name]["wilcoxon_p"] < thr and conf[name]["aic_site_mean"] > 0
        holds = ok
        L.append("- %s%s: p = %.3g (threshold %.4f), mean dAIC/site %+.4f -> %s"
                 % (name, " (primary)" if name == "cold" else "", conf[name]["wilcoxon_p"], thr,
                    conf[name]["aic_site_mean"], "SUPPORTED" if ok else "not supported"))
    return "\n".join(L) + "\n"


def report2(a, sel, fin, allrows):
    L = ["# Report 2: extras", ""]
    # 1. paper scale against WAG_HSSP
    L += ["## AIC and BIC per site relative to WAG_HSSP (the paper's reference)", "",
          "S1: all models' weights fixed. S2: the C10-type models re-optimise 9 weights (AICc and BIC penalties); "
          "WAG_HSSP is its S1 run. KH counts in S1 at p < 0.01 (the paper's Table 1 threshold).", "",
          "| test set | model | S1 mean dAIC/site | S1 mean dBIC/site | S1 KH p<0.01 better / worse | S2 mean dAIC/site | S2 mean dBIC/site |",
          "|---|---|---|---|---|---|---|"]
    for ts in ("tb", "hog"):
        wag = sc(a.scores, "S1", ts, "WAGHSSP")
        for m in ("C10", "C10pf", "warm_sel", "cold_sel"):
            r1 = summarize(sc(a.scores, "S1", ts, m), wag)
            r2 = summarize(sc(a.scores, "S2", ts, m), wag, df_new=9)
            L.append("| %s | %s | %+.4f | %+.4f | %d / %d | %+.4f | %+.4f |"
                     % (TESTS[ts], NAMES[m], r1["aic_site_mean"], r1["bic_site_mean"], r1["kh01"][0], r1["kh01"][1],
                        r2["aic_site_mean"], r2["bic_site_mean"]))
    # 2. C10 -> C20 gap (S3)
    L += ["", "## Fraction of the C10pf -> C20pf gap closed (S3, mean dAIC/site)", "",
          "| test set | C20pf vs C10pf | warm vs C10pf | cold vs C10pf | warm closes | cold closes |", "|---|---|---|---|---|---|"]
    for ts in ("tb", "hog"):
        g = summarize(sc(a.scores, "S3", ts, "C20pf"), sc(a.scores, "S3", ts, "C10pf"), boot=200)["aic_site_mean"]
        w = summarize(sc(a.scores, "S3", ts, "warm_sel"), sc(a.scores, "S3", ts, "C10pf"), boot=200)["aic_site_mean"]
        c = summarize(sc(a.scores, "S3", ts, "cold_sel"), sc(a.scores, "S3", ts, "C10pf"), boot=200)["aic_site_mean"]
        L.append("| %s | %+.4f | %+.4f | %+.4f | %.0f%% | %.0f%% |" % (TESTS[ts], g, w, c, 100 * w / g, 100 * c / g))
    # 3. lucky optimum
    L += ["", "## Lucky-optimum check (S1, TreeBase-57, vs C10pf)", "",
          "The selected model and the median finalist must both clear F = max(0.008, range of the finalists' mean dAIC/site).", ""]
    for arm in ("warm", "cold"):
        means = []
        for k, r in enumerate(fin[arm]):
            m = "%s_sel" % arm if k == 0 else "%s_%s" % (arm, r[0])
            means.append((r[0], summarize(sc(a.scores, "S1", "tb", m), sc(a.scores, "S1", "tb", "C10pf"), boot=200)["aic_site_mean"]))
        vals = [v for _, v in means]
        F = max(0.008, max(vals) - min(vals)) if len(vals) > 1 else 0.008
        med = sorted(vals)[len(vals) // 2]
        ok = vals[0] >= F and med >= F
        L.append("- %s: %s; F = %.4f; selected %+.4f, median %+.4f -> %s"
                 % (arm, ", ".join("%s %+.4f" % mv for mv in means), F, vals[0], med, "PASS" if ok else "not met"))
    # 4. profile matching
    L += ["", "## Profile matching of the selected models against the published C10", ""]
    _, pub, wpub, _ = read_mdef(a.c10pf)
    wpub = wpub / wpub.sum()
    for arm in ("warm", "cold"):
        _, p, w, rates = read_mdef(os.path.join(a.runs, arm, "finalists", sel[arm] + ".nex"))
        w = w / w.sum()
        pairs, mnn = match_profiles(p, pub)
        L += ["### %s (%s): mutual nearest neighbours %d/10" % (arm, sel[arm], sum(mnn)), "",
              "| new class | weight | top amino acids | C10 class | C10 weight | C10 top amino acids | sym-KL |", "|---|---|---|---|---|---|---|"]
        top = lambda q: " ".join("%s %.2f" % (AA[k], q[k]) for k in np.argsort(-q)[:4])
        for (i, j, c), m in zip(pairs, mnn):
            L.append("| %d | %.3f | %s | %d | %.3f | %s | %.3f%s |" % (i + 1, w[i], top(p[i]), j + 1, wpub[j], top(pub[j]), c, "" if m else " (not MNN)"))
        L.append("")
    # 5. training diagnostics
    L += ["## Training diagnostics", ""]
    for arm in ("warm", "cold"):
        L += ["### %s" % arm, "", "| chain | start | status | rounds | penalised training lnL |", "|---|---|---|---|---|"]
        for r in sorted(allrows[arm], key=lambda r: -float(r[4]) if r[4] != "nan" else 1e300):
            L.append("| %s | %s | %s | %s | %s |" % tuple(r))
        h = os.path.join(a.runs, arm, "halving.json")
        if os.path.exists(h):
            plan = json.load(open(h))
            L += ["", "Halving: " + "; ".join("round %s kept %s" % (k, ", ".join(v["kept"])) for k, v in plan.items())]
        L.append("")
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", required=True)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--c10pf", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    fin, allrows, sel = {}, {}, {}
    for arm in ("warm", "cold"):
        fin[arm], allrows[arm] = top3(a.runs, arm)
        sel[arm] = fin[arm][0][0]
    with open(os.path.join(a.out, "report1.md"), "w") as f:
        f.write(report1(a, sel))
    with open(os.path.join(a.out, "report2.md"), "w") as f:
        f.write(report2(a, sel, fin, allrows))
    print("wrote %s/report1.md and report2.md" % a.out)


if __name__ == "__main__":
    main()
