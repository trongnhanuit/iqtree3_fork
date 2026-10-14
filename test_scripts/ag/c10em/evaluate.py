#!/usr/bin/env python3
"""
Compare models scored by `em.py score`, and match profiles of two model files.

  evaluate.py compare NEW/scores.tsv BASE/scores.tsv [--df-new M] [--df-base M] [--boot 10000] [--seed 1]
  evaluate.py match A.nex B.nex

compare: per alignment the information-criterion gain of NEW over BASE, with
the paper's second-order AIC  AIC = -2 lnL + 2m + 2m(m+1)/(N-m-1)  and
BIC = -2 lnL + m ln N, where N is the alignment length and m the number of
model-specific parameters (--df-new / --df-base, e.g. 9 for re-optimised C10
weights; with m = 0 both gains are 2 dlnL). Reports:
  - per-alignment dlnL: median, quartiles, number > 0, total;
  - number of alignments with dAIC > 10 either way;
  - mean and median dAIC per site and dBIC per site (unweighted over alignments),
    a 95% bootstrap CI of the mean dAIC/site (alignments resampled) and a
    one-sided Wilcoxon signed-rank p-value (H1: new better);
  - Kishino-Hasegawa wins/losses at p < 0.01 and p < 0.05 when both score
    directories kept site log-likelihoods (em.py score --wsl): two-sided RELL
    test on the summed site differences (10,000 centred bootstrap replicates).

match: exact minimum-cost assignment of the profiles by symmetric KL, with
mutual-nearest-neighbour flags and both weights.
"""
import os
import sys
import math
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_mdef, match_profiles


def read_scores(path):
    out = {}
    with open(path) as f:
        next(f)
        for l in f:
            i, nsite, lnl, _ = l.rstrip("\n").split("\t")
            if lnl:
                out[i] = (int(nsite), float(lnl))
    return out


def wilcoxon_greater(d):
    """one-sided p-value of the signed-rank test that the median of d is > 0"""
    d = np.asarray([x for x in d if x != 0.0])
    n = len(d)
    if n == 0:
        return 1.0
    a = np.abs(d)
    order = np.argsort(a)
    ranks = np.empty(n)
    i = 0
    while i < n:       # average ranks over ties
        j = i
        while j + 1 < n and a[order[j + 1]] == a[order[i]]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1
        i = j + 1
    w = ranks[d > 0].sum()
    mean = n * (n + 1) / 4.0
    _, counts = np.unique(a, return_counts=True)
    var = n * (n + 1) * (2 * n + 1) / 24.0 - np.sum(counts ** 3 - counts) / 48.0
    z = (w - mean - 0.5) / math.sqrt(var)
    return 0.5 * math.erfc(z / math.sqrt(2))


def aicc(lnl, m, n):
    return -2 * lnl + 2 * m + (2.0 * m * (m + 1) / (n - m - 1) if m else 0.0)


def bic(lnl, m, n):
    return -2 * lnl + m * math.log(n)


def site_lnl(path):
    for l in open(path):
        if l.startswith("Site_Lh"):
            return np.array([float(x) for x in l.split()[1:]])
    return None


def kh_rell(x, rng, boot=10000):
    """two-sided KH p-value by RELL: centred bootstrap of the summed site differences"""
    n = len(x)
    reps = np.empty(boot)
    for k in range(0, boot, 1000):     # chunks keep memory small for long alignments
        b = min(1000, boot - k)
        idx = rng.randint(n, size=(b, n))
        reps[k:k + b] = x[idx].sum(axis=1)
    obs = x.sum()
    return float(np.mean(np.abs(reps - obs) >= abs(obs)))


def kh_counts(dir_a, dir_b, ids, seed=1):
    """{alpha: (better, worse)} and the number of alignments tested"""
    rng = np.random.RandomState(seed)
    res = {0.01: [0, 0], 0.05: [0, 0]}
    n = 0
    for i in ids:
        pa, pb = os.path.join(dir_a, "sitelh", i + ".sitelh"), os.path.join(dir_b, "sitelh", i + ".sitelh")
        if not (os.path.exists(pa) and os.path.exists(pb)):
            continue
        x = site_lnl(pa) - site_lnl(pb)
        n += 1
        if not np.any(x):
            continue
        p = kh_rell(x, rng)
        for a in res:
            if p < a:
                res[a][0 if x.sum() > 0 else 1] += 1
    return res, n


def summarize(new, base, df_new=0, df_base=0, boot=10000, seed=1):
    """all comparison metrics of NEW over BASE as a dict (see the module docstring)"""
    a, b = read_scores(new), read_scores(base)
    ids = sorted(set(a) & set(b))
    N = np.array([a[i][0] for i in ids], float)
    dl = np.array([a[i][1] - b[i][1] for i in ids])
    daic = np.array([aicc(b[i][1], df_base, b[i][0]) - aicc(a[i][1], df_new, a[i][0]) for i in ids])
    dbic = np.array([bic(b[i][1], df_base, b[i][0]) - bic(a[i][1], df_new, a[i][0]) for i in ids])
    d = daic / N
    rng = np.random.RandomState(seed)
    boots = np.array([d[rng.randint(len(d), size=len(d))].mean() for _ in range(boot)])
    lo, hi = np.percentile(boots, [2.5, 97.5])
    q1, q3 = np.percentile(dl, [25, 75])
    kh, nkh = kh_counts(os.path.dirname(os.path.abspath(new)), os.path.dirname(os.path.abspath(base)), ids, seed)
    return dict(n=len(ids), missing_new=len(set(b) - set(a)), missing_base=len(set(a) - set(b)),
                dlnl_median=float(np.median(dl)), dlnl_q1=float(q1), dlnl_q3=float(q3), dlnl_pos=int((dl > 0).sum()),
                dlnl_total=float(dl.sum()), daic10_new=int((daic > 10).sum()), daic10_base=int((daic < -10).sum()),
                aic_site_mean=float(d.mean()), aic_site_median=float(np.median(d)), aic_site_ci=(float(lo), float(hi)),
                bic_site_mean=float((dbic / N).mean()), wilcoxon_p=wilcoxon_greater(d),
                kh_n=nkh, kh01=tuple(kh[0.01]), kh05=tuple(kh[0.05]))


def compare(args):
    s = summarize(args.new, args.base, args.df_new, args.df_base, args.boot, args.seed)
    print("alignments compared: %d (missing: new %d, base %d)" % (s["n"], s["missing_new"], s["missing_base"]))
    print("dlnL per alignment  : median %+.2f, quartiles %+.2f .. %+.2f, > 0 in %d, total %+.2f"
          % (s["dlnl_median"], s["dlnl_q1"], s["dlnl_q3"], s["dlnl_pos"], s["dlnl_total"]))
    print("dAIC > 10 new / base: %d / %d" % (s["daic10_new"], s["daic10_base"]))
    print("mean dAIC/site      : %+.5f   (95%% bootstrap CI %+.5f .. %+.5f); median %+.5f"
          % (s["aic_site_mean"], s["aic_site_ci"][0], s["aic_site_ci"][1], s["aic_site_median"]))
    print("mean dBIC/site      : %+.5f" % s["bic_site_mean"])
    print("Wilcoxon one-sided p: %.3g" % s["wilcoxon_p"])
    if s["kh_n"]:
        print("KH (RELL) better / worse: p<0.01 %d / %d, p<0.05 %d / %d, of %d"
              % (s["kh01"] + s["kh05"] + (s["kh_n"],)))


def match(args):
    _, pa, wa, _ = read_mdef(args.a)
    _, pb, wb, _ = read_mdef(args.b)
    pairs, mnn = match_profiles(pa, pb)
    print("class_a\tclass_b\tsym_kl\tmutual_nn\tw_a\tw_b")
    for (i, j, c), m in zip(pairs, mnn):
        print("%d\t%d\t%.4f\t%s\t%.4f\t%.4f" % (i + 1, j + 1, c, "yes" if m else "no", wa[i] / wa.sum(), wb[j] / wb.sum()))
    print("mutual nearest neighbours: %d/%d; total matched sym-KL %.4f" % (sum(mnn), len(mnn), sum(c for _, _, c in pairs)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compare")
    c.add_argument("new")
    c.add_argument("base")
    c.add_argument("--df-new", type=int, default=0)
    c.add_argument("--df-base", type=int, default=0)
    c.add_argument("--boot", type=int, default=10000)
    c.add_argument("--seed", type=int, default=1)
    m = sub.add_parser("match")
    m.add_argument("a")
    m.add_argument("b")
    args = ap.parse_args()
    {"compare": compare, "match": match}[args.cmd](args)


if __name__ == "__main__":
    main()
