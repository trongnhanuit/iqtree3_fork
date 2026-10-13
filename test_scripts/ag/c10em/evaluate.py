#!/usr/bin/env python3
"""
Compare models scored by `em.py score`, and match profiles of two model files.

  evaluate.py compare NEW/scores.tsv BASE/scores.tsv [--boot 10000] [--seed 1]
  evaluate.py match A.nex B.nex

compare: per alignment d_i = 2 (lnL_new - lnL_base) / nsite_i, i.e. the
AIC-per-site gain when both models have the same number of parameters (as
for fixed-weight C10 variants). Reports the unweighted mean and median of d,
wins/losses, a one-sided Wilcoxon signed-rank p-value (H1: new better;
normal approximation with tie and zero corrections) and a 95% bootstrap CI
of the mean (alignments resampled).

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


def compare(args):
    a, b = read_scores(args.new), read_scores(args.base)
    ids = sorted(set(a) & set(b))
    d = np.array([2 * (a[i][1] - b[i][1]) / a[i][0] for i in ids])
    rng = np.random.RandomState(args.seed)
    boots = np.array([d[rng.randint(len(d), size=len(d))].mean() for _ in range(args.boot)])
    lo, hi = np.percentile(boots, [2.5, 97.5])
    print("alignments compared: %d (missing: new %d, base %d)" % (len(ids), len(set(b) - set(a)), len(set(a) - set(b))))
    print("mean dAIC/site      : %+.5f   (95%% bootstrap CI %+.5f .. %+.5f)" % (d.mean(), lo, hi))
    print("median dAIC/site    : %+.5f" % np.median(d))
    print("new better / worse  : %d / %d" % ((d > 0).sum(), (d < 0).sum()))
    print("Wilcoxon one-sided p: %.3g" % wilcoxon_greater(d))
    print("total dlnL          : %+.3f" % sum(a[i][1] - b[i][1] for i in ids))


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
    c.add_argument("--boot", type=int, default=10000)
    c.add_argument("--seed", type=int, default=1)
    m = sub.add_parser("match")
    m.add_argument("a")
    m.add_argument("b")
    args = ap.parse_args()
    {"compare": compare, "match": match}[args.cmd](args)


if __name__ == "__main__":
    main()
