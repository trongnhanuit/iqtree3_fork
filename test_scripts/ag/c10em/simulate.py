#!/usr/bin/env python3
"""
Simulate protein alignments from a fixed mixture with AliSim (pipeline tests, gate M4).

  simulate.py --bin IQTREE --mdef FILE -m NAME --out DIR [--n 50] [--taxa 40]
              [--sites 250] [--alpha 0.8] [--seed 1] [-j 24]
              [--trees-from INPUTS.tsv --data DATA [--skip N]]

Default: each alignment gets its own random Yule-Harding tree, --sites and
--alpha. With --trees-from (a round's inputs.tsv: id, alpha, tree), alignment k
reuses the k-th listed tree (after --skip), its alpha and the site count of
that alignment in DATA/index.tsv, i.e. realistic shapes and branch lengths.
Writes DIR/aln/sim_<k>.phy and DIR/tree/sim_<k>.tree, ready for
prep_data.py --aln-dir DIR/aln --tree-dir DIR/tree.
"""
import os
import sys
import shutil
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", required=True)
    ap.add_argument("--mdef", required=True)
    ap.add_argument("-m", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument("--taxa", type=int, default=40)
    ap.add_argument("--sites", type=int, default=250)
    ap.add_argument("--alpha", type=float, default=0.8)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--trees-from")
    ap.add_argument("--data")
    ap.add_argument("--skip", type=int, default=0)
    ap.add_argument("-j", type=int, default=24)
    args = ap.parse_args()
    out = os.path.abspath(args.out)
    for sub in ("aln", "tree", "work"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    specs = None
    if args.trees_from:
        nsite = {}
        with open(os.path.join(args.data, "index.tsv")) as f:
            next(f)
            for l in f:
                v = l.split("\t")
                nsite[v[0]] = int(v[2])
        rows = [l.rstrip("\n").split("\t") for l in open(args.trees_from)][args.skip:args.skip + args.n]
        specs = [(r[2], float(r[1]), nsite[r[0]]) for r in rows]
        if len(specs) < args.n:
            sys.exit("only %d trees available in %s" % (len(specs), args.trees_from))

    def one(k):
        pre = os.path.join(out, "work", "sim_%d" % k)
        if specs:
            with open(pre + ".in.tree", "w") as f:
                f.write(specs[k][0] + "\n")
            tree, alpha, sites = pre + ".in.tree", specs[k][1], specs[k][2]
        else:
            tree, alpha, sites = "RANDOM{yh/%d}" % args.taxa, args.alpha, args.sites
        cmd = [os.path.abspath(args.bin), "--alisim", pre, "-t", tree,
               "-mdef", os.path.abspath(args.mdef), "-m", args.m,
               "-a", repr(alpha), "--seqtype", "AA", "--length", str(sites),
               "-seed", str(args.seed + k), "-redo", "-quiet"]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=dict(os.environ, OMP_NUM_THREADS="1"))
        if not os.path.exists(pre + ".phy"):
            return k, False
        shutil.move(pre + ".phy", os.path.join(out, "aln", "sim_%d.phy" % k))
        if specs:
            shutil.move(tree, os.path.join(out, "tree", "sim_%d.tree" % k))
        elif os.path.exists(pre + ".treefile"):
            shutil.move(pre + ".treefile", os.path.join(out, "tree", "sim_%d.tree" % k))
        else:
            return k, False
        return k, True

    with ThreadPoolExecutor(max_workers=args.j) as pool:
        res = list(pool.map(one, range(args.n)))
    shutil.rmtree(os.path.join(out, "work"), ignore_errors=True)
    bad = [k for k, ok in res if not ok]
    print("simulated %d alignments (%d failed) -> %s" % (args.n - len(bad), len(bad), out))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
