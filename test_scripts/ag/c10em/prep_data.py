#!/usr/bin/env python3
"""
Prepare protein alignments (and trees) for the C10 re-estimation pipeline.

  prep_data.py --aln-dir DIR --out DATA [--tree-dir DIR] [--bin IQTREE]
               [--tree-model WAG+G4] [--holdout 0.10] [--seed 1] [-j 24]

For every alignment file (FASTA or PHYLIP) in --aln-dir:
  - taxon names: characters outside [A-Za-z0-9_.-] become '_';
  - B, Z, J, U, O and other non-amino-acid letters become X ('.' and '~' become '-');
  - sequences without any amino acid and all-gap columns are dropped;
  - fewer than 4 remaining taxa: skipped (listed in DATA/skipped.tsv).
Trees: the file in --tree-dir with the same base name (pruned to the kept
taxa), else one IQ-TREE search per alignment under --tree-model with
-keep-ident (needs --bin). Trees are frozen: the pipeline never changes
their topology.

Writes DATA/aln/<id>.phy, DATA/tree/<id>.tree, DATA/index.tsv (id, ntaxa,
nsite), DATA/all.txt, DATA/train.txt and DATA/heldout.txt (a size-stratified
holdout: one alignment per block of 1/holdout alignments sorted by size),
DATA/prep.json and DATA/sha256.txt.
"""
import os
import re
import sys
import json
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import AA, read_alignment, write_phylip, parse_newick, write_newick, leaves, prune, sha256

AAX = set(AA) | {"X", "-"}


def clean_name(s):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", s)


def clean_seq(s):
    s = s.upper().replace(".", "-").replace("~", "-").replace("?", "-")
    return "".join(c if c in AAX else "X" for c in s)


def prep_one(path, tree_dir):
    """-> (id, names, seqs, newick or None, reason or None)"""
    base = os.path.basename(path)
    aid = clean_name(base.split(".")[0])
    try:
        names0, seqs0 = read_alignment(path)
    except Exception as e:
        return aid, None, None, None, "unreadable: %s" % e
    names, seqs, seen = [], {}, set()
    for n in names0:
        c = clean_name(n)
        if c in seen:
            return aid, None, None, None, "duplicate taxon name %s" % c
        seen.add(c)
        s = clean_seq(seqs0[n])
        if any(ch in AA for ch in s):
            names.append(c)
            seqs[c] = s
    if len(names) < 4:
        return aid, None, None, None, "fewer than 4 taxa with data"
    cols = [j for j in range(len(seqs[names[0]])) if any(seqs[n][j] != "-" for n in names)]
    seqs = {n: "".join(seqs[n][j] for j in cols) for n in names}
    newick = None
    if tree_dir:
        cand = [f for f in os.listdir(tree_dir) if f.split(".")[0] == base.split(".")[0]]
        text = open(os.path.join(tree_dir, cand[0])).read().strip() if cand else ""
        if text:    # an empty tree file counts as missing (a tree is then built)
            t = parse_newick(text)
            for leaf in _leaf_nodes(t):
                leaf.name = clean_name(leaf.name)
            t = prune(t, set(names))
            if t is None or set(leaves(t)) != set(names):
                return aid, None, None, None, "tree taxa do not match the alignment"
            newick = write_newick(t)
    return aid, names, seqs, newick, None


def _leaf_nodes(n):
    return [n] if not n.children else [x for c in n.children for x in _leaf_nodes(c)]


def build_tree(args, aid, out):
    """one IQ-TREE search under --tree-model (frozen afterwards)"""
    d = os.path.join(out, "treebuild", aid)
    os.makedirs(d, exist_ok=True)
    pre = os.path.join(d, "t")
    subprocess.run([args.bin, "-s", os.path.join(out, "aln", aid + ".phy"), "-m", args.tree_model, "-keep-ident",
                    "-nt", "1", "-seed", str(args.seed), "-quiet", "-redo", "-pre", pre],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=dict(os.environ, OMP_NUM_THREADS="1"))
    if not os.path.exists(pre + ".treefile"):
        return aid, None
    return aid, open(pre + ".treefile").read().strip()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tree-dir")
    ap.add_argument("--bin")
    ap.add_argument("--tree-model", default="WAG+G4")
    ap.add_argument("--holdout", type=float, default=0.10)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("-j", type=int, default=24)
    args = ap.parse_args()
    out = os.path.abspath(args.out)
    for sub in ("aln", "tree"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    files = sorted(os.path.join(args.aln_dir, f) for f in os.listdir(args.aln_dir) if not f.startswith("."))
    with ThreadPoolExecutor(max_workers=args.j) as pool:
        res = list(pool.map(lambda p: prep_one(p, args.tree_dir), files))
    kept, skipped, need_tree = {}, [], []
    for aid, names, seqs, newick, reason in res:
        if reason:
            skipped.append((aid, reason))
            continue
        if aid in kept:
            skipped.append((aid, "duplicate alignment id"))
            continue
        write_phylip(os.path.join(out, "aln", aid + ".phy"), names, seqs)
        kept[aid] = (len(names), len(seqs[names[0]]))
        if newick:
            with open(os.path.join(out, "tree", aid + ".tree"), "w") as f:
                f.write(newick + "\n")
        else:
            need_tree.append(aid)
    if need_tree:
        if not args.bin:
            sys.exit("%d alignments have no tree; give --tree-dir or --bin to build %s trees" % (len(need_tree), args.tree_model))
        print("building %d %s trees ..." % (len(need_tree), args.tree_model), flush=True)
        with ThreadPoolExecutor(max_workers=args.j) as pool:
            for aid, nwk in pool.map(lambda a: build_tree(args, a, out), need_tree):
                if nwk is None:
                    skipped.append((aid, "tree search failed"))
                    del kept[aid]
                    os.remove(os.path.join(out, "aln", aid + ".phy"))
                    continue
                with open(os.path.join(out, "tree", aid + ".tree"), "w") as f:
                    f.write(nwk + "\n")
    ids = sorted(kept)
    with open(os.path.join(out, "index.tsv"), "w") as f:
        f.write("id\tntaxa\tnsite\n")
        for i in ids:
            f.write("%s\t%d\t%d\n" % (i, kept[i][0], kept[i][1]))
    held = set()
    if args.holdout > 0:
        rng = np.random.RandomState(args.seed)
        by_size = sorted(ids, key=lambda i: (kept[i][0] * kept[i][1], i))
        block = max(1, int(round(1.0 / args.holdout)))
        for b in range(0, len(by_size), block):
            chunk = by_size[b:b + block]
            if len(chunk) == block:
                held.add(chunk[rng.randint(len(chunk))])
    for name, lst in (("all", ids), ("train", [i for i in ids if i not in held]), ("heldout", sorted(held))):
        with open(os.path.join(out, name + ".txt"), "w") as f:
            f.write("".join(i + "\n" for i in lst))
    with open(os.path.join(out, "skipped.tsv"), "w") as f:
        f.write("".join("%s\t%s\n" % s for s in skipped))
    with open(os.path.join(out, "prep.json"), "w") as f:
        json.dump(dict(aln_dir=os.path.abspath(args.aln_dir), tree_dir=args.tree_dir, tree_model=args.tree_model,
                       built_trees=len(need_tree), holdout=args.holdout, seed=args.seed, n=len(ids),
                       n_train=len(ids) - len(held), n_heldout=len(held), n_skipped=len(skipped)), f, indent=1)
    with open(os.path.join(out, "sha256.txt"), "w") as f:
        for i in ids:
            for sub, ext in (("aln", ".phy"), ("tree", ".tree")):
                p = os.path.join(out, sub, i + ext)
                f.write("%s  %s\n" % (sha256(p), os.path.relpath(p, out)))
    print("prepared %d alignments (%d train, %d held out), skipped %d -> %s"
          % (len(ids), len(ids) - len(held), len(held), len(skipped), out))


if __name__ == "__main__":
    main()
