#!/usr/bin/env python3
"""
Independent oracle for the --analytical-gradients engine (NumPy only).

Rebuilds the phylogenetic likelihood from scratch (dense transition
matrices via eigendecomposition, Felsenstein pruning in ordinary state
space, no scaling tricks) from the model parameters and tree that IQ-TREE
wrote with --ag-dump-gradient, then:

  1. checks that its own log-likelihood equals IQ-TREE's (precondition:
     the two programs must agree on the model before gradients mean
     anything), and
  2. compares IQ-TREE's analytic dlogL/dt for every branch against central
     finite differences of the oracle's log-likelihood.

Supported here (Stage 1): reversible DNA/protein models, optional mixture
components with their own exchangeabilities/frequencies and weights,
discrete rate categories (+G/+R) and +I, unfused. Small alignments only.

Usage:
  oracle.py <alignment.phy> <prefix.aggrad.tsv> [--tol 1e-6] [--selftest]

--selftest perturbs one analytic value before comparing and must then
report FAILED; a passing self-test means the oracle cannot see anything.

Exit code 0 = all comparisons within tolerance, 1 = failure, 2 = usage.
"""
import sys
import math
import numpy as np

DNA = "ACGT"
AA = "ARNDCQEGHILKMFPSTWYV"
DNA_AMBIG = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "M": "AC", "K": "GT", "S": "CG", "W": "AT",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
    "N": "ACGT", "X": "ACGT", "-": "ACGT", "?": "ACGT", ".": "ACGT", "~": "ACGT", "O": "ACGT",
}
AA_AMBIG = {"B": "ND", "Z": "QE", "J": "IL", "X": AA, "-": AA, "?": AA, ".": AA, "~": AA, "*": AA, "U": AA, "O": AA}


def read_phylip(path):
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f if l.strip()]
    ntax, nsite = [int(x) for x in lines[0].split()[:2]]
    seqs = {}
    order = []
    for l in lines[1:]:
        parts = l.split(None, 1)
        if len(parts) == 2 and parts[0] not in seqs and len(seqs) < ntax:
            seqs[parts[0]] = parts[1].replace(" ", "").upper()
            order.append(parts[0])
        else:  # interleaved continuation
            name = order[len([n for n in order if len(seqs[n]) >= nsite]) % ntax] if False else None
            # simple interleaved handling: append to sequences in order of insufficient length
            for n in order:
                if len(seqs[n]) < nsite:
                    seqs[n] += l.replace(" ", "").upper()
                    break
    for n in order:
        assert len(seqs[n]) == nsite, "sequence %s has length %d, expected %d" % (n, len(seqs[n]), nsite)
    return order, seqs, nsite


def read_dump(path):
    meta = {"mix": []}
    rows = []
    with open(path) as f:
        for l in f:
            l = l.rstrip("\n")
            if not l:
                continue
            if l.startswith("#"):
                key, _, rest = l[1:].partition("\t")
                vals = rest.split("\t") if rest else []
                if key == "mix":
                    # "#mix <m> weight <w> reversible <0|1>"
                    meta["mix"].append({"weight": float(vals[2]), "reversible": int(vals[4])})
                elif key == "exchangeabilities":
                    meta["mix"][int(vals[0])]["exch"] = np.array([float(x) for x in vals[1:]])
                elif key == "freqs":
                    meta["mix"][int(vals[0])]["freqs"] = np.array([float(x) for x in vals[1:]])
                elif key in ("cat_rates", "cat_props"):
                    meta[key] = np.array([float(x) for x in vals])
                elif key == "tree":
                    meta["tree"] = rest.strip()
                else:
                    meta[key] = vals[0] if vals else ""
            elif l.startswith("id\t"):
                continue
            else:
                idx, length, grad, label = l.split("\t")
                rows.append((int(idx), float(length), float(grad), label))
    return meta, rows


# ---------------- tree ----------------

class Node:
    __slots__ = ("name", "children", "length", "parent", "leaves")

    def __init__(self, name=None):
        self.name = name
        self.children = []
        self.length = 0.0
        self.parent = None
        self.leaves = None


def parse_newick(s):
    s = s.strip()
    if s.endswith(";"):
        s = s[:-1]
    pos = 0

    def parse_node():
        nonlocal pos
        node = Node()
        if s[pos] == "(":
            pos += 1
            while True:
                child = parse_node()
                child.parent = node
                node.children.append(child)
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
        # name
        start = pos
        while pos < len(s) and s[pos] not in ":,)(;":
            pos += 1
        node.name = s[start:pos].strip() or None
        if pos < len(s) and s[pos] == ":":
            pos += 1
            start = pos
            while pos < len(s) and s[pos] not in ",)(;":
                pos += 1
            node.length = float(s[start:pos])
        return node

    root = parse_node()
    return root


def all_nodes(root):
    out = []
    stack = [root]
    while stack:
        n = stack.pop()
        out.append(n)
        stack.extend(n.children)
    return out


def leaf_sets(root):
    for n in reversed(all_nodes(root)):
        if not n.children:
            n.leaves = frozenset([n.name])
        else:
            n.leaves = frozenset().union(*[c.leaves for c in n.children])


# ---------------- model ----------------

def build_Q(exch, freqs, nst):
    """IQ-TREE convention: exchangeabilities in upper-triangular order,
    Q_ij = r_ij * pi_j, rows sum to zero, scaled so the mean rate is 1."""
    R = np.zeros((nst, nst))
    k = 0
    for i in range(nst):
        for j in range(i + 1, nst):
            R[i, j] = R[j, i] = exch[k]
            k += 1
    assert k == len(exch), "exchangeability count %d does not match nstates %d" % (len(exch), nst)
    Q = R * freqs[None, :]
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    mean_rate = -(freqs * np.diag(Q)).sum()
    return Q / mean_rate


def eig_reversible(Q, freqs):
    """Symmetric eigendecomposition: Q = U diag(lam) U^-1 with U = D^-1/2 W, U^-1 = W^T D^1/2."""
    d = np.sqrt(freqs)
    S = (d[:, None] * Q) / d[None, :]
    S = 0.5 * (S + S.T)
    lam, W = np.linalg.eigh(S)
    U = W / d[:, None]
    Uinv = W.T * d[None, :]
    return lam, U, Uinv


def expm_rev(lam, U, Uinv, t):
    return (U * np.exp(lam * t)[None, :]) @ Uinv


# ---------------- likelihood ----------------

def tip_vector(ch, seqtype, nst):
    if seqtype == "DNA":
        states = DNA_AMBIG.get(ch, None)
        alphabet = DNA
    else:
        states = AA_AMBIG.get(ch, ch if ch in AA else AA)
        alphabet = AA
    if states is None:
        raise ValueError("unknown DNA character %r" % ch)
    v = np.zeros(nst)
    for s in states:
        v[alphabet.index(s)] = 1.0
    return v


def loglik(root, seqs, nsite, seqtype, meta):
    nst = int(meta["nstates"])
    comps = meta["mix"]
    cat_rates = meta["cat_rates"]
    cat_props = meta["cat_props"]
    pinv = float(meta["pinvar"])
    nodes = all_nodes(root)
    post = list(reversed(nodes))  # children before parents
    total = 0.0
    # precompute eigen systems
    eigs = []
    for comp in comps:
        Q = build_Q(comp["exch"], comp["freqs"], nst)
        eigs.append((comp, eig_reversible(Q, comp["freqs"])))
    # per site
    site_like = np.zeros(nsite)
    for comp, (lam, U, Uinv) in eigs:
        w = comp["weight"]
        pi = comp["freqs"]
        for c in range(len(cat_rates)):
            r, p = cat_rates[c], cat_props[c]
            P = {}
            for n in nodes:
                if n is not root:
                    P[id(n)] = expm_rev(lam, U, Uinv, r * n.length)
            partial = {}
            for n in post:
                if not n.children:
                    L = np.stack([tip_vector(seqs[n.name][s], seqtype, nst) for s in range(nsite)])  # nsite x nst
                else:
                    L = np.ones((nsite, nst))
                    for ch in n.children:
                        L = L * (partial[id(ch)] @ P[id(ch)].T)
                partial[id(n)] = L
            site_like += w * p * (partial[id(root)] @ pi)
    if pinv > 0.0:
        # constant-site term: p_inv * sum over states compatible with every taxon of pi_x
        # (IQ-TREE uses the mixture-averaged frequencies)
        pi_bar = sum(comp["weight"] * comp["freqs"] for comp in comps)
        for s in range(nsite):
            mask = np.ones(nst)
            for name in seqs:
                mask = mask * tip_vector(seqs[name][s], seqtype, nst)
            site_like[s] += pinv * (mask * pi_bar).sum()
    return float(np.log(site_like).sum())


# ---------------- main ----------------

def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    aln_path, dump_path = argv[1], argv[2]
    tol = 1e-6
    selftest = False
    for i, a in enumerate(argv):
        if a == "--tol":
            tol = float(argv[i + 1])
        if a == "--selftest":
            selftest = True

    order, seqs, nsite = read_phylip(aln_path)
    meta, rows = read_dump(dump_path)
    seqtype = meta.get("seqtype", "DNA")
    if int(meta.get("fused", "0")):
        print("oracle: fused mixture-rate models are not supported")
        return 2
    root = parse_newick(meta["tree"])
    # IQ-TREE prints an unrooted tree; the root of the parse is a trifurcation. Fine
    # for a reversible model: the likelihood does not depend on the root position.
    leaf_sets(root)
    missing = [n.name for n in all_nodes(root) if not n.children and n.name not in seqs]
    if missing:
        print("oracle: taxa in tree but not alignment:", missing)
        return 1

    # 1. likelihood agreement
    ll = loglik(root, seqs, nsite, seqtype, meta)
    ll_iq = float(meta["logl"])
    rel = abs(ll - ll_iq) / max(1.0, abs(ll_iq))
    print("oracle: logl oracle=%.10f iqtree=%.10f rel_diff=%.2e" % (ll, ll_iq, rel))
    if rel > 1e-8:
        print("oracle: FAILED precondition: log-likelihoods disagree; model or tree mismatch")
        return 1

    # 2. branch gradients by central differences on the oracle likelihood
    all_leaves = frozenset(order)
    by_side = {}
    for n in all_nodes(root):
        if n is root:
            continue
        by_side[n.leaves] = n
        by_side[all_leaves - n.leaves] = n
    n_fail = 0
    max_rel = 0.0
    gmax = max(abs(r[2]) for r in rows) if rows else 1.0
    for k, (idx, length, grad, label) in enumerate(rows):
        side = frozenset(label.split("|"))
        node = by_side.get(side)
        if node is None:
            print("oracle: branch %d (%s) not found in tree" % (idx, label))
            n_fail += 1
            continue
        h = 1e-5 * max(1.0, length)
        saved = node.length
        node.length = saved + h
        lp = loglik(root, seqs, nsite, seqtype, meta)
        node.length = max(saved - h, 1e-12)
        lm = loglik(root, seqs, nsite, seqtype, meta)
        node.length = saved
        fd = (lp - lm) / (2 * h)
        analytic = grad
        if selftest and k == 0:
            analytic = grad * 1.01 + 1e-3
        err = abs(analytic - fd) / max(abs(analytic), abs(fd), 1e-6 * gmax)
        max_rel = max(max_rel, err)
        status = "PASS" if err <= tol else "FAIL"
        if status == "FAIL":
            n_fail += 1
        print("branch %3d len=%.6f analytic=%.10g oracle_fd=%.10g rel=%.2e %s" % (idx, length, analytic, fd, err, status))
    print("oracle: %d branches, max_rel=%.2e, n_fail=%d -> %s" % (len(rows), max_rel, n_fail, "PASSED" if n_fail == 0 else "FAILED"))
    if selftest:
        # a perturbed gradient must be detected
        return 0 if n_fail > 0 else 1
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
