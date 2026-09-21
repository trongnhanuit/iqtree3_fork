#!/usr/bin/env python3
"""
Independent oracle for the --analytical-gradients engine (NumPy only).

Rebuilds the phylogenetic likelihood from scratch (dense transition
matrices via eigendecomposition, Felsenstein pruning in ordinary state
space, no scaling tricks) from the model parameters and tree that IQ-TREE
wrote with --ag-dump-gradient, then:

  1. checks that its own log-likelihood equals IQ-TREE's (precondition:
     the two programs must agree on the model before gradients mean
     anything);
  2. compares IQ-TREE's analytic dlogL/dt for every branch against central
     finite differences of the oracle's log-likelihood;
  3. compares the natural-parameter gradients (#nat_* lines) against finite
     differences of the oracle: dR per exchangeability entry (Q renormalised
     to mean rate 1 as IQ-TREE does), dpi per raw frequency entry (pi
     renormalised to sum 1), dw per raw mixture weight (no renormalisation),
     drate/dprop per raw category rate/proportion, dpinv_direct (p_inv in
     the invariant-site term only);
  4. compares dlogL/dQ per component (#dlogl_dQ lines, every entry of Q
     treated as independent) against finite differences of an explicit-Q
     likelihood (general eigendecomposition, no reversibility assumed).
     E[N_ab] = Q_ab dlogL/dQ_ab and E[T_a] = dlogL/dQ_aa follow from this.

Supported: reversible DNA/protein models, optional mixture components with
their own exchangeabilities/frequencies and weights, discrete rate
categories (+G/+R) and +I, unfused. Small alignments only.

Usage:
  oracle.py <alignment.phy> <prefix.aggrad.tsv> [--tol 1e-6] [--selftest]

--selftest perturbs one branch gradient and one natural gradient before
comparing and must then report FAILED; a passing self-test means the
oracle cannot see anything.

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
        else:  # interleaved continuation: append to the first sequence still short
            for n in order:
                if len(seqs[n]) < nsite:
                    seqs[n] += l.replace(" ", "").upper()
                    break
    for n in order:
        assert len(seqs[n]) == nsite, "sequence %s has length %d, expected %d" % (n, len(seqs[n]), nsite)
    return order, seqs, nsite


def read_dump(path):
    meta = {"mix": [], "nat": {}}
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
                elif key == "dlogl_dQ":
                    meta["mix"][int(vals[0])]["dQ"] = np.array([float(x) for x in vals[1:]])
                elif key in ("nat_dR", "nat_dpi"):
                    meta["nat"].setdefault(key, {})[int(vals[0])] = np.array([float(x) for x in vals[1:]])
                elif key in ("nat_dw", "nat_drate", "nat_dprop"):
                    meta["nat"][key] = np.array([float(x) for x in vals])
                elif key in ("nat_dpinv_direct", "nat_dpinv", "nat_dalpha"):
                    meta["nat"][key] = float(vals[0])
                elif key in ("cat_rates", "cat_props"):
                    meta[key] = np.array([float(x) for x in vals])
                elif key == "tree":
                    meta["tree"] = rest.strip()
                elif key == "theta":
                    meta["theta"] = vals
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

    return parse_node()


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


def reroot(node):
    """Make `node` the root: reverse the parent links on the path to the old
    root, moving each branch length to the node it now hangs from."""
    path = []
    n = node
    while n is not None:
        path.append(n)
        n = n.parent
    lengths = [n.length for n in path]      # branch (path[i], path[i+1]) has length lengths[i]
    for i in range(len(path) - 1):
        child, parent = path[i], path[i + 1]
        parent.children.remove(child)
        child.children.append(parent)
        parent.parent = child
        parent.length = lengths[i]
    node.parent = None
    node.length = 0.0
    return node


# ---------------- model ----------------

def build_Q(exch, freqs, nst):
    """IQ-TREE convention: exchangeabilities in upper-triangular order,
    Q_ij = r_ij * pi_j with pi renormalised to sum 1, rows sum to zero,
    scaled so the mean rate is 1."""
    pi = freqs / freqs.sum()
    R = np.zeros((nst, nst))
    k = 0
    for i in range(nst):
        for j in range(i + 1, nst):
            R[i, j] = R[j, i] = exch[k]
            k += 1
    assert k == len(exch), "exchangeability count %d does not match nstates %d" % (len(exch), nst)
    Q = R * pi[None, :]
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    mean_rate = -(pi * np.diag(Q)).sum()
    return Q / mean_rate


def eig_reversible(Q, freqs):
    """Symmetric eigendecomposition: Q = U diag(lam) U^-1 with U = D^-1/2 W, U^-1 = W^T D^1/2."""
    d = np.sqrt(freqs / freqs.sum())
    S = (d[:, None] * Q) / d[None, :]
    S = 0.5 * (S + S.T)
    lam, W = np.linalg.eigh(S)
    U = W / d[:, None]
    Uinv = W.T * d[None, :]
    return lam, U, Uinv


def eig_general(Q):
    """General (possibly non-reversible) eigendecomposition for an explicit Q."""
    lam, V = np.linalg.eig(Q)
    return lam, V, np.linalg.inv(V)


def expm_eig(lam, U, Uinv, t):
    P = (U * np.exp(lam * t)[None, :]) @ Uinv
    return np.real(P) if np.iscomplexobj(P) else P


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


class Oracle:
    """Caches the tip vectors and the constant-site masks for one alignment."""

    def __init__(self, root, seqs, nsite, seqtype, nst):
        self.root = root
        self.nodes = all_nodes(root)
        self.post = list(reversed(self.nodes))
        self.nsite = nsite
        self.nst = nst
        self.tips = {n.name: np.stack([tip_vector(seqs[n.name][s], seqtype, nst) for s in range(nsite)])
                     for n in self.nodes if n.name in seqs}
        mask = np.ones((nsite, nst))
        for name in seqs:
            mask = mask * np.stack([tip_vector(seqs[name][s], seqtype, nst) for s in range(nsite)])
        self.const_mask = mask

    def loglik(self, meta, Q_override=None):
        """Log-likelihood of the model in `meta`; Q_override maps a component
        index to an explicit Q matrix (used for the dlogL/dQ check)."""
        nst = self.nst
        comps = meta["mix"]
        cat_rates = meta["cat_rates"]
        cat_props = meta["cat_props"]
        pinv = float(meta["pinvar"])
        site_like = np.zeros(self.nsite)
        for m, comp in enumerate(comps):
            if Q_override is not None and m in Q_override:
                lam, U, Uinv = eig_general(Q_override[m])
            else:
                Q = build_Q(comp["exch"], comp["freqs"], nst)
                lam, U, Uinv = eig_reversible(Q, comp["freqs"])
            w = comp["weight"]
            pi = comp["freqs"] / comp["freqs"].sum()
            for c in range(len(cat_rates)):
                r, p = cat_rates[c], cat_props[c]
                P = {id(n): expm_eig(lam, U, Uinv, r * n.length) for n in self.nodes if n is not self.root}
                partial = {}
                for n in self.post:
                    # a leaf keeps its tip vector even when it is the root (after reroot)
                    L = self.tips[n.name] if n.name in self.tips else np.ones((self.nsite, nst))
                    for ch in n.children:
                        L = L * (partial[id(ch)] @ P[id(ch)].T)
                    partial[id(n)] = L
                site_like += w * p * (partial[id(self.root)] @ pi)
        if pinv > 0.0:
            # constant-site term: p_inv * sum over states compatible with every
            # taxon of pi_bar_x, pi_bar = weighted average of the (normalised)
            # component frequencies with the raw mixture weights
            pi_bar = sum(comp["weight"] * comp["freqs"] / comp["freqs"].sum() for comp in comps)
            site_like += pinv * (self.const_mask @ pi_bar)
        return float(np.log(site_like).sum())


# ---------------- comparison helpers ----------------

class Report:
    def __init__(self, tol):
        self.tol = tol
        self.n_fail = 0
        self.n = 0
        self.max_rel = 0.0

    def compare(self, what, analytic, fd, gmax):
        err = abs(analytic - fd) / max(abs(analytic), abs(fd), 1e-6 * gmax)
        self.max_rel = max(self.max_rel, err)
        self.n += 1
        status = "PASS" if err <= self.tol else "FAIL"
        if status == "FAIL":
            self.n_fail += 1
        print("%-26s analytic=%.10g oracle_fd=%.10g rel=%.2e %s" % (what, analytic, fd, err, status))


def central_fd(fun, x0, h):
    """Richardson-extrapolated central difference (steps h and h/2): the
    truncation error is O(h^4), so h can be large enough that rounding in
    the log-likelihood (about 1e-16 |logL| / h) stays far below the tolerance."""
    d1 = (fun(x0 + h) - fun(x0 - h)) / (2 * h)
    d2 = (fun(x0 + h / 2) - fun(x0 - h / 2)) / h
    return (4 * d2 - d1) / 3


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
    nst = int(meta["nstates"])
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
    # dlogL/dQ (every entry independent, so not reversible) depends on where
    # the tree is rooted: root the oracle where IQ-TREE's engine does
    if "root_side" in meta:
        side = frozenset(meta["root_side"].split("|"))
        all_leaves = frozenset(order)
        target = None
        for n in all_nodes(root):
            if n is root:
                continue
            if n.leaves == side:
                target = n
            elif n.leaves == all_leaves - side:
                target = n.parent
        if target is None:
            print("oracle: root side %s not found in tree" % meta["root_side"])
            return 1
        root = reroot(target)
        leaf_sets(root)
    orc = Oracle(root, seqs, nsite, seqtype, nst)

    # 1. likelihood agreement
    ll = orc.loglik(meta)
    ll_iq = float(meta["logl"])
    rel = abs(ll - ll_iq) / max(1.0, abs(ll_iq))
    print("oracle: logl oracle=%.10f iqtree=%.10f rel_diff=%.2e" % (ll, ll_iq, rel))
    if rel > 1e-8:
        print("oracle: FAILED precondition: log-likelihoods disagree; model or tree mismatch")
        return 1
    rep = Report(tol)

    # 2. branch gradients by central differences on the oracle likelihood
    all_leaves = frozenset(order)
    by_side = {}
    for n in all_nodes(root):
        if n is root:
            continue
        by_side[n.leaves] = n
        by_side[all_leaves - n.leaves] = n
    gmax = max(abs(r[2]) for r in rows) if rows else 1.0
    for k, (idx, length, grad, label) in enumerate(rows):
        side = frozenset(label.split("|"))
        node = by_side.get(side)
        if node is None:
            print("oracle: branch %d (%s) not found in tree" % (idx, label))
            rep.n_fail += 1
            continue
        h = 1e-4 * max(1.0, length)
        saved = node.length

        def f(x, node=node):
            node.length = x
            return orc.loglik(meta)
        fd = central_fd(f, saved, h)
        node.length = saved
        analytic = grad * 1.01 + 1e-3 if (selftest and k == 0) else grad
        rep.compare("branch %d len=%.5f" % (idx, length), analytic, fd, gmax)
    n_branch = rep.n

    # 3. natural-parameter gradients
    nat = meta["nat"]
    if nat:
        comps = meta["mix"]

        def fd_on(arr, k, h, what, analytic, gmax):
            saved = arr[k]

            def f(x):
                arr[k] = x
                return orc.loglik(meta)
            fd = central_fd(f, saved, h)
            arr[k] = saved
            rep.compare(what, analytic, fd, gmax)

        for m, comp in enumerate(comps):
            if "nat_dR" in nat and m in nat["nat_dR"]:
                dR = nat["nat_dR"][m]
                if selftest and m == 0 and len(dR):
                    dR = dR.copy()
                    dR[0] = dR[0] * 1.01 + 1e-3
                g = max(np.abs(dR).max(), 1e-300)
                entries = range(len(dR)) if len(dR) <= 30 else sorted(set([0, 1, 5, 17, 42, 63, 99, 128, 150, 177, 188, 189]))
                for k in entries:
                    if k >= len(dR):
                        continue
                    fd_on(comp["exch"], k, 1e-4 * max(1.0, abs(comp["exch"][k])), "dR[%d][%d]" % (m, k), dR[k], g)
            if "nat_dpi" in nat and m in nat["nat_dpi"]:
                dpi = nat["nat_dpi"][m]
                g = max(np.abs(dpi).max(), 1e-300)
                for k in range(nst):
                    if comp["freqs"][k] <= 1e-10:
                        continue
                    fd_on(comp["freqs"], k, 1e-4 * comp["freqs"][k], "dpi[%d][%d]" % (m, k), dpi[k], g)
        if "nat_dw" in nat and len(comps) > 1:
            dw = nat["nat_dw"]
            g = max(np.abs(dw).max(), 1e-300)
            for m, comp in enumerate(comps):
                saved = comp["weight"]

                def f(x, comp=comp):
                    comp["weight"] = x
                    return orc.loglik(meta)
                fd = central_fd(f, saved, 1e-4 * saved)
                comp["weight"] = saved
                rep.compare("dw[%d]" % m, dw[m], fd, g)
        for key, arr_key in (("nat_drate", "cat_rates"), ("nat_dprop", "cat_props")):
            if key in nat and len(meta[arr_key]) > 1:
                vals = nat[key]
                g = max(np.abs(vals).max(), 1e-300)
                for c in range(len(vals)):
                    fd_on(meta[arr_key], c, 1e-4 * max(1.0, abs(meta[arr_key][c])), "%s[%d]" % (key[4:], c), vals[c], g)
        if "nat_dpinv_direct" in nat and float(meta["pinvar"]) > 0.0:
            saved = float(meta["pinvar"])

            def f(x):
                meta["pinvar"] = x
                return orc.loglik(meta)
            # the initial p_inv can be 1e-6; the invariant term is linear in p, so an absolute step is safe
            fd = central_fd(f, saved, 1e-4 * max(saved, 1e-2))
            meta["pinvar"] = saved
            rep.compare("dpinv_direct", nat["nat_dpinv_direct"], fd, max(abs(nat["nat_dpinv_direct"]), 1e-300))

    # 4. dlogL/dQ per component: every entry independent, explicit-Q likelihood
    for m, comp in enumerate(meta["mix"]):
        if "dQ" not in comp:
            continue
        D = comp["dQ"].reshape(nst, nst)
        Q0 = build_Q(comp["exch"], comp["freqs"], nst)
        g = max(np.abs(D).max(), 1e-300)
        if nst <= 4:
            entries = [(a, b) for a in range(nst) for b in range(nst)]
        else:
            rng = np.random.RandomState(12345)
            entries = [(a, a) for a in (0, 7, 19)] + [tuple(rng.randint(0, nst, 2)) for _ in range(9)]
        for a, b in entries:
            h = 1e-4 * max(1.0, abs(Q0[a, b]))

            def f(x):
                Q = Q0.copy()
                Q[a, b] = x
                return orc.loglik(meta, {m: Q})
            fd = central_fd(f, Q0[a, b], h)
            rep.compare("dlogl_dQ[%d][%d,%d]" % (m, a, b), D[a, b], fd, g)

    print("oracle: %d branches + %d parameter derivatives, max_rel=%.2e, n_fail=%d -> %s"
          % (n_branch, rep.n - n_branch, rep.max_rel, rep.n_fail, "PASSED" if rep.n_fail == 0 else "FAILED"))
    if selftest:
        # the perturbed gradients must be detected
        return 0 if rep.n_fail > 0 else 1
    return 0 if rep.n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
