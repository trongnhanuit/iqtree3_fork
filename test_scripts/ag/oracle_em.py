#!/usr/bin/env python3
"""
Independent check of --ag-dump-em (<prefix>.agem.tsv), NumPy only.

Supports mixtures whose components have all exchangeabilities equal to 1
(F81 / POISSON with user frequencies, e.g. JC/POISSON+G4+FMIX{...}) and
discrete rate categories, without +I. Each class's Q is rebuilt from the
dumped frequencies and rate field (tns), then:

  1. the log-likelihood must equal IQ-TREE's (precondition);
  2. class mass per rate category and root-state counts are recomputed by
     plain pruning, rooted where IQ-TREE's engine was (#root_side);
  3. subin_b = sum_{a!=b} Q_ab dlogL/dQ_ab by a central difference that scales
     the off-diagonal entries of column b of Q, and dwell_a = dlogL/dQ_aa by a
     central difference on Q_aa (explicit-Q likelihood, entries independent);
  4. identities: sum of all masses = nsite, sum_b root = mass per class.

Usage:
  oracle_em.py <alignment.phy> <prefix.agem.tsv> [--tol 1e-6] [--selftest]

--selftest perturbs one subin entry first and must then report FAILED.
Exit code 0 = all within tolerance, 1 = failure, 2 = usage.
"""
import sys
import numpy as np
from oracle import (read_phylip, parse_newick, leaf_sets, reroot, all_nodes,
                    eig_general, expm_eig, tip_vector, central_fd, Report)


def read_agem(path):
    meta = {"classes": [], "end": False}
    with open(path) as f:
        for l in f:
            l = l.rstrip("\n")
            if not l.startswith("#"):
                continue
            key, _, rest = l[1:].partition("\t")
            vals = rest.split("\t") if rest else []
            if key == "END":
                meta["end"] = True
            elif key == "class":
                # "#class <m> weight <w> tns <t>"
                meta["classes"].append({"weight": float(vals[2]), "tns": float(vals[4])})
            elif key in ("freq", "mass", "root", "subin", "dwell"):
                meta["classes"][int(vals[0])][key] = np.array([float(x) for x in vals[1:]])
            elif key in ("cat_rates", "cat_props"):
                meta[key] = np.array([float(x) for x in vals])
            elif key == "tree":
                meta["tree"] = rest.strip()
            else:
                meta[key] = vals[0] if vals else ""
    return meta


def f81_Q(pi, tns):
    """All exchangeabilities 1: Q_ab = pi_b, rows sum to zero, mean rate tns."""
    S = len(pi)
    Q = np.tile(pi, (S, 1))
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    return Q * (tns / -(pi * np.diag(Q)).sum())


class EmOracle:
    def __init__(self, root, seqs, nsite, seqtype, nst):
        self.root = root
        self.nodes = all_nodes(root)
        self.post = list(reversed(self.nodes))
        self.nsite = nsite
        self.nst = nst
        self.tips = {n.name: np.stack([tip_vector(seqs[n.name][s], seqtype, nst) for s in range(nsite)])
                     for n in self.nodes if n.name in seqs}

    def class_terms(self, meta, Qs):
        """per (m, r): root partial (nsite x S) of class (m, r) times w_m p_r, without pi"""
        out = {}
        for m, cl in enumerate(meta["classes"]):
            lam, U, Uinv = eig_general(Qs[m])
            for r, (rate, prop) in enumerate(zip(meta["cat_rates"], meta["cat_props"])):
                P = {id(n): expm_eig(lam, U, Uinv, rate * n.length) for n in self.nodes if n is not self.root}
                partial = {}
                for n in self.post:
                    L = self.tips[n.name] if n.name in self.tips else np.ones((self.nsite, self.nst))
                    for ch in n.children:
                        L = L * (partial[id(ch)] @ P[id(ch)].T)
                    partial[id(n)] = L
                out[(m, r)] = cl["weight"] * prop * partial[id(self.root)]
        return out

    def loglik(self, meta, Qs):
        terms = self.class_terms(meta, Qs)
        site = sum(terms[(m, r)] @ meta["classes"][m]["freq"] for (m, r) in terms)
        return float(np.log(site).sum()), terms, site


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    aln_path, dump_path = argv[1], argv[2]
    tol = float(argv[argv.index("--tol") + 1]) if "--tol" in argv else 1e-6
    selftest = "--selftest" in argv

    order, seqs, nsite = read_phylip(aln_path)
    meta = read_agem(dump_path)
    if not meta["end"]:
        print("oracle_em: FAILED: no #END line (truncated dump)")
        return 1
    nst = int(meta["nstates"])
    seqtype = "DNA" if nst == 4 else "AA"
    classes = meta["classes"]
    if selftest:
        classes[0]["subin"] = classes[0]["subin"].copy()
        classes[0]["subin"][0] = classes[0]["subin"][0] * 1.01 + 1e-3

    # root the tree where IQ-TREE's engine did (dlogL/dQ and root counts depend on it)
    root = parse_newick(meta["tree"])
    leaf_sets(root)
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
    if target is None and root.leaves == side:
        target = root
    if target is None:
        print("oracle_em: root side %s not found in tree" % meta["root_side"])
        return 1
    root = reroot(target)
    leaf_sets(root)
    orc = EmOracle(root, seqs, nsite, seqtype, nst)

    # 1. likelihood agreement
    Qs = [f81_Q(cl["freq"], cl["tns"]) for cl in classes]
    ll, terms, site = orc.loglik(meta, Qs)
    ll_iq = float(meta["logl"])
    rel = abs(ll - ll_iq) / max(1.0, abs(ll_iq))
    print("oracle_em: logl oracle=%.10f iqtree=%.10f rel_diff=%.2e" % (ll, ll_iq, rel))
    if rel > 1e-8:
        print("oracle_em: FAILED precondition: log-likelihoods disagree; model or tree mismatch")
        return 1
    rep = Report(tol)

    # 2. class mass and root-state counts by pruning
    for m, cl in enumerate(classes):
        ref_mass = np.array([((terms[(m, r)] @ cl["freq"]) / site).sum() for r in range(len(meta["cat_rates"]))])
        ref_root = sum((terms[(m, r)] * cl["freq"][None, :]) / site[:, None]
                       for r in range(len(meta["cat_rates"]))).sum(axis=0)
        g = max(ref_mass.max(), 1e-300)
        for r in range(len(ref_mass)):
            rep.compare("mass[%d][%d]" % (m, r), cl["mass"][r], ref_mass[r], g)
        g = max(ref_root.max(), 1e-300)
        for b in range(nst):
            rep.compare("root[%d][%d]" % (m, b), cl["root"][b], ref_root[b], g)

    # 3. subin and dwell by central differences on the explicit-Q likelihood
    for m, cl in enumerate(classes):
        Q0 = Qs[m]
        gs = max(np.abs(cl["subin"]).max(), 1e-300)
        for b in range(nst):
            M = np.zeros_like(Q0)
            M[:, b] = Q0[:, b]
            M[b, b] = 0.0

            def f(eps, M=M):
                Q = list(Qs)
                Q[m] = Q0 + eps * M
                return orc.loglik(meta, Q)[0]
            rep.compare("subin[%d][%d]" % (m, b), cl["subin"][b], central_fd(f, 0.0, 1e-4), gs)
        gd = max(np.abs(cl["dwell"]).max(), 1e-300)
        for a in range(nst):
            def f(x, a=a):
                Q = list(Qs)
                Qa = Q0.copy()
                Qa[a, a] = x
                Q[m] = Qa
                return orc.loglik(meta, Q)[0]
            rep.compare("dwell[%d][%d]" % (m, a), cl["dwell"][a], central_fd(f, Q0[a, a], 1e-4 * max(1.0, abs(Q0[a, a]))), gd)

    # 4. identities
    total = sum(cl["mass"].sum() for cl in classes)
    id_mass = abs(total - float(meta["nsite"])) / float(meta["nsite"])
    id_root = max(abs(cl["root"].sum() - cl["mass"].sum()) / max(cl["mass"].sum(), 1e-300) for cl in classes)
    id_ok = id_mass <= 1e-10 and id_root <= 1e-8
    print("oracle_em: identities sum(mass)=nsite rel=%.2e, sum(root)=mass rel=%.2e -> %s"
          % (id_mass, id_root, "PASS" if id_ok else "FAIL"))
    if not id_ok:
        rep.n_fail += 1

    print("oracle_em: %d statistics, max_rel=%.2e, n_fail=%d -> %s"
          % (rep.n, rep.max_rel, rep.n_fail, "PASSED" if rep.n_fail == 0 else "FAILED"))
    if selftest:
        return 0 if rep.n_fail > 0 else 1
    return 0 if rep.n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
