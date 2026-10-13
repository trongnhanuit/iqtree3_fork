#!/usr/bin/env python3
"""
Pooled M-step of the C10 re-estimation EM (design doc section 16).

Paper semantics: every class shares one event rate u per unit branch length
(FMIX rate field rho_k proportional to phi(pi_k) = 1 - sum(pi_k^2)), so with
u held fixed the profile update is exact and closed form up to one 1-D root:

  pi_b = (n_b + eps) / (mu + h_b),   n_b = sum_aln (root_b + subin_b),
                                     h_b = u * sum_{a != b} dwell_a,
  w_k  = (M_k + eps_w) / (N + K eps_w),

with mu the unique root of sum_b pi_b = 1 on (-min h, sum(n+eps) - min h].
Trees are then multiplied by c = sum w' phi(pi') / sum w phi(pi), which keeps
u*t (event time) fixed, so the next round starts at a likelihood no lower
than this round's final one (the GEM inequality the driver checks).

CLI (debugging): mstep.py <model.nex> <stats.npz> <out.nex> [--mode full|weights]
"""
import sys
import numpy as np
from common import phi, write_mdef, read_mdef

EPS = 0.1      # Dir(1+EPS) prior on each profile
EPS_W = 0.1    # Dir(1+EPS_W) prior on the weights


def pool(dumps):
    """sum the per-alignment statistics (dumps in alignment-id order)"""
    K = len(dumps[0]["classes"])
    S = len(dumps[0]["classes"][0]["freq"])
    P = {"mass": np.zeros(K), "root": np.zeros((K, S)), "subin": np.zeros((K, S)),
         "dwell": np.zeros((K, S)), "nsite": 0.0, "logl": 0.0, "logl_init": 0.0}
    for d in dumps:
        for k, cl in enumerate(d["classes"]):
            P["mass"][k] += cl["mass"].sum()
            P["root"][k] += cl["root"]
            P["subin"][k] += cl["subin"]
            P["dwell"][k] += cl["dwell"]
        P["nsite"] += d["nsite"]
        P["logl"] += d["logl"]
        P["logl_init"] += d["logl_init"]
    return P


def log_prior(profiles, weights):
    w = np.asarray(weights) / np.sum(weights)
    return EPS * sum(float(np.sum(np.log(p))) for p in profiles) + EPS_W * float(np.sum(np.log(w)))


def solve_profile(n, h):
    """argmax sum (n_b+EPS) log pi_b - sum pi_b h_b on the simplex (strictly concave)"""
    a = n + EPS
    hmin = float(h.min())
    lo, hi = -hmin, float(a.sum()) - hmin
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if np.sum(a / (mid + h)) > 1.0:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 1e-15 * max(1.0, abs(hi)):
            break
    pi = a / (hi + h)
    return pi / pi.sum()


def mstep(P, profiles, weights, mode="full"):
    """-> (profiles', weights', rates', c); rates are paper semantics, c rescales the trees"""
    K = len(profiles)
    w = np.asarray(weights, float) / np.sum(weights)
    phis = np.array([phi(p) for p in profiles])
    u = 1.0 / float(np.dot(w, phis))
    new_w = (P["mass"] + EPS_W) / (P["nsite"] + K * EPS_W)
    if mode == "weights":
        new_p = [np.asarray(p, float) / np.sum(p) for p in profiles]
    else:
        new_p = []
        for k in range(K):
            T = P["dwell"][k]
            new_p.append(solve_profile(P["root"][k] + P["subin"][k], u * (T.sum() - T)))
    new_phis = np.array([phi(p) for p in new_p])
    c = float(np.dot(new_w, new_phis)) / float(np.dot(w, phis))
    return new_p, new_w, new_phis, c


def q_function(P, profiles, weights, u):
    """expected complete-data log-likelihood + log-prior at fixed u (should not drop across an M-step)"""
    w = np.asarray(weights) / np.sum(weights)
    q = float(np.dot(P["mass"], np.log(w)))
    for k, p in enumerate(profiles):
        T = P["dwell"][k]
        q += float(np.dot(P["root"][k] + P["subin"][k], np.log(p))) - u * float(np.dot(p, T.sum() - T))
    return q + log_prior(profiles, weights)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(2)
    name, profiles, weights, _ = read_mdef(sys.argv[1])
    z = np.load(sys.argv[2])
    P = {k: z[k] for k in z.files}
    mode = sys.argv[sys.argv.index("--mode") + 1] if "--mode" in sys.argv else "full"
    new_p, new_w, rates, c = mstep(P, profiles, weights, mode)
    write_mdef(sys.argv[3], name, new_p, new_w, rates)
    print("tree scale factor c = %.12g" % c)
