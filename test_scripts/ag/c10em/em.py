#!/usr/bin/env python3
"""
Driver of the C10 re-estimation EM (design doc section 16).

  em.py train  RUNDIR --data DATA --bin IQTREE [--option warm|cold|weights] [-j 24]
               [--chains N] [--halving 30:6,60:3|none] [--max-rounds 300] [--seed 1]
               [--stop-nats X] [--dry-run]
  em.py status RUNDIR
  em.py score  OUT --data DATA --list LIST --bin IQTREE [-j 24] (--mdef FILE -m NAME | -m MODEL)
               [--search] [--opt-weights] [--wsl]

train runs one arm in RUNDIR:
  warm     3 chains from the published C10 (paper rates) and 2 jittered copies
  cold     12 chains without C-series information (6 k-means++, 6 Dirichlet),
           halved to 6 at round 30 and to 3 at round 60 by projected final value
  weights  1 chain from the published C10 re-estimating only the weights (baseline C10pf-w)
Each round runs one IQ-TREE job per training alignment (fixed tree, branch
lengths and alpha optimised, --ag-dump-em), pools the dumps in alignment-id
order and applies mstep.py. All chains run in lockstep on one shared pool of
single-threaded jobs. Rerun the same command to resume. Results:
RUNDIR/best.nex, RUNDIR/finalists/*.nex, RUNDIR/summary.tsv.

score fits a model to every alignment of LIST and writes OUT/scores.tsv:
on its DATA tree (branch lengths and alpha), or with --search on a tree
IQ-TREE searches for that model; --opt-weights also re-optimises the
mixture weights (otherwise they stay as given).
"""
import os
import sys
import json
import time
import fcntl
import shutil
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (read_alignment, site_counts, published_c10, paper_rates, phi, write_mdef,
                    read_mdef, read_agem, echo_matches, scale_tree, match_profiles, sha256, read_list)
import mstep as M

NAME = "EMC"              # model name inside the round files
JOB_TIMEOUT = 1800        # seconds per IQ-TREE job
MONO_TOL = 0.01           # nats: GEM inequality tolerance (numerical only)
STOP_REMAIN = 4e-4        # nats per site: projected remaining gain to stop (100 nats on HSSP)
STOP_LAMBDA = 0.995       # raw median gain ratio must be below this
STOP_ROUNDS = 3           # consecutive rounds meeting both
HALVE_MARGIN = 50.0       # nats: chains this close to the cut-off survive halving
SAME_BASIN = 0.05         # max matched symmetric KL for two finalists to count as one basin


def log(msg):
    print(time.strftime("%H:%M:%S ") + msg, flush=True)


# ---------------- data ----------------

class Data:
    def __init__(self, path, list_name="train.txt"):
        self.path = os.path.abspath(path)
        self.size = {}
        with open(os.path.join(self.path, "index.tsv")) as f:
            next(f)
            for l in f:
                v = l.split("\t")
                self.size[v[0]] = int(v[1]) * int(v[2])
        excl = os.path.join(self.path, "exclude.txt")
        self.excluded = set(read_list(excl)) if os.path.exists(excl) else set()
        # a list name is looked up in DATA first, else taken as a path
        in_data = os.path.join(self.path, list_name)
        self.list_path = in_data if os.path.exists(in_data) else os.path.abspath(list_name)
        self.ids = sorted(i for i in read_list(self.list_path) if i not in self.excluded)

    def aln(self, i):
        return os.path.join(self.path, "aln", i + ".phy")

    def tree(self, i):
        with open(os.path.join(self.path, "tree", i + ".tree")) as f:
            return f.read().strip()


# ---------------- one IQ-TREE job ----------------

def run_job(job):
    """job: dict(bin, aln, tree, alpha, mdef, model, dir, check) -> (dump or None, error text)"""
    d = job["dir"]
    os.makedirs(d, exist_ok=True)
    pre = os.path.join(d, "j")
    old = read_agem(pre + ".agem.tsv")
    if old is not None and job["check"](old):
        return old, ""                        # resume: finished and matches this round's model
    if os.path.exists(pre + ".agem.tsv"):
        os.remove(pre + ".agem.tsv")
    with open(os.path.join(d, "in.tree"), "w") as f:
        f.write(job["tree"] + "\n")
    cmd = [job["bin"], "-s", job["aln"], "-te", os.path.join(d, "in.tree"), "-mdef", job["mdef"],
           "-m", job["model"], "-me", "1e-4", "-seed", "1", "-keep-ident", "-nt", "1",
           "--no-outfiles", "-quiet", "-redo", "-pre", pre, "--ag-dump-em"]
    if job.get("alpha"):
        cmd += ["-a", repr(job["alpha"]), "-optfromgiven"]
    env = dict(os.environ, OMP_NUM_THREADS="1")
    err = ""
    for attempt in (1, 2):
        try:
            p = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, env=env,
                               timeout=JOB_TIMEOUT, universal_newlines=True)
            rc, err = p.returncode, p.stderr[-2000:]
        except subprocess.TimeoutExpired:
            rc, err = -1, "timeout"
        dump = read_agem(pre + ".agem.tsv")
        if rc == 0 and dump is not None and job["check"](dump):
            return dump, ""
        if rc < 0 and err != "timeout":
            break                              # killed: do not retry
    return None, "exit %s: %s" % (rc, err.strip().splitlines()[-1] if err.strip() else "no dump")


def run_jobs(jobs, nproc):
    """run all jobs (longest first) on one pool; results in the input order"""
    order = sorted(range(len(jobs)), key=lambda i: -jobs[i]["size"])
    results = [None] * len(jobs)
    pool = ThreadPoolExecutor(max_workers=nproc)
    try:
        for i, r in zip(order, pool.map(run_job, [jobs[i] for i in order])):
            results[i] = r
    except KeyboardInterrupt:
        pool.shutdown(wait=False, cancel_futures=True)
        raise
    pool.shutdown()
    return results


# ---------------- chains ----------------

class Chain:
    def __init__(self, rundir, name):
        self.name = name
        self.dir = os.path.join(rundir, name)
        with open(os.path.join(self.dir, "chain.json")) as f:
            self.cfg = json.load(f)

    def rdir(self, r):
        return os.path.join(self.dir, "round_%03d" % r)

    def current_round(self):
        r = 0
        while os.path.exists(os.path.join(self.rdir(r + 1), "model.nex")):
            r += 1
        return r

    def status(self):
        p = os.path.join(self.dir, "STATUS")
        return open(p).read().strip() if os.path.exists(p) else "running"

    def set_status(self, s):
        with open(os.path.join(self.dir, "STATUS"), "w") as f:
            f.write(s + "\n")

    def summaries(self):
        out = []
        r = 0
        while os.path.exists(os.path.join(self.rdir(r), "summary.json")):
            with open(os.path.join(self.rdir(r), "summary.json")) as f:
                out.append(json.load(f))
            r += 1
        return out

    def model(self, r):
        return read_mdef(os.path.join(self.rdir(r), "model.nex"), NAME)

    def inputs(self, r):
        """id -> (alpha or None, newick)"""
        out = {}
        with open(os.path.join(self.rdir(r), "inputs.tsv")) as f:
            for l in f:
                i, a, t = l.rstrip("\n").split("\t")
                out[i] = (float(a) if a else None, t)
        return out


def new_chain(rundir, name, cfg, profiles, weights, data):
    c = os.path.join(rundir, name)
    r0 = os.path.join(c, "round_000")
    os.makedirs(r0, exist_ok=True)
    with open(os.path.join(c, "chain.json"), "w") as f:
        json.dump(cfg, f, indent=1)
    with open(os.path.join(r0, "inputs.tsv"), "w") as f:
        for i in data.ids:
            f.write("%s\t\t%s\n" % (i, data.tree(i)))
    write_mdef(os.path.join(r0, "model.nex"), NAME, profiles, weights, paper_rates(profiles))


def projected(gains):
    """(raw median ratio of the last 5 gains, projected remaining gain) or (None, inf)"""
    if len(gains) < 6:
        return None, float("inf")
    g = gains[-6:]
    ratios = [g[i + 1] / g[i] if g[i] != 0 else 0.0 for i in range(5)]
    lam = float(np.median(ratios))
    lc = min(max(lam, 0.0), 0.999999)
    return lam, gains[-1] * lc / (1.0 - lc)


def converged(summ, stop_nats=None):
    limit = stop_nats if stop_nats else max(1.0, STOP_REMAIN * summ[-1]["nsite"])
    objs = [s["obj"] for s in summ]
    gains = [objs[i] - objs[i - 1] for i in range(1, len(objs))]
    ok = 0
    for n in range(len(gains) - STOP_ROUNDS + 1, len(gains) + 1):
        if n < 6:
            return False
        lam, rem = projected(gains[:n])
        ok += lam is not None and lam < STOP_LAMBDA and rem < limit
    return ok == STOP_ROUNDS


# ---------------- starts ----------------

def cold_starts(data, kind, seed, K):
    rng = np.random.RandomState(seed)
    counts = []
    for i in data.ids:
        _, seqs = read_alignment(data.aln(i))
        counts.append(site_counts(seqs))
    C = np.vstack(counts)
    pbar = C.sum(axis=0) / C.sum()
    if kind == "dirichlet":
        profiles = [0.95 * rng.dirichlet(10 * pbar) + 0.05 * pbar for _ in range(K)]
    else:   # k-means++ seeding with KL distance on smoothed site compositions
        X = C[(C > 0).sum(axis=1) >= 3]
        X = (X + 2 * pbar) / (X.sum(axis=1, keepdims=True) + 2)
        centers = [X[rng.randint(len(X))]]
        dmin = np.full(len(X), np.inf)
        for _ in range(K - 1):
            c = centers[-1]
            dmin = np.minimum(dmin, np.sum(X * np.log(X / c), axis=1))
            centers.append(X[rng.choice(len(X), p=dmin / dmin.sum())])
        profiles = centers
    return [p / p.sum() for p in profiles], np.full(K, 1.0 / K)


def jitter(profiles, weights, seed):
    rng = np.random.RandomState(seed)
    ps = [p * np.exp(0.3 * rng.standard_normal(len(p))) for p in profiles]
    return [p / p.sum() for p in ps], rng.dirichlet(200 * np.asarray(weights) / np.sum(weights))


# ---------------- train ----------------

def preflight(data, binary, nproc, rundir):
    """one pass at the published C10 over all alignments; failures -> DATA/exclude.txt (once)"""
    excl = os.path.join(data.path, "exclude.txt")
    if os.path.exists(excl):
        return
    ids = read_list(os.path.join(data.path, "all.txt"))
    pdir = os.path.join(rundir, "preflight")
    os.makedirs(pdir, exist_ok=True)
    profiles, w = published_c10()
    rates = paper_rates(profiles)
    mdef = os.path.join(pdir, "c10pf.nex")
    write_mdef(mdef, NAME, profiles, w, rates)
    check = lambda d: echo_matches(d, profiles, w, rates)
    jobs = [dict(bin=binary, aln=data.aln(i), tree=data.tree(i), alpha=None, mdef=mdef, model=NAME,
                 dir=os.path.join(pdir, i), check=check, size=data.size[i]) for i in ids]
    log("preflight: %d alignments at the published C10" % len(jobs))
    res = run_jobs(jobs, nproc)
    bad = [(i, e) for i, (d, e) in zip(ids, res) if d is None]
    with open(excl + ".tmp", "w") as f:
        for i, e in bad:
            f.write("%s\t%s\n" % (i, e))
    os.replace(excl + ".tmp", excl)
    log("preflight: %d excluded (DATA/exclude.txt)" % len(bad))
    if len(bad) > 0.005 * len(ids):
        sys.exit("more than 0.5%% of alignments failed the preflight; see %s" % excl)
    shutil.rmtree(pdir)


def init_chains(args, data):
    K = args.K
    pub_p, pub_w = published_c10()
    names = []
    if args.option == "weights":
        specs = [("weights1", "c10pf", 0)]
    elif args.option == "warm":
        n = args.chains or 3
        specs = [("warm1", "c10pf", 0)] + [("warm%d" % (k + 1), "jitter", args.seed + k) for k in range(1, n)]
    else:
        n = args.chains or 12
        specs = [("cold_km%d" % (k + 1), "kmeans", args.seed + k) for k in range(n // 2)]
        specs += [("cold_dir%d" % (k + 1), "dirichlet", args.seed + 100 + k) for k in range(n - n // 2)]
    for name, start, seed in specs:
        names.append(name)
        if os.path.exists(os.path.join(args.rundir, name, "chain.json")):
            continue
        if start == "c10pf":
            p, w = pub_p, pub_w
        elif start == "jitter":
            p, w = jitter(pub_p, pub_w, seed)
        else:
            p, w = cold_starts(data, start, seed, K)
        cfg = dict(start=start, seed=seed, mode="weights" if args.option == "weights" else "full")
        new_chain(args.rundir, name, cfg, p, w, data)
    return names


def reduce_chain(ch, r, dumps, data):
    """pool round r, check, M-step, write round r+1 (model.nex last); returns new status or None"""
    name, profiles, weights, rates = ch.model(r)
    P = M.pool(dumps)
    lp = M.log_prior(profiles, weights)
    summ = ch.summaries()
    obj, obj_init = P["logl"] + lp, P["logl_init"] + lp
    warn = sum(d["logl"] < d["logl_init"] - max(1e-3, 1e-3 * abs(d["logl"])) for d in dumps)
    s = dict(round=r, nsite=P["nsite"], logl=P["logl"], logl_init=P["logl_init"], log_prior=lp, obj=obj, obj_init=obj_init,
             gain=obj - summ[-1]["obj"] if summ else None,
             mstep_gain=obj_init - summ[-1]["obj"] if summ else None, job_warnings=int(warn))
    status = None
    if summ and obj_init < summ[-1]["obj"] - MONO_TOL:
        status = "failed: EM inequality violated at round %d (%.4f < %.4f)" % (r, obj_init, summ[-1]["obj"])
    new_p, new_w, new_rates, c = M.mstep(P, profiles, weights, ch.cfg["mode"])
    s["scale"] = c
    s["weights"] = [float(x) for x in np.asarray(weights) / np.sum(weights)]
    nxt = ch.rdir(r + 1)
    os.makedirs(nxt, exist_ok=True)
    with open(os.path.join(nxt, "inputs.tsv"), "w") as f:
        for i, d in zip(data.ids, dumps):
            f.write("%s\t%r\t%s\n" % (i, d["alpha"], scale_tree(d["tree"], c)))
    with open(os.path.join(ch.rdir(r), "summary.json"), "w") as f:
        json.dump(s, f, indent=1)
    if status:
        return status
    write_mdef(os.path.join(nxt, "model.nex"), NAME, new_p, new_w, new_rates)
    shutil.rmtree(os.path.join(ch.rdir(r), "jobs"), ignore_errors=True)
    return None


def halve(rundir, chains, r, keep):
    """keep the best `keep` live chains by projected final objective (+ margin), one per basin"""
    path = os.path.join(rundir, "halving.json")
    plan = json.load(open(path)) if os.path.exists(path) else {}
    if str(r) not in plan:
        live = [c for c in chains if c.status() == "running"]
        scored = []
        for c in live:
            objs = [s["obj"] for s in c.summaries()]
            gains = [objs[i] - objs[i - 1] for i in range(1, len(objs))]
            lam, _ = projected(gains)
            lc = min(max(lam if lam is not None else 0.0, 0.0), 0.95)
            scored.append((objs[-1] + (gains[-1] * lc / (1 - lc) if gains else 0.0), c))
        scored.sort(key=lambda x: -x[0])
        kept, cutoff = [], None
        for v, c in scored:
            p_c = c.model(c.current_round())[1]
            if any(max(x[2] for x in match_profiles(p_c, k.model(k.current_round())[1])[0]) < SAME_BASIN
                   for _, k in kept):
                continue    # same basin as a better chain: promote the next one instead
            if len(kept) < keep or v >= cutoff - HALVE_MARGIN:
                kept.append((v, c))
                if len(kept) == keep:
                    cutoff = v
        plan[str(r)] = {"kept": [c.name for _, c in kept], "projected": {c.name: v for v, c in scored}}
        with open(path + ".tmp", "w") as f:
            json.dump(plan, f, indent=1)
        os.replace(path + ".tmp", path)
    for c in chains:
        if c.status() == "running" and c.name not in plan[str(r)]["kept"]:
            c.set_status("pruned at round %d" % r)
    log("round %d: kept %s" % (r, ", ".join(plan[str(r)]["kept"])))


def select(rundir, chains):
    done = [c for c in chains if c.status() in ("converged", "not converged")]
    os.makedirs(os.path.join(rundir, "finalists"), exist_ok=True)
    rows = []
    for c in chains:
        summ = c.summaries()
        rows.append((c.name, c.cfg["start"], c.status(), len(summ), summ[-1]["obj"] if summ else float("nan")))
    for c in done:
        shutil.copy(os.path.join(c.rdir(c.current_round()), "model.nex"),
                    os.path.join(rundir, "finalists", c.name + ".nex"))
    if done:
        best = max(done, key=lambda c: c.summaries()[-1]["obj"])
        shutil.copy(os.path.join(best.rdir(best.current_round()), "model.nex"), os.path.join(rundir, "best.nex"))
        log("best chain: %s" % best.name)
    with open(os.path.join(rundir, "summary.tsv"), "w") as f:
        f.write("chain\tstart\tstatus\trounds\tfinal_penalised_logl\n")
        for row in rows:
            f.write("%s\t%s\t%s\t%d\t%.6f\n" % row)


def train(args):
    os.makedirs(args.rundir, exist_ok=True)
    lockf = open(os.path.join(args.rundir, ".lock"), "w")
    try:
        fcntl.flock(lockf, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        sys.exit("another em.py is running in %s" % args.rundir)
    binary = os.path.abspath(args.bin)
    data = Data(args.data)
    if not os.path.exists(os.path.join(data.path, "exclude.txt")):
        if args.dry_run:
            log("dry run: would run the preflight first")
        else:
            preflight(data, binary, args.j, args.rundir)
            data = Data(args.data)
    halving = {}
    if args.halving != "none" and args.option == "cold":
        halving = {int(a): int(b) for a, b in (x.split(":") for x in args.halving.split(","))}
    names = init_chains(args, data) if not args.dry_run else []
    cfg_path = os.path.join(args.rundir, "config.json")
    if not os.path.exists(cfg_path):
        with open(cfg_path, "w") as f:
            json.dump(dict(option=args.option, chains=names, halving=halving, data=data.path,
                           train_list_sha256=sha256(data.list_path), n_train=len(data.ids),
                           binary=binary, binary_sha256=sha256(binary), seed=args.seed,
                           max_rounds=args.max_rounds), f, indent=1)
    if args.dry_run:
        log("dry run: option %s, %d training alignments, %d jobs per chain-round, J=%d"
            % (args.option, len(data.ids), len(data.ids), args.j))
        return
    chains = [Chain(args.rundir, n) for n in names]
    while True:
        live = [c for c in chains if c.status() == "running"]
        if not live:
            break
        r = min(c.current_round() for c in live)
        step = [c for c in live if c.current_round() == r]
        jobs, owners = [], []
        for c in step:
            name, profiles, weights, rates = c.model(r)
            check = (lambda p, w, q: (lambda d: echo_matches(d, p, w, q)))(profiles, weights, rates)
            mdef = os.path.join(c.rdir(r), "model.nex")
            for i, (alpha, tree) in c.inputs(r).items():
                jobs.append(dict(bin=binary, aln=data.aln(i), tree=tree, alpha=alpha, mdef=mdef, model=NAME,
                                 dir=os.path.join(c.rdir(r), "jobs", i), check=check, size=data.size[i]))
                owners.append(c)
        t0 = time.time()
        res = run_jobs(jobs, args.j)
        for c in step:
            mine = [(j, x) for j, x, o in zip(jobs, res, owners) if o is c]
            fails = [(os.path.basename(j["dir"]), x[1]) for j, x in mine if x[0] is None]
            if fails:
                c.set_status("failed: %d jobs in round %d, e.g. %s: %s" % (len(fails), r, fails[0][0], fails[0][1]))
                continue
            st = reduce_chain(c, r, [x[0] for _, x in mine], data)
            if st:
                c.set_status(st)
            elif converged(c.summaries(), args.stop_nats):
                c.set_status("converged")
            elif r + 1 >= args.max_rounds:
                c.set_status("not converged")
            s = c.summaries()[-1]
            log("%s round %d: logl %.4f gain %s (%.0f s)" % (c.name, r, s["logl"],
                "%.4f" % s["gain"] if s["gain"] is not None else "-", time.time() - t0))
        if r + 1 in halving:
            halve(args.rundir, [c for c in chains if c.cfg["start"] in ("kmeans", "dirichlet")], r + 1, halving[r + 1])
    select(args.rundir, chains)


# ---------------- status / score ----------------

def status(args):
    cfg = json.load(open(os.path.join(args.rundir, "config.json")))
    print("%-12s %-30s %6s %18s %12s %10s" % ("chain", "status", "round", "penalised logl", "last gain", "remaining"))
    for n in cfg["chains"]:
        c = Chain(args.rundir, n)
        summ = c.summaries()
        objs = [s["obj"] for s in summ]
        gains = [objs[i] - objs[i - 1] for i in range(1, len(objs))]
        _, rem = projected(gains)
        print("%-12s %-30s %6d %18.4f %12s %10s" % (n, c.status()[:30], len(summ), objs[-1] if objs else float("nan"),
              "%.4f" % gains[-1] if gains else "-", "%.1f" % rem if rem != float("inf") else "-"))


def score(args):
    data = Data(args.data, args.list)
    os.makedirs(args.out, exist_ok=True)
    binary = os.path.abspath(args.bin)

    if args.wsl:
        os.makedirs(os.path.join(args.out, "sitelh"), exist_ok=True)

    def one(i):
        d = os.path.join(args.out, "runs", i)
        os.makedirs(d, exist_ok=True)
        cmd = [binary, "-s", data.aln(i), "-m", args.m, "-seed", "1",
               "-keep-ident", "-nt", "1", "-quiet", "-redo", "-pre", os.path.join(d, "s")]
        if not args.search:     # fixed topology: the DATA tree; otherwise IQ-TREE's own tree search
            with open(os.path.join(d, "in.tree"), "w") as f:
                f.write(data.tree(i) + "\n")
            cmd += ["-te", os.path.join(d, "in.tree")]
        if args.mdef:
            cmd += ["-mdef", os.path.abspath(args.mdef)]
        if args.opt_weights:    # re-optimise the mixture weights on this alignment
            cmd += ["-mwopt"]
        if args.wsl:
            cmd += ["-wsl"]
        try:
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=args.timeout,
                           env=dict(os.environ, OMP_NUM_THREADS="1"))
        except subprocess.TimeoutExpired:
            pass
        if args.wsl and os.path.exists(os.path.join(d, "s.sitelh")):
            shutil.copy(os.path.join(d, "s.sitelh"), os.path.join(args.out, "sitelh", i + ".sitelh"))
        lnl = alpha = None
        try:
            for l in open(os.path.join(d, "s.iqtree")):
                if l.startswith("Log-likelihood of the tree:"):
                    lnl = float(l.split(":")[1].split()[0])
                if "Gamma shape alpha:" in l:
                    alpha = float(l.split(":")[1].split()[0])
        except OSError:
            pass
        return i, lnl, alpha

    with ThreadPoolExecutor(max_workers=args.j) as pool:
        res = list(pool.map(one, sorted(data.ids, key=lambda i: -data.size[i])))
    res.sort()
    nsite = {}
    with open(os.path.join(data.path, "index.tsv")) as f:
        next(f)
        for l in f:
            v = l.split("\t")
            nsite[v[0]] = int(v[2])
    with open(os.path.join(args.out, "scores.tsv"), "w") as f:
        f.write("id\tnsite\tlogl\talpha\n")
        for i, lnl, a in res:
            f.write("%s\t%d\t%s\t%s\n" % (i, nsite[i], "" if lnl is None else repr(lnl), "" if a is None else repr(a)))
    with open(os.path.join(args.out, "score.json"), "w") as f:
        json.dump(dict(model=args.m, mdef=args.mdef and os.path.abspath(args.mdef), search=args.search,
                       opt_weights=args.opt_weights, list=data.list_path, binary=binary), f, indent=1)
    shutil.rmtree(os.path.join(args.out, "runs"), ignore_errors=True)
    print("scored %d alignments (%d failed) -> %s" % (len(res), sum(r[1] is None for r in res),
                                                    os.path.join(args.out, "scores.tsv")))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("train")
    t.add_argument("rundir")
    t.add_argument("--data", required=True)
    t.add_argument("--bin", required=True)
    t.add_argument("--option", choices=["warm", "cold", "weights"], default="warm")
    t.add_argument("-j", type=int, default=24)
    t.add_argument("--chains", type=int, default=0)
    t.add_argument("--halving", default="30:6,60:3")
    t.add_argument("--max-rounds", type=int, default=300)
    t.add_argument("--seed", type=int, default=1)
    t.add_argument("--K", type=int, default=10)
    t.add_argument("--stop-nats", type=float, help="absolute projected-gain stop (simulations: 1)")
    t.add_argument("--dry-run", action="store_true")
    s = sub.add_parser("status")
    s.add_argument("rundir")
    c = sub.add_parser("score")
    c.add_argument("out")
    c.add_argument("--data", required=True)
    c.add_argument("--list", required=True)
    c.add_argument("--bin", required=True)
    c.add_argument("-j", type=int, default=24)
    c.add_argument("--mdef")
    c.add_argument("-m", required=True)
    c.add_argument("--search", action="store_true", help="optimise the topology (default: fixed DATA tree)")
    c.add_argument("--opt-weights", action="store_true", help="re-optimise mixture weights per alignment (-mwopt)")
    c.add_argument("--wsl", action="store_true", help="keep site log-likelihoods in OUT/sitelh/")
    c.add_argument("--timeout", type=int, default=4 * 3600)
    args = ap.parse_args()
    {"train": train, "status": status, "score": score}[args.cmd](args)


if __name__ == "__main__":
    main()
