#!/usr/bin/env python3
"""
Benchmark of --analytical-gradients against the default optimiser.

Arms separate the starting point from the optimiser (design doc, section
13), so speed and quality are not conflated:

  old-default   default optimiser, IQ-TREE's usual start
  old-init2     default optimiser, -init_nucl_freq 2 (distinct profile starts)   [--full]
  old-c10warm   default optimiser started from C10 profiles (-mfopt)            [--full, LG+Fk only]
  old-em        default optimiser with -optalg_qmix EM                          [--full]
  new-warm      --analytical-gradients (C-series warm start, EM axes, polish)
  new-cold      --analytical-gradients --ag-start cold
  new-multi     --analytical-gradients --ag-multistart -1 (automatic count)   [--full]

Datasets are AliSim simulations with recorded truth (random Yule-Harding
trees, seeds from --seeds) plus the repository's real alignments. Every
run records the final log-likelihood, wall time, peak memory, tree length,
the optimiser's evaluation counts and, for simulated mixtures, the RMSE of
the recovered weights and profiles against the truth after matching the
classes.

Usage:
  bench.py <iqtree_binary> <out_dir> [--quick | --full | --thorough | --mini] [--big] [-j N]
           [--threads 1,8] [--repeats N] [--seeds 101,102] [--only name,name] [--dry-run]

--thorough: LG+F10, GTR20+F12 and GTR20+C60 (linked exchangeabilities, with
and without -mwopt) under +R8, +I+G4 and +I+R10 variants, 32 taxa x 3000 sites
(C60: 24 x 2000), three replicates, on the true tree, every analytic setting
(warm/cold start, cascade, multi-start, without EM) against the default
optimiser. The F10/F12 profiles and weights are drawn at random per seed
(Dirichlet around the LG/WAG frequencies) so that no arm starts at the truth;
the -mwopt C60 alignments are simulated with random C60 weights, the
fixed-weight ones with the built-in weights. Arms that only differ in the
profile start (cold, multi-start) are skipped for fixed-profile C-series
models, where they are identical to new-warm. -j defaults to 50% of the
cores.

--mini: LG+F2+I+G4 (16 taxa x 1500 sites) and WAG+F4+R4 (18 taxa x 2000
sites), two seeds, small enough that the default (numerical) optimiser
actually finishes in reasonable time; same truth methodology and arms as
--thorough. -j defaults to 25% of the cores (a --thorough run holds 50%, so
both together stay within 75% of the machine).

In the --thorough and --mini tiers every simulated dataset also gets a
`truth` row: the log-likelihood of the simulation model itself on the true
tree with fixed branch lengths (-blfix), i.e. no optimisation at all. It is
a reference for the lnL columns, not an optimiser arm (no speed-up, excluded
from the "best arm" used for dlnL).

Outputs in <out_dir>: results.tsv, summary.md, manifest.json, plots (if
matplotlib is available). Re-generate the report with bench_report.py, or
re-parse the run directories with `bench.py --rescore <out_dir>`.
"""
import sys
import os
import re
import json
import math
import time
import random
import shutil
import socket
import itertools
import subprocess
from concurrent.futures import ThreadPoolExecutor

RUN_TIMEOUT = 12 * 3600   # seconds per run; a run past it is recorded with status exit124
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import bench_report  # noqa: E402

P1 = "0.18/0.1/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.02/0.02/0.02/0.02/0.02/0.02"
P2 = "0.02/0.02/0.02/0.02/0.02/0.02/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.1/0.18"
P3 = "/".join(["0.05"] * 20)
P4 = "0.1/0.1/0.1/0.1/0.02/0.02/0.02/0.02/0.02/0.02/0.02/0.02/0.02/0.02/0.02/0.02/0.02/0.02/0.1/0.1"
R4 = "+R4{0.2,0.3,0.3,0.8,0.3,1.2,0.2,2.5}"
R8 = "+R8{0.1,0.1,0.1,0.3,0.1,0.6,0.15,0.9,0.15,1.2,0.15,1.6,0.1,2.5,0.05,4}"

# equilibrium frequencies of the LG and WAG matrices (order ARNDCQEGHILKMFPSTWYV), the centres of the random profiles
BASE_FREQ = {
    "LG":  [0.079066, 0.055941, 0.041977, 0.053052, 0.012937, 0.040767, 0.071586, 0.057337, 0.022355, 0.062157,
            0.099081, 0.064600, 0.022951, 0.042302, 0.044040, 0.061197, 0.053287, 0.012066, 0.034155, 0.069147],
    "WAG": [0.0866279, 0.043972, 0.0390894, 0.0570451, 0.0193078, 0.0367281, 0.0580589, 0.0832518, 0.0244313, 0.048466,
            0.086209, 0.0620286, 0.0195027, 0.0384319, 0.0457631, 0.0695179, 0.0610127, 0.0143859, 0.0352742, 0.0708956],
}


def mix4(base, weights):
    return "MIX{%s}" % ",".join("%s+F{%s}:1:%s" % (base, p, w) for p, w in zip((P1, P2, P3, P4), weights))


def cseries(k):
    """profiles (name -> list) and weights of the built-in C<k> mixture, read from model/modelmixture.cpp"""
    src = open(os.path.join(ROOT, "model", "modelmixture.cpp")).read()
    prof = {m.group(1): [float(x) for x in m.group(2).split()]
            for m in re.finditer(r"^frequency (C%dpi\d+) = ([0-9. ]+);" % k, src, re.M)}
    line = re.search(r"^model C%d = [^;]+;" % k, src, re.M).group(0)
    weights = {n: float(w) for n, w in re.findall(r"(C%dpi\d+):1:([0-9.]+)" % k, line)}
    return prof, weights


def explicit_mix(base, comps):
    """MIX{base+F{p1/p2/..}:1:w,...} for a list of (weight, profile)"""
    return "MIX{%s}" % ",".join("%s+F{%s}:1:%.6f" % (base, "/".join("%.6f" % v for v in p), w) for w, p in comps)


def cmix(base, k, n):
    """MIX of the first n profiles of C<k> under exchangeabilities `base`, weights renormalised"""
    prof, w = cseries(k)
    names = ["C%dpi%d" % (k, i + 1) for i in range(n)]
    tot = sum(w[x] for x in names)
    return explicit_mix(base, [(w[x] / tot, prof[x]) for x in names])


def dirichlet(rng, alphas, floor):
    """one Dirichlet draw (gamma normalisation), entries floored and renormalised"""
    x = [rng.gammavariate(a, 1.0) for a in alphas]
    s = sum(x)
    x = [max(v / s, floor) for v in x]
    s = sum(x)
    return [v / s for v in x]


def random_profile_mix(base, k, seed):
    """k random profiles around the base frequencies (Dirichlet, concentration 8, floor 1e-3) with random
    weights (Dirichlet(3), floor 0.02); the draw depends on (base, k, seed) only, so the rate variants
    of one family share the truth at a given seed"""
    rng = random.Random("profiles-%s-%d-%d" % (base, k, seed))
    comps = [dirichlet(rng, [8.0 * f for f in BASE_FREQ[base]], 1e-3) for _ in range(k)]
    weights = dirichlet(rng, [3.0] * k, 0.02)
    return explicit_mix(base, list(zip(weights, comps)))


def random_c60_mix(base, seed):
    """the 60 C60 profiles under `base` with random weights (Dirichlet(1), floor 0.002)"""
    prof, _ = cseries(60)
    rng = random.Random("c60-weights-%s-%d" % (base, seed))
    w = dirichlet(rng, [1.0] * 60, 0.002)
    return explicit_mix(base, [(w[i], prof["C60pi%d" % (i + 1)]) for i in range(60)])


IG = "+I{0.15}+G4{0.5}"
IR10 = "+I{0.1}+R10{0.2,0.05,0.16,0.15,0.14,0.3,0.12,0.5,0.1,0.75,0.08,1.0,0.07,1.4,0.06,1.9,0.04,2.6,0.03,3.8}"

# name -> dict(sim=(alisim model | callable(seed) -> model, ntaxa, length) | real=path | alias_of=name,
#              model=iqtree model args, tier)
# tier "thorough": three model families with several rate-heterogeneity variants (32 taxa x 3000 sites)
C60_INF = "--gtr20-model LG --link-exchange-rates"
DATASETS = {
    "th_lg_f10_r8":       dict(sim=(lambda s: random_profile_mix("LG", 10, s) + R8, 32, 3000), model="LG+F10+R8", tier="thorough"),
    "th_lg_f10_ig":       dict(sim=(lambda s: random_profile_mix("LG", 10, s) + IG, 32, 3000), model="LG+F10+I+G4", tier="thorough"),
    "th_lg_f10_ir10":     dict(sim=(lambda s: random_profile_mix("LG", 10, s) + IR10, 32, 3000), model="LG+F10+I+R10", tier="thorough"),
    "th_gtr20_f12_r8":    dict(sim=(lambda s: random_profile_mix("WAG", 12, s) + R8, 32, 3000), model="GTR20+F12+R8 --gtr20-model LG", tier="thorough"),
    "th_gtr20_f12_ig":    dict(sim=(lambda s: random_profile_mix("WAG", 12, s) + IG, 32, 3000), model="GTR20+F12+I+G4 --gtr20-model LG", tier="thorough"),
    "th_gtr20_c60_r8":    dict(sim=("WAG+C60" + R8, 24, 2000), model="GTR20+C60+R8 " + C60_INF, tier="thorough"),
    "th_gtr20_c60_r8_mwopt": dict(sim=(lambda s: random_c60_mix("WAG", s) + R8, 24, 2000), model="GTR20+C60+R8 " + C60_INF + " -mwopt", tier="thorough"),
    "th_gtr20_c60_ig":    dict(sim=("WAG+C60" + IG, 24, 2000), model="GTR20+C60+I+G4 " + C60_INF, tier="thorough"),
    "th_gtr20_c60_ig_mwopt": dict(sim=(lambda s: random_c60_mix("WAG", s) + IG, 24, 2000), model="GTR20+C60+I+G4 " + C60_INF + " -mwopt", tier="thorough"),
    # mini tier: small alignments so the default (numerical) optimiser actually
    # finishes, giving a real old-vs-new speed/quality comparison rather than a
    # timeout; 2- and 4-profile mixtures, random Dirichlet truth per seed like
    # the thorough tier's F10/F12 datasets.
    "mini_lg_f2_ig":  dict(sim=(lambda s: random_profile_mix("LG", 2, s) + IG, 16, 1500), model="LG+F2+I+G4", tier="mini"),
    "mini_wag_f4_r4": dict(sim=(lambda s: random_profile_mix("WAG", 4, s) + R4, 18, 2000), model="WAG+F4+R4", tier="mini"),
    # GTR20 needs more sites than LG/WAG for its 189 exchangeabilities to be informed; simulated under WAG so the LG start is not the truth
    "mini_gtr20_f2_r4": dict(sim=(lambda s: random_profile_mix("WAG", 2, s) + R4, 20, 2500), model="GTR20+F2+R4 --gtr20-model LG", tier="mini"),
    "sim_dna_gtr_g4":   dict(sim=("GTR{1.5,3,0.8,1.2,2.5}+F{0.3,0.2,0.2,0.3}+G4{0.7}", 20, 5000), model="GTR+FO+G4", tier="quick"),
    "sim_lg_f4_r4":     dict(sim=(mix4("LG", (0.25, 0.25, 0.25, 0.25)) + R4, 20, 4000), model="LG+F4+R4", tier="quick"),
    "real_aa_example_lg_f2_g4": dict(real="example/aa_example.phy", model="LG+F2+G4", tier="quick"),
    "sim_wag_f4_r8":    dict(sim=(mix4("WAG", (0.4, 0.3, 0.2, 0.1)) + R8, 20, 4000), model="WAG+F4+R8", tier="full"),
    "sim_lg_c10_r8":    dict(sim=("LG+C10" + R8, 50, 5000), model="LG+C10+R8", tier="full"),
    "sim_gtr20_c10_r8": dict(sim=("LG+C10" + R8, 40, 5000), model="GTR20+C10+R8 --gtr20-model LG", tier="full"),
    "real_turtle_lg_f4_r4":  dict(real="test_scripts/test_data/turtle_aa.fasta", model="LG+F4+R4", tier="full"),
    "real_turtle_lg_c10_r4": dict(real="test_scripts/test_data/turtle_aa.fasta", model="LG+C10+R4", tier="full"),
    "real_m126_lg_f10_r8":   dict(real="test_scripts/test_data/prot_M126_27_269.phy", model="LG+F10+R8", tier="full"),
    "real_m126_gtr20_f10_r8": dict(real="test_scripts/test_data/prot_M126_27_269.phy", model="GTR20+F10+R8 --gtr20-model LG", tier="full"),
    "sim_gtr20_c10_r8_big": dict(sim=("LG+C10" + R8, 100, 20000), model="GTR20+C10+R8 --gtr20-model LG", tier="big"),
    "sim_gtr20_f60_r8":     dict(sim=("LG+C60" + R8, 50, 20000), model="GTR20+F60+R8 --gtr20-model LG", tier="big"),
}

ARMS = {
    "old-default": ("", "quick"),
    "new-warm":    ("--analytical-gradients --ag-stats", "quick"),
    "new-cold":    ("--analytical-gradients --ag-stats --ag-start cold", "quick"),
    "old-init2":   ("-init_nucl_freq 2", "full"),
    "old-em":      ("-optalg_qmix EM", "full"),
    "old-c10warm": ("__C10WARM__", "full"),
    "new-multi":   ("--analytical-gradients --ag-stats --ag-multistart -1", "full"),
    # thorough tier: every analytic setting against the default optimiser
    "new-cascade":    ("--analytical-gradients --ag-stats --ag-cascade on", "thorough"),
    "new-noem":       ("--analytical-gradients --ag-stats --ag-em-axes none", "thorough"),   # "none": no W/R/F letters -> no EM step
    "new-cold-multi": ("--analytical-gradients --ag-stats --ag-start cold --ag-multistart -1", "thorough"),
}
THOROUGH_ARMS = ["old-default", "new-warm", "new-cascade", "new-cold", "new-multi", "new-noem", "new-cold-multi"]
MINI_ARMS = THOROUGH_ARMS   # same set; the point of this tier is small enough data that old-default finishes
# arms that only change the start of estimated (+FO) profiles; identical to new-warm when the profiles are fixed
PROFILE_START_ARMS = ("new-cold", "new-multi", "new-cold-multi")


def sim_model(d, seed):
    """the AliSim model string of a simulated dataset at a seed"""
    m = d["sim"][0]
    return m(seed) if callable(m) else m


def c10warm_model(model):
    """LG+Fk+Rn  ->  MIX{LG+FC10pi1,...,LG+FC10pik}+Rn with -mfopt (k <= 10), else None"""
    m = re.match(r"^([A-Za-z0-9]+)\+F(\d+)(\+.*)?$", model)
    if not m or int(m.group(2)) > 10:
        return None
    base, k, rest = m.group(1), int(m.group(2)), m.group(3) or ""
    return "MIX{%s}%s" % (",".join("%s+FC10pi%d" % (base, i + 1) for i in range(k)), rest)


def truth_from_alisim(model):
    """weights and profiles of an explicit MIX{...+F{...}:1:w,...}[+R..] or <base>+C<k>[+..] simulation model"""
    mc = re.match(r"^[A-Za-z0-9]+\+C(\d+)(\+.*)?$", model)
    if mc:
        k = int(mc.group(1))
        prof, w = cseries(k)
        return [(w["C%dpi%d" % (k, i + 1)], prof["C%dpi%d" % (k, i + 1)]) for i in range(k)]
    if not model.startswith("MIX{"):
        return None
    depth, end = 0, -1
    for i, ch in enumerate(model):   # the matching brace of MIX{ (profiles and +R{} have braces of their own)
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                end = i
                break
    if end < 0:
        return None
    comps = []
    for part in model[4:end].split(","):
        mm = re.match(r".*\+F\{([^}]*)\}:1:([0-9.]+)$", part)
        if not mm:
            return None
        comps.append((float(mm.group(2)), [float(x) for x in mm.group(1).split("/")]))
    return comps


def read_truth(prefix):
    """the recorded truth of a simulated dataset (<prefix>.truth.json), or None"""
    path = prefix + ".truth.json"
    if not os.path.exists(path):
        return None
    t = json.load(open(path))
    return [(c[0], c[1]) for c in t["truth"]] if t.get("truth") is not None else None


def parse_iqtree_report(path):
    out = {"logl": float("nan"), "tree_length": float("nan"), "comps": []}
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith("Log-likelihood of the tree:"):
            out["logl"] = float(re.search(r"tree: (-?[0-9.]+)", line).group(1))
        elif line.startswith("Total tree length"):
            out["tree_length"] = float(re.search(r": ([0-9.]+)", line).group(1))
        else:
            mm = re.match(r"\s*\d+\s+\S+\s+([0-9.]+)\s+([0-9.]+)\s+(\S+)", line)
            if mm and re.search(r"\+F", mm.group(3)):
                mp = re.search(r"\+FO?\{([^}]*)\}", mm.group(3))
                out["comps"].append((float(mm.group(2)), [float(x) for x in mp.group(1).split(",")] if mp else None))
    return out


def rmse_vs_truth(comps, truth):
    """match classes (all permutations up to 8 classes, greedy beyond) and return (rmse weights, rmse profiles)"""
    if not comps or not truth or len(comps) != len(truth):
        return float("nan"), float("nan")
    k = len(truth)
    if any(c[1] is None for c in comps):
        # fixed profiles (C-series): classes keep their order, only the weights are estimated
        return math.sqrt(sum((comps[i][0] - truth[i][0]) ** 2 for i in range(k)) / k), float("nan")

    def cost(perm):
        return sum(sum((a - b) ** 2 for a, b in zip(comps[perm[i]][1], truth[i][1])) for i in range(k))
    if k <= 8:
        best = min(itertools.permutations(range(k)), key=cost)
    else:
        best, used = [], set()
        for i in range(k):
            j = min((j for j in range(k) if j not in used), key=lambda j: sum((a - b) ** 2 for a, b in zip(comps[j][1], truth[i][1])))
            used.add(j)
            best.append(j)
    rw = math.sqrt(sum((comps[best[i]][0] - truth[i][0]) ** 2 for i in range(k)) / k)
    rp = math.sqrt(cost(best) / (k * len(truth[0][1])))
    return rw, rp


def parse_ag_stats(path):
    if not os.path.exists(path):
        return float("nan"), float("nan")
    lh = gr = 0
    for line in open(path):
        if line.startswith("AG stats:"):
            lh += int(re.search(r"likelihood_evaluations=(\d+)", line).group(1))
            gr += int(re.search(r"gradient_evaluations=(\d+)", line).group(1))
    return (lh, gr) if lh else (float("nan"), float("nan"))


def run(cmd, cwd, log):
    """run a command, return (exit, wall seconds, peak RSS MB)"""
    timebin = shutil.which("time") if os.path.exists("/usr/bin/time") else None
    tfile = os.path.join(cwd, "time.txt")
    if timebin:
        full = [timebin, "-f", "%e %M", "-o", tfile] + cmd
    else:
        full = cmd
    t0 = time.time()
    with open(log, "w") as lf:
        # own process group, so a timeout kills IQ-TREE itself and not only the time wrapper
        proc = subprocess.Popen(full, cwd=cwd, stdout=lf, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            rc = proc.wait(timeout=RUN_TIMEOUT)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(proc.pid, 9)
            except OSError:
                pass
            proc.wait()
            rc = 124
    wall, rss = time.time() - t0, float("nan")
    if timebin and os.path.exists(tfile):
        try:
            w, m = open(tfile).read().split()[-2:]
            wall, rss = float(w), float(m) / 1024.0
        except ValueError:
            pass
    return rc, wall, rss


def rescore(out_dir):
    """bench.py --rescore <out_dir>: re-parse every run directory (report, stats, RMSE) and rewrite results.tsv/summary.md"""
    results_path = os.path.join(out_dir, "results.tsv")
    rows = bench_report.read_results(results_path)
    for r in rows:
        rdir = os.path.join(out_dir, "runs", "%s_s%s_%s_t%s_%s_r%s" % (r["dataset"], r["seed"], r["tree_mode"], r["threads"], r["arm"], r["repeat"]))
        rep_ = parse_iqtree_report(os.path.join(rdir, "run.iqtree"))
        lh, gr = parse_ag_stats(os.path.join(rdir, "run.stdout"))
        d = DATASETS.get(r["dataset"], {})
        truth = None
        if "sim" in d:
            prefix = os.path.join(out_dir, "data", "%s_s%s" % (r["dataset"], r["seed"]))
            truth = read_truth(prefix)
            if truth is None:
                truth = truth_from_alisim(sim_model(d, int(r["seed"])))
        rw, rp = rmse_vs_truth(rep_["comps"], truth)
        r.update(logl="%.4f" % rep_["logl"], tree_length="%.4f" % rep_["tree_length"], lh_evals=lh, grad_evals=gr,
                 rmse_weights="%.5f" % rw, rmse_profiles="%.5f" % rp)
    header = ["dataset", "seed", "tree_mode", "threads", "arm", "repeat", "status", "logl", "wall_s", "peak_rss_mb",
              "tree_length", "lh_evals", "grad_evals", "rmse_weights", "rmse_profiles", "command"]
    with open(results_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r[h]) for h in header) + "\n")
    manifest = None
    mpath = os.path.join(out_dir, "manifest.json")
    if os.path.exists(mpath):
        manifest = json.load(open(mpath))
    bench_report.write_summary(rows, os.path.join(out_dir, "summary.md"), manifest)
    print("rescored", len(rows), "runs")
    return 0


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    if argv[1] == "--rescore":
        return rescore(os.path.abspath(argv[2]))
    binary = os.path.abspath(argv[1])
    out_dir = os.path.abspath(argv[2])
    mode, big, jobs, threads, repeats, seeds, only, dry = "quick", False, 4, [1], 1, [101], None, False
    i = 3
    while i < len(argv):
        a = argv[i]
        if a == "--quick": mode = "quick"
        elif a == "--full": mode, threads, repeats = "full", [1, 8], 3
        elif a == "--thorough": mode, seeds, jobs = "thorough", [101, 102, 103], max(1, min(28, int((os.cpu_count() or 4) * 0.5)))
        elif a == "--mini": mode, seeds, jobs = "mini", [101, 102], max(1, min(14, int((os.cpu_count() or 4) * 0.25)))
        elif a == "--big": big = True
        elif a == "--dry-run": dry = True
        elif a == "-j": i += 1; jobs = int(argv[i])
        elif a == "--threads": i += 1; threads = [int(x) for x in argv[i].split(",")]
        elif a == "--repeats": i += 1; repeats = int(argv[i])
        elif a == "--seeds": i += 1; seeds = [int(x) for x in argv[i].split(",")]
        elif a == "--only": i += 1; only = set(argv[i].split(","))
        else:
            print("unknown argument", a); return 2
        i += 1
    os.makedirs(out_dir, exist_ok=True)
    tiers = {"quick"} if mode == "quick" else {"quick", "full"} if mode == "full" else {"thorough"} if mode == "thorough" else {"mini"}
    if big:
        tiers.add("big")
    datasets = {n: d for n, d in DATASETS.items() if d["tier"] in tiers and (only is None or n in only)}
    if mode in ("thorough", "mini"):
        arms = {n: ARMS[n][0] for n in (THOROUGH_ARMS if mode == "thorough" else MINI_ARMS)}
    else:
        arms = {n: a for n, (a, t) in ARMS.items() if t in tiers}

    # 1. simulate
    sim_dir = os.path.join(out_dir, "data")
    os.makedirs(sim_dir, exist_ok=True)
    inputs = {}   # (dataset, seed) -> (alignment, true tree or None, truth or None)
    for name, d in datasets.items():
        if "real" in d:
            inputs[(name, 0)] = (os.path.join(ROOT, d["real"]), None, None)
            continue
        if "alias_of" in d:
            continue
        _, ntaxa, length = d["sim"]
        for seed in seeds:
            model = sim_model(d, seed)
            prefix = os.path.join(sim_dir, "%s_s%d" % (name, seed))
            truth = truth_from_alisim(model)
            if dry:
                print("sim: %s seed %d  %d x %d  %s" % (name, seed, ntaxa, length, model[:100] + ("..." if len(model) > 100 else "")))
            elif not os.path.exists(prefix + ".phy"):
                rc, _, _ = run([binary, "--alisim", prefix, "-t", "RANDOM{yh,%d}" % ntaxa, "-rlen", "0.01", "0.1", "0.5",
                                "-m", model, "--length", str(length), "-seed", str(seed), "-redo", "-quiet"], sim_dir, prefix + ".alisim.log")
                if rc != 0:
                    print("simulation failed:", name, seed); return 1
            if not dry and not os.path.exists(prefix + ".truth.json"):
                with open(prefix + ".truth.json", "w") as f:
                    json.dump({"dataset": name, "seed": seed, "ntaxa": ntaxa, "length": length, "model": model, "truth": truth}, f)
            inputs[(name, seed)] = (prefix + ".phy", prefix + ".treefile", truth)
    for name, d in datasets.items():
        if "alias_of" in d:
            for seed in seeds:
                inputs[(name, seed)] = inputs[(d["alias_of"], seed)]

    # 2. the run list
    tasks = []
    for (name, seed), (aln, tree, truth) in inputs.items():
        d = datasets[name]
        modes = ["te"] if tree else []
        if not tree or mode == "full":
            modes.append("search")
        if tree and truth is not None and mode in ("thorough", "mini"):
            for thr in threads:
                tasks.append((name, seed, "te", thr, "truth", 0, aln, tree, truth, sim_model(d, seed), "-blfix"))
        for tmode in modes:
            for thr in threads:
                for arm, extra in arms.items():
                    model = d["model"]
                    if arm == "old-em" and not re.search(r"\+(F\d+|C\d+)", model.split()[0]):
                        continue   # -optalg_qmix EM only applies to mixtures
                    if arm in PROFILE_START_ARMS and not re.search(r"\+F\d+", model.split()[0]):
                        continue   # no estimated profiles: cold start / multi-start are no-ops, the run equals new-warm
                    if extra == "__C10WARM__":
                        cm = c10warm_model(model.split()[0])
                        if cm is None:
                            continue
                        extra = "-mfopt"
                        model = cm + model[len(model.split()[0]):]
                    for rep in range(repeats):
                        tasks.append((name, seed, tmode, thr, arm, rep, aln, tree, truth, model, extra))

    if dry:
        for t in tasks:
            print("run: %s s%d %s t%d %-15s -m %s %s" % (t[0], t[1], t[2], t[3], t[4], t[9], t[10]))
        print("%d runs, %d datasets, %d arms, %d parallel (dry run)" % (len(tasks), len(datasets), len(arms), jobs))
        return 0

    results_path = os.path.join(out_dir, "results.tsv")
    header = ["dataset", "seed", "tree_mode", "threads", "arm", "repeat", "status", "logl", "wall_s", "peak_rss_mb",
              "tree_length", "lh_evals", "grad_evals", "rmse_weights", "rmse_profiles", "command"]
    with open(results_path, "w") as f:
        f.write("\t".join(header) + "\n")
    manifest = {"binary": binary, "mode": mode, "repeats": repeats, "seeds": seeds, "threads": threads,
                "host": socket.gethostname(), "datasets": sorted(datasets), "arms": sorted(arms), "runs": len(tasks),
                "git": subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip() or "?"}
    with open(os.path.join(out_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    print("%d runs, %d datasets, %d arms, %d parallel" % (len(tasks), len(datasets), len(arms), jobs))

    def one(task):
        name, seed, tmode, thr, arm, rep, aln, tree, truth, model, extra = task
        rdir = os.path.join(out_dir, "runs", "%s_s%d_%s_t%d_%s_r%d" % (name, seed, tmode, thr, arm, rep))
        os.makedirs(rdir, exist_ok=True)
        cmd = [binary, "-s", aln, "-m"] + model.split() + ["-nt", str(thr), "-seed", str(1 + rep), "--prefix", "run", "-redo", "-lk", "FMA"]
        if tmode == "te":
            cmd += ["-te", tree]
        cmd += extra.split()
        rc, wall, rss = run(cmd, rdir, os.path.join(rdir, "run.stdout"))
        rep_ = parse_iqtree_report(os.path.join(rdir, "run.iqtree"))
        lh, gr = parse_ag_stats(os.path.join(rdir, "run.stdout"))
        rw, rp = rmse_vs_truth(rep_["comps"], truth)
        row = [name, seed, tmode, thr, arm, rep, "ok" if rc == 0 else "exit%d" % rc, "%.4f" % rep_["logl"], "%.2f" % wall,
               "%.1f" % rss, "%.4f" % rep_["tree_length"], lh, gr, "%.5f" % rw, "%.5f" % rp, " ".join(cmd)]
        line = "\t".join(str(x) for x in row)
        with open(results_path, "a") as f:
            f.write(line + "\n")
        print("done: %s %s %s t%d %s r%d  lnL=%s  wall=%.1fs" % (name, seed, tmode, thr, arm, rep, row[7], wall), flush=True)

    with ThreadPoolExecutor(max_workers=jobs) as ex:
        list(ex.map(one, tasks))

    bench_report.write_summary(bench_report.read_results(results_path), os.path.join(out_dir, "summary.md"), manifest)
    print("wrote", os.path.join(out_dir, "summary.md"))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
