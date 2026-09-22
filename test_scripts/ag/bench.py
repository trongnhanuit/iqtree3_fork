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
  bench.py <iqtree_binary> <out_dir> [--quick | --full | --thorough] [--big] [-j N]
           [--threads 1,8] [--repeats N] [--seeds 101,102] [--only name,name]

--thorough: LG+F10, GTR20+F12 and GTR20+C60 (with and without -mwopt) under
+R8, +I+G4 and +I+R10 variants, 32 taxa x 3000 sites (C60: 24 x 2000), three replicates, on
the true tree, every analytic setting (warm/cold start, cascade, multi-start,
without EM) against the default optimiser and its EM variant; -j defaults
to 70% of the cores.

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


def cmix(base, k, n):
    """MIX of the first n profiles of C<k> under exchangeabilities `base`, weights renormalised"""
    prof, w = cseries(k)
    names = ["C%dpi%d" % (k, i + 1) for i in range(n)]
    tot = sum(w[x] for x in names)
    return "MIX{%s}" % ",".join("%s+F{%s}:1:%.6f" % (base, "/".join("%g" % v for v in prof[x]), w[x] / tot) for x in names)


IG = "+I{0.15}+G4{0.7}"
IR10 = "+I{0.1}+R10{0.1,0.1,0.1,0.25,0.1,0.45,0.1,0.7,0.1,0.95,0.1,1.2,0.1,1.5,0.1,1.9,0.1,2.5,0.1,3.6}"

# name -> dict(sim=(alisim model, ntaxa, length) | real=path | alias_of=name, model=iqtree model args, tiers)
# tier "thorough": the user's three model families with several rate-heterogeneity variants (32 taxa x 3000 sites)
DATASETS = {
    "th_lg_f10_r8":       dict(sim=("LG+C10" + R8, 32, 3000), model="LG+F10+R8", tier="thorough"),
    "th_lg_f10_ig":       dict(sim=("LG+C10" + IG, 32, 3000), model="LG+F10+I+G4", tier="thorough"),
    "th_lg_f10_ir10":     dict(sim=("LG+C10" + IR10, 32, 3000), model="LG+F10+I+R10", tier="thorough"),
    "th_gtr20_f12_r8":    dict(sim=(cmix("WAG", 20, 12) + R8, 32, 3000), model="GTR20+F12+R8 --gtr20-model LG", tier="thorough"),
    "th_gtr20_f12_ig":    dict(sim=(cmix("WAG", 20, 12) + IG, 32, 3000), model="GTR20+F12+I+G4 --gtr20-model LG", tier="thorough"),
    "th_gtr20_c60_r8":    dict(sim=("WAG+C60" + R8, 24, 2000), model="GTR20+C60+R8 --gtr20-model LG", tier="thorough"),
    "th_gtr20_c60_r8_mwopt": dict(alias_of="th_gtr20_c60_r8", model="GTR20+C60+R8 --gtr20-model LG -mwopt", tier="thorough"),
    "th_gtr20_c60_ig":    dict(sim=("WAG+C60" + IG, 24, 2000), model="GTR20+C60+I+G4 --gtr20-model LG", tier="thorough"),
    "th_gtr20_c60_ig_mwopt": dict(alias_of="th_gtr20_c60_ig", model="GTR20+C60+I+G4 --gtr20-model LG -mwopt", tier="thorough"),
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
    # thorough tier: every analytic setting against both default-optimiser variants
    "new-cascade":    ("--analytical-gradients --ag-stats --ag-cascade on", "thorough"),
    "new-noem":       ("--analytical-gradients --ag-stats --ag-em-axes none", "thorough"),   # "none": no W/R/F letters -> no EM step
    "new-cold-multi": ("--analytical-gradients --ag-stats --ag-start cold --ag-multistart -1", "thorough"),
}
THOROUGH_ARMS = ["old-default", "old-em", "new-warm", "new-cascade", "new-cold", "new-multi", "new-noem", "new-cold-multi"]


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
        try:
            rc = subprocess.call(full, cwd=cwd, stdout=lf, stderr=subprocess.STDOUT, timeout=RUN_TIMEOUT)
        except subprocess.TimeoutExpired:
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
        truth = truth_from_alisim(d["sim"][0]) if "sim" in d else None
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
    mode, big, jobs, threads, repeats, seeds, only = "quick", False, 4, [1], 1, [101], None
    i = 3
    while i < len(argv):
        a = argv[i]
        if a == "--quick": mode = "quick"
        elif a == "--full": mode, threads, repeats = "full", [1, 8], 3
        elif a == "--thorough": mode, seeds, jobs = "thorough", [101, 102, 103], max(1, int((os.cpu_count() or 4) * 0.7))
        elif a == "--big": big = True
        elif a == "-j": i += 1; jobs = int(argv[i])
        elif a == "--threads": i += 1; threads = [int(x) for x in argv[i].split(",")]
        elif a == "--repeats": i += 1; repeats = int(argv[i])
        elif a == "--seeds": i += 1; seeds = [int(x) for x in argv[i].split(",")]
        elif a == "--only": i += 1; only = set(argv[i].split(","))
        else:
            print("unknown argument", a); return 2
        i += 1
    os.makedirs(out_dir, exist_ok=True)
    tiers = {"quick"} if mode == "quick" else {"quick", "full"} if mode == "full" else {"thorough"}
    if big:
        tiers.add("big")
    datasets = {n: d for n, d in DATASETS.items() if d["tier"] in tiers and (only is None or n in only)}
    if mode == "thorough":
        arms = {n: ARMS[n][0] for n in THOROUGH_ARMS}
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
        model, ntaxa, length = d["sim"]
        for seed in seeds:
            prefix = os.path.join(sim_dir, "%s_s%d" % (name, seed))
            if not os.path.exists(prefix + ".phy"):
                rc, _, _ = run([binary, "--alisim", prefix, "-t", "RANDOM{yh,%d}" % ntaxa, "-rlen", "0.01", "0.1", "0.5",
                                "-m", model, "--length", str(length), "-seed", str(seed), "-redo", "-quiet"], sim_dir, prefix + ".alisim.log")
                if rc != 0:
                    print("simulation failed:", name, seed); return 1
            inputs[(name, seed)] = (prefix + ".phy", prefix + ".treefile", truth_from_alisim(model))
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
        for tmode in modes:
            for thr in threads:
                for arm, extra in arms.items():
                    model = d["model"]
                    if arm == "old-em" and not re.search(r"\+(F\d+|C\d+)", model.split()[0]):
                        continue   # -optalg_qmix EM only applies to mixtures
                    if extra == "__C10WARM__":
                        cm = c10warm_model(model.split()[0])
                        if cm is None:
                            continue
                        extra = "-mfopt"
                        model = cm + model[len(model.split()[0]):]
                    for rep in range(repeats):
                        tasks.append((name, seed, tmode, thr, arm, rep, aln, tree, truth, model, extra))

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
