#!/usr/bin/env python3
"""
Report generator for the --analytical-gradients benchmark (bench.py).

Reads results.tsv (one row per run) and writes summary.md: for every
dataset / tree mode / thread count a table of the arms with the final
log-likelihood, its difference to the best arm, the median wall time, the
speed-up against the old default arm, the peak memory, the optimiser's own
evaluation counts, and (for simulated data) the RMSE of the recovered
mixture weights and profiles against the truth. Plots are written when
matplotlib is importable.

Usage: bench_report.py <results.tsv> [summary.md]
"""
import sys
import os
import json
import math
from collections import defaultdict


def read_results(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            if not line.strip():
                continue
            vals = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, vals)))
    return rows


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def median(xs):
    xs = sorted(x for x in xs if not math.isnan(x))
    if not xs:
        return float("nan")
    n = len(xs)
    return xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])


def fmt(x, digits=2):
    return "-" if math.isnan(x) else ("%.*f" % (digits, x))


def summarise(rows):
    """group by (dataset, tree_mode, threads, arm); repeats collapse to medians"""
    groups = defaultdict(list)
    for r in rows:
        groups[(r["dataset"], r["tree_mode"], r["threads"], r["arm"])].append(r)
    out = {}
    for key, rs in groups.items():
        out[key] = {
            "logl": median([fnum(r["logl"]) for r in rs]),
            "wall": median([fnum(r["wall_s"]) for r in rs]),
            "rss": median([fnum(r["peak_rss_mb"]) for r in rs]),
            "treelen": median([fnum(r["tree_length"]) for r in rs]),
            "lh_evals": median([fnum(r["lh_evals"]) for r in rs]),
            "grad_evals": median([fnum(r["grad_evals"]) for r in rs]),
            "rmse_w": median([fnum(r["rmse_weights"]) for r in rs]),
            "rmse_pi": median([fnum(r["rmse_profiles"]) for r in rs]),
            "n": len(rs),
            "failed": sum(1 for r in rs if r.get("status", "ok") != "ok"),
        }
    return out


def write_summary(rows, path, manifest=None):
    summ = summarise(rows)
    blocks = defaultdict(list)
    for (ds, mode, thr, arm), v in summ.items():
        blocks[(ds, mode, thr)].append((arm, v))
    lines = ["# Analytical-gradients benchmark", ""]
    if manifest:
        lines.append("Binary: `%s`  " % manifest.get("binary", "?"))
        lines.append("Git: `%s`  " % manifest.get("git", "?"))
        lines.append("Mode: %s, repeats: %s, host: %s" % (manifest.get("mode"), manifest.get("repeats"), manifest.get("host")))
        lines.append("")
    lines.append("Columns: lnL = median final log-likelihood over repeats; dlnL = difference to the best optimiser arm "
                 "(0 = best); wall = median seconds; speedup = old-default wall / this wall; RSS = peak memory (MB); "
                 "evals = the optimiser's likelihood/gradient evaluation counts (new arms only); RMSE = recovered "
                 "mixture weights / profiles against the simulation truth (simulated data only). The `truth` row, "
                 "where present, is the simulation model evaluated on the true tree with fixed branch lengths "
                 "(no optimisation): a reference for lnL, not an arm.")
    lines.append("")
    for (ds, mode, thr) in sorted(blocks):
        arms = blocks[(ds, mode, thr)]
        best = max(v["logl"] for a, v in arms if a != "truth" and not math.isnan(v["logl"]))
        old = [v for a, v in arms if a == "old-default"]
        old_wall = old[0]["wall"] if old else float("nan")
        lines.append("## %s, tree: %s, threads: %s" % (ds, mode, thr))
        lines.append("")
        lines.append("| arm | lnL | dlnL | wall (s) | speedup | RSS (MB) | evals lh/grad | RMSE w | RMSE pi | n |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|")
        for arm, v in sorted(arms, key=lambda av: -av[1]["logl"] if not math.isnan(av[1]["logl"]) else 1e9):
            speed = old_wall / v["wall"] if arm != "truth" and v["wall"] and not math.isnan(old_wall) else float("nan")
            evals = "-" if math.isnan(v["lh_evals"]) else "%d/%d" % (v["lh_evals"], v["grad_evals"])
            note = " (%d failed)" % v["failed"] if v["failed"] else ""
            lines.append("| %s%s | %s | %s | %s | %s | %s | %s | %s | %s | %d |" % (
                arm, note, fmt(v["logl"], 3), fmt(v["logl"] - best, 3), fmt(v["wall"], 1), fmt(speed, 2),
                fmt(v["rss"], 0), evals, fmt(v["rmse_w"], 4), fmt(v["rmse_pi"], 4), v["n"]))
        lines.append("")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    try:
        plot(summ, os.path.dirname(os.path.abspath(path)))
    except Exception as e:  # matplotlib missing or headless problems: the tables are the deliverable
        lines.append("(plots skipped: %s)" % e)
    return summ


def plot(summ, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    blocks = defaultdict(list)
    for (ds, mode, thr, arm), v in summ.items():
        blocks[(ds, mode, thr)].append((arm, v))
    for (ds, mode, thr), arms in blocks.items():
        arms = sorted(arms)
        names = [a for a, _ in arms]
        best = max(v["logl"] for a, v in arms if a != "truth" and not math.isnan(v["logl"]))
        fig, ax = plt.subplots(1, 2, figsize=(10, 3.5))
        ax[0].bar(names, [v["logl"] - best for _, v in arms])
        ax[0].set_ylabel("lnL - best")
        ax[1].bar(names, [v["wall"] for _, v in arms])
        ax[1].set_ylabel("wall (s)")
        for a in ax:
            a.tick_params(axis="x", rotation=45)
        fig.suptitle("%s, tree %s, %s threads" % (ds, mode, thr))
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "plot_%s_%s_t%s.png" % (ds, mode, thr)))
        plt.close(fig)


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    results = argv[1]
    out = argv[2] if len(argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(results)), "summary.md")
    manifest = None
    mpath = os.path.join(os.path.dirname(os.path.abspath(results)), "manifest.json")
    if os.path.exists(mpath):
        with open(mpath) as f:
            manifest = json.load(f)
    write_summary(read_results(results), out, manifest)
    print("wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
