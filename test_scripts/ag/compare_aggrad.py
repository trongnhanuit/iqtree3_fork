#!/usr/bin/env python3
"""Compare two --ag-dump-gradient (.aggrad.tsv) files.

Checks the branch-length gradients (the per-edge table's dlogl_dt column),
the raw per-mixture rows (#dlogl_dQ and every #nat_* row -- where SIMD lane
reordering actually shows up, section 4/5 of the design doc), and the final
unconstrained-space gradient (#theta) -- together, everything the engine
computes and the optimiser consumes. Exits 0 and prints PASS/max_rel, or
exits 1 and prints the worst mismatch.
"""
import sys

RAW_TAGS = ("#dlogl_dQ", "#nat_dR", "#nat_dpi", "#nat_dw", "#nat_drate",
            "#nat_dprop", "#nat_dpinv_direct", "#nat_dpinv", "#nat_dalpha")


def parse(path):
    logl = None
    theta = {}
    dt = {}
    raw = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#logl"):
                logl = float(line.split("\t")[1])
            elif line.startswith("#theta"):
                for tok in line[len("#theta"):].split():
                    name, rest = tok.split("=", 1)
                    _, grad = rest.split(":")
                    theta[name] = float(grad)
            elif any(line.startswith(tag + "\t") for tag in RAW_TAGS):
                tag, rest = line.split("\t", 1)
                for i, v in enumerate(rest.split("\t")):
                    raw[f"{tag[1:]}[{i}]"] = float(v)
            elif line.startswith("id\t"):
                continue
            elif not line.startswith("#"):
                parts = line.split("\t")
                if len(parts) >= 3:
                    dt[parts[0]] = float(parts[2])
    return logl, theta, dt, raw


def worst_rel(a, b):
    worst_rel, worst_key, worst_abs = 0.0, None, 0.0
    for k in a:
        if k not in b:
            continue
        d = abs(a[k] - b[k])
        rel = d / max(abs(a[k]), abs(b[k]), 1e-12)
        if rel > worst_rel:
            worst_rel, worst_key, worst_abs = rel, k, d
    return worst_rel, worst_key, worst_abs


def main():
    if len(sys.argv) != 3:
        print("usage: compare_aggrad.py <a.aggrad.tsv> <b.aggrad.tsv>", file=sys.stderr)
        return 2
    logl_a, theta_a, dt_a, raw_a = parse(sys.argv[1])
    logl_b, theta_b, dt_b, raw_b = parse(sys.argv[2])
    if logl_a != logl_b:
        print(f"FAIL logl differs: {logl_a} vs {logl_b}")
        return 1
    candidates = [
        (*worst_rel(theta_a, theta_b), "theta"),
        (*worst_rel(dt_a, dt_b), "dt"),
        (*worst_rel(raw_a, raw_b), "raw"),
    ]
    rel, key, abs_, kind = max(candidates, key=lambda c: c[0])
    tol = 1e-9
    if rel > tol:
        print(f"FAIL max_rel={rel:.3e} (> {tol:.0e}) worst={kind}:{key} abs={abs_:.3e}")
        return 1
    print(f"PASS max_rel={rel:.3e} worst={kind}:{key}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
