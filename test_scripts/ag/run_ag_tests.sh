#!/bin/bash
# Test driver for the --analytical-gradients work.
#
# Usage:
#   test_scripts/ag/run_ag_tests.sh <iqtree_binary> [out_dir] [--suite gradcheck|oracle|threads|all] [-j N]
#
# Suites:
#   gradcheck  --ag-gradient-check-only on a fixed set of models; PASS = exit 0
#              (every analytic entry within tolerance of its reference, and the
#              per-edge likelihood self-check passes)
#   oracle     AliSim-simulated small datasets; IQ-TREE's raw gradient dump is
#              compared with the independent NumPy oracle (oracle.py), including
#              the oracle self-test that must fail on a perturbed gradient
#   threads    the same model at -nt 1 and -nt 4 must give equal analytic gradients
#              (relative 1e-9); also documents that lnL agrees
#
# Cases run in parallel, at most N processes at a time (default: half the cores).
# Exit code 0 = all cases passed, 1 = a failure, 2 = usage.

set -u

BIN="${1:-}"
OUT_DIR="${2:-ag_test_out}"
SUITE="all"
MAXJOBS=0
argv=("$@")
for ((k=0; k<${#argv[@]}; k++)); do
    case "${argv[$k]}" in
        --suite) SUITE="${argv[$((k+1))]:-all}" ;;
        -j) MAXJOBS="${argv[$((k+1))]:-0}" ;;
    esac
done
if [ -z "$BIN" ] || [ ! -x "$BIN" ]; then
    echo "usage: $0 <iqtree_binary> [out_dir] [--suite gradcheck|oracle|threads|all] [-j N]" >&2
    exit 2
fi
if [ "$MAXJOBS" -le 0 ]; then
    NCPU=$( (nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2) )
    MAXJOBS=$(( NCPU / 2 )); [ "$MAXJOBS" -lt 1 ] && MAXJOBS=1
fi
BIN=$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
EX="$ROOT/example"
WD="$ROOT/test_scripts/test_data"
SEED=1
mkdir -p "$OUT_DIR"; OUT_DIR=$(cd "$OUT_DIR" && pwd)
REP="$OUT_DIR/report"; mkdir -p "$REP"; rm -f "$REP"/*

# ---- gradcheck cases: "<id>|<iqtree arguments>" (check runs at the initial point) ----
GRAD=(
  "g_dna_gtr_g4|-s $EX/example.phy -m GTR+F+G4 -te $HERE/data/example_gtr_g.nwk"
  "g_dna_gtr_fo_r4|-s $EX/example.phy -m GTR+FO+R4 -te $HERE/data/example_gtr_g.nwk"
  "g_dna_gtr_i_g4|-s $EX/example.phy -m GTR+F+I+G4 -te $HERE/data/example_gtr_g.nwk"
  "g_dna_hky_f|-s $EX/example.phy -m HKY{2.0}+F -te $HERE/data/example_gtr_g.nwk"
  "g_dna_jc|-s $EX/example.phy -m JC -te $HERE/data/example_gtr_g.nwk"
  "g_dna_f81|-s $EX/example.phy -m F81+F -te $HERE/data/example_gtr_g.nwk"
  "g_aa_lg_f2_g4|-s $EX/aa_example.phy -m LG+F2+G4 -te $HERE/data/aa_example_lg.nwk"
  "g_aa_lg_c10_r4|-s $WD/turtle_aa.fasta -m LG+C10+R4"
  "g_aa_mix_unlinked|-s $EX/aa_example.phy -m MIX{LG+FO,WAG+FO}+G4 -te $HERE/data/aa_example_lg.nwk"
  "g_dna_safe_r4|-s $EX/example.phy -m GTR+F+R4 -safe -te $HERE/data/example_gtr_g.nwk"
  "g_aa_safe_c10|-s $WD/turtle_aa.fasta -m LG+C10+G4 -safe"
)

# ---- oracle cases: "<id>|<alisim model>|<newick>|<iqtree model>" ----
ORACLE=(
  "o_dna_gtr_g4|GTR{1.5,3,0.8,1.2,2.5}+F{0.3,0.2,0.2,0.3}+G4{0.7}|(A:0.1,B:0.2,(C:0.15,D:0.05):0.1);|GTR+F+G4"
  "o_dna_hky_i_g4|HKY{2.5}+F{0.35,0.15,0.2,0.3}+I{0.2}+G4{0.5}|(A:0.05,B:0.3,(C:0.1,D:0.25):0.15);|HKY+F+I+G4"
  "o_aa_lg_g4|LG+G4{0.8}|(A:0.2,(B:0.1,C:0.3):0.1,D:0.15);|LG+G4"
  "o_dna_5tax_r3|GTR{1,2,1,1,3}+F{0.25,0.25,0.25,0.25}+R3{0.3,0.2,0.4,0.8,0.3,2.0}|((A:0.1,B:0.2):0.05,(C:0.15,D:0.05):0.1,E:0.3);|GTR+F+R3"
)

run_grad() {   # id args -> report/<id>.txt, marker FAIL
    local id="$1" args="$2" dir="$OUT_DIR/$1"
    mkdir -p "$dir"
    # shellcheck disable=SC2086
    ( cd "$dir" && "$BIN" $args -nt 1 -seed $SEED --prefix "$id" -redo --analytical-gradients --ag-gradient-check-only > "$id.stdout" 2>&1 )
    local rc=$?
    { echo "== $id (exit $rc)"; grep -E "^AG:|GRADCHECK" "$dir/$id.stdout"; } > "$REP/$id.txt"
    [ "$rc" = "0" ] || { touch "$REP/$id.FAIL"; grep -E "ERROR|FAIL" "$dir/$id.stdout" | head -3 >> "$REP/$id.txt"; }
}

run_oracle() {   # id alisim_model newick iqtree_model
    local id="$1" amodel="$2" nwk="$3" imodel="$4" dir="$OUT_DIR/$1" rc=0
    mkdir -p "$dir"
    {
        cd "$dir" || exit 1
        echo "$nwk" > tree.nwk
        "$BIN" --alisim sim -t tree.nwk -m "$amodel" --length 400 -seed 7 -redo -quiet > alisim.stdout 2>&1 || { echo "  alisim failed"; exit 1; }
        "$BIN" -s sim.phy -m "$imodel" -te tree.nwk -nt 1 -seed $SEED --prefix "$id" -redo --analytical-gradients --ag-gradient-check-only --ag-dump-gradient > "$id.stdout" 2>&1
        rc=$?
        grep -E "^AG:|GRADCHECK" "$id.stdout"
        [ "$rc" = "0" ] || { echo "  iqtree gradient check failed (exit $rc)"; exit 1; }
        python3 "$HERE/oracle.py" sim.phy "$id.aggrad.tsv" --tol 1e-6 > oracle.txt 2>&1 || { echo "  oracle FAILED:"; tail -4 oracle.txt; exit 1; }
        tail -1 oracle.txt
        python3 "$HERE/oracle.py" sim.phy "$id.aggrad.tsv" --selftest > selftest.txt 2>&1 || { echo "  oracle self-test did not detect a perturbed gradient"; exit 1; }
        echo "  oracle self-test: perturbed gradient detected"
    } > "$REP/$id.txt" 2>&1
    rc=$?
    sed -i.bak "1i\\
== $id (exit $rc)" "$REP/$id.txt" 2>/dev/null || { (echo "== $id (exit $rc)"; cat "$REP/$id.txt") > "$REP/$id.tmp" && mv "$REP/$id.tmp" "$REP/$id.txt"; }
    rm -f "$REP/$id.txt.bak"
    [ "$rc" = "0" ] || touch "$REP/$id.FAIL"
}

run_threads() {   # compare analytic columns at -nt 1 vs -nt 4
    local id="t_threads" dir="$OUT_DIR/$id" rc=0
    mkdir -p "$dir"
    {
        cd "$dir" || exit 1
        for nt in 1 4; do
            "$BIN" -s "$EX/aa_example.phy" -m LG+F2+G4 -te "$HERE/data/aa_example_lg.nwk" -nt $nt -seed $SEED --prefix "nt$nt" -redo --analytical-gradients --ag-gradient-check-only > "nt$nt.stdout" 2>&1 || { echo "  run at -nt $nt failed"; exit 1; }
        done
        python3 - <<'PY' || exit 1
import math
def load(f):
    rows = [l.rstrip("\n").split("\t") for l in open(f)][1:]
    return {r[3]: (float(r[4]), float(r[5])) for r in rows}
a, b = load("nt1.gradcheck.tsv"), load("nt4.gradcheck.tsv")
gmax = max(abs(v[1]) for v in a.values())
worst = 0.0
for k in a:
    ga, gb = a[k][1], b[k][1]
    worst = max(worst, abs(ga-gb)/max(abs(ga), abs(gb), 1e-6*gmax))
print("  -nt 1 vs -nt 4: %d branches, max rel diff %.2e -> %s" % (len(a), worst, "PASS" if worst <= 1e-9 else "FAIL"))
raise SystemExit(0 if worst <= 1e-9 else 1)
PY
    } > "$REP/$id.txt" 2>&1
    rc=$?
    (echo "== $id (exit $rc)"; cat "$REP/$id.txt") > "$REP/$id.tmp" && mv "$REP/$id.tmp" "$REP/$id.txt"
    [ "$rc" = "0" ] || touch "$REP/$id.FAIL"
}

throttle() { while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$MAXJOBS" ]; do sleep 1; done; }

ORDER=()
if [ "$SUITE" = "gradcheck" ] || [ "$SUITE" = "all" ]; then
    for e in "${GRAD[@]}"; do id="${e%%|*}"; ORDER+=("$id"); throttle; run_grad "$id" "${e#*|}" & done
fi
if [ "$SUITE" = "oracle" ] || [ "$SUITE" = "all" ]; then
    for e in "${ORACLE[@]}"; do
        IFS='|' read -r id amodel nwk imodel <<< "$e"
        ORDER+=("$id"); throttle; run_oracle "$id" "$amodel" "$nwk" "$imodel" &
    done
fi
if [ "$SUITE" = "threads" ] || [ "$SUITE" = "all" ]; then
    ORDER+=("t_threads"); throttle; run_threads &
fi
wait

fail=0
for id in "${ORDER[@]}"; do
    cat "$REP/$id.txt"
    [ -f "$REP/$id.FAIL" ] && fail=1
done
echo "cases: ${#ORDER[@]}, failures: $(ls "$REP"/*.FAIL 2>/dev/null | wc -l | tr -d ' ')"
if [ "$fail" = "0" ]; then echo "AG TESTS PASSED"; else echo "AG TESTS FAILED"; fi
exit $fail
