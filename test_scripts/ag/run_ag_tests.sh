#!/bin/bash
# Test driver for the --analytical-gradients work.
#
# Usage:
#   test_scripts/ag/run_ag_tests.sh <iqtree_binary> [out_dir] [--suite gradcheck|oracle|threads|quality|robust|all] [-j N] [--full]
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
#   quality    the live optimiser: the same command with the flag off and on; the
#              final log-likelihood with the flag must be >= the default's minus a
#              tolerance (times are reported); -Q must use the new path per
#              partition and -p must fall back with the warning; a simulated
#              two-profile mixture must be recovered within absolute tolerances,
#              and warm, cold and multi-start must reach the same optimum on it
#   robust     integration behaviour: abort/resume from the checkpoint, checkpoint
#              cross-compatibility with the default path, ModelFinder, PMSF's
#              site-specific pass, +I+G restarts, fault injection, concurrent
#              partitions at 4 threads
#
# Cases run in parallel, at most N processes at a time (default: half the cores).
# Exit code 0 = all cases passed, 1 = a failure, 2 = usage.

set -u

BIN="${1:-}"
OUT_DIR="${2:-ag_test_out}"
SUITE="all"
MAXJOBS=0
FULL=0
argv=("$@")
for ((k=0; k<${#argv[@]}; k++)); do
    case "${argv[$k]}" in
        --suite) SUITE="${argv[$((k+1))]:-all}" ;;
        -j) MAXJOBS="${argv[$((k+1))]:-0}" ;;
        --full) FULL=1 ;;
    esac
done
if [ -z "$BIN" ] || [ ! -x "$BIN" ]; then
    echo "usage: $0 <iqtree_binary> [out_dir] [--suite gradcheck|oracle|threads|quality|robust|all] [-j N] [--full]" >&2
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
  "g_dna_gtr_fo_r4_i|-s $EX/example.phy -m GTR+FO+R4+I -te $HERE/data/example_gtr_g.nwk"
  "g_dna_link_mix|-s $EX/example.phy -m MIX{GTR+FO,GTR+FO}+G4 --link-exchange-rates -te $HERE/data/example_gtr_g.nwk"
  "g_aa_gtr20_link|-s $EX/aa_example.phy -m GTR20+F2+R2 --gtr20-model LG -te $HERE/data/aa_example_lg.nwk"
  "g_aa_lg_f10|-s $EX/aa_example.phy -m LG+F10 -te $HERE/data/aa_example_lg.nwk"
  "g_dna_legacy_eigen|-s $EX/example.phy -m GTR+FO+G4 --eigen -te $HERE/data/example_gtr_g.nwk"
  "g_selftest|-s $EX/example.phy -m GTR+FO+R3+I -te $HERE/data/example_gtr_g.nwk --ag-selftest"
)

# ---- oracle cases: "<id>|<alisim model>|<newick>|<iqtree model>" ----
ORACLE=(
  "o_dna_gtr_g4|GTR{1.5,3,0.8,1.2,2.5}+F{0.3,0.2,0.2,0.3}+G4{0.7}|(A:0.1,B:0.2,(C:0.15,D:0.05):0.1);|GTR+F+G4"
  "o_dna_hky_i_g4|HKY{2.5}+F{0.35,0.15,0.2,0.3}+I{0.2}+G4{0.5}|(A:0.05,B:0.3,(C:0.1,D:0.25):0.15);|HKY+F+I+G4"
  "o_aa_lg_g4|LG+G4{0.8}|(A:0.2,(B:0.1,C:0.3):0.1,D:0.15);|LG+G4"
  "o_dna_5tax_r3|GTR{1,2,1,1,3}+F{0.25,0.25,0.25,0.25}+R3{0.3,0.2,0.4,0.8,0.3,2.0}|((A:0.1,B:0.2):0.05,(C:0.15,D:0.05):0.1,E:0.3);|GTR+F+R3"
)

# ---- quality cases: "<id>|<tolerance>|<iqtree arguments>" (run with the flag off and on) ----
QUALITY=(
  "q_dna_gtrfo_g4_search|0.1|-s $EX/example.phy -m GTR+FO+G4"
  "q_dna_gtr_i_g4_search|0.1|-s $EX/example.phy -m GTR+F+I+G4"
  "q_aa_lg_f2_g4_te|0.1|-s $EX/aa_example.phy -m LG+F2+G4 -te $HERE/data/aa_example_lg.nwk"
  "q_dna_mix_link_te|0.1|-s $EX/example.phy -m MIX{GTR+FO,GTR+FO}+G4 --link-exchange-rates -te $HERE/data/example_gtr_g.nwk"
)
QUALITY_FULL=(
  "q_aa_lg_f4_r4_search|0.1|-s $WD/turtle_aa.fasta -m LG+F4+R4"
  "q_aa_lg_c10_r4_search|0.1|-s $WD/turtle_aa.fasta -m LG+C10+R4"
)
[ "$FULL" = "1" ] && QUALITY+=("${QUALITY_FULL[@]}")

# Sanitizer mode (AG_SANITIZER=1, binary built with -fsanitize=address,undefined and
# UBSAN_OPTIONS=halt_on_error=0): IQ-TREE's existing code has UBSan findings of its
# own (e.g. an uninitialised bool read in ModelMarkov's constructor), so a case fails
# only on (a) any AddressSanitizer report, (b) an undefined-behaviour report whose
# location is in the new files, or (c) a failed gradient check. Other reports are
# listed as pre-existing and ignored.
SAN_NEW_FILES='phylogradient|gradientoptimizer|modelparammap'

run_grad() {   # id args -> report/<id>.txt, marker FAIL
    local id="$1" args="$2" dir="$OUT_DIR/$1"
    mkdir -p "$dir"
    # shellcheck disable=SC2086
    ( cd "$dir" && "$BIN" $args -nt 1 -seed $SEED --prefix "$id" -redo --analytical-gradients --ag-gradient-check-only > "$id.stdout" 2>&1 )
    local rc=$?
    { echo "== $id (exit $rc)"; grep -E "^AG:|GRADCHECK|SELFTEST" "$dir/$id.stdout"; } > "$REP/$id.txt"
    if [ "${AG_SANITIZER:-0}" = "1" ]; then
        local asan ubsan_new ubsan_old
        asan=$(grep -c "ERROR: AddressSanitizer" "$dir/$id.stdout")
        ubsan_new=$(grep -E "runtime error" "$dir/$id.stdout" | grep -c -E "$SAN_NEW_FILES")
        ubsan_old=$(grep -E "runtime error" "$dir/$id.stdout" | grep -v -c -E "$SAN_NEW_FILES")
        echo "  sanitizer: asan_reports=$asan ubsan_in_new_files=$ubsan_new ubsan_pre_existing=$ubsan_old" >> "$REP/$id.txt"
        if [ "$asan" != "0" ] || [ "$ubsan_new" != "0" ]; then
            touch "$REP/$id.FAIL"
            grep -E "ERROR: AddressSanitizer|runtime error" "$dir/$id.stdout" | grep -E "AddressSanitizer|$SAN_NEW_FILES" | head -3 >> "$REP/$id.txt"
        fi
        # the gradient check itself must still pass (GRADCHECK line with n_fail=0)
        grep -q "GRADCHECK.* n_fail=0 .*edge_lnl_check=PASS.*identities=PASS" "$dir/$id.stdout" || touch "$REP/$id.FAIL"
        grep -q "SELFTEST.*FAIL" "$dir/$id.stdout" && touch "$REP/$id.FAIL"
        return
    fi
    [ "$rc" = "0" ] || { touch "$REP/$id.FAIL"; grep -E "ERROR|FAIL" "$dir/$id.stdout" | head -3 >> "$REP/$id.txt"; }
}

run_oracle() {   # id alisim_model newick iqtree_model
    local id="$1" amodel="$2" nwk="$3" imodel="$4" dir="$OUT_DIR/$1" rc=0
    mkdir -p "$dir"
    (
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
    ) > "$REP/$id.txt" 2>&1
    rc=$?
    sed -i.bak "1i\\
== $id (exit $rc)" "$REP/$id.txt" 2>/dev/null || { (echo "== $id (exit $rc)"; cat "$REP/$id.txt") > "$REP/$id.tmp" && mv "$REP/$id.tmp" "$REP/$id.txt"; }
    rm -f "$REP/$id.txt.bak"
    [ "$rc" = "0" ] || touch "$REP/$id.FAIL"
}

run_threads() {   # compare analytic columns at -nt 1 vs -nt 4
    local id="t_threads"
    local dir="$OUT_DIR/$id" rc=0
    mkdir -p "$dir"
    (
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
print("  -nt 1 vs -nt 4: %d gradient entries, max rel diff %.2e -> %s" % (len(a), worst, "PASS" if worst <= 1e-9 else "FAIL"))
raise SystemExit(0 if worst <= 1e-9 else 1)
PY
    ) > "$REP/$id.txt" 2>&1
    rc=$?
    (echo "== $id (exit $rc)"; cat "$REP/$id.txt") > "$REP/$id.tmp" && mv "$REP/$id.tmp" "$REP/$id.txt"
    [ "$rc" = "0" ] || touch "$REP/$id.FAIL"
}

run_quality() {   # id tol args: flag off vs flag on, same binary
    local id="$1" tol="$2" args="$3" dir="$OUT_DIR/$1" rc=0
    mkdir -p "$dir"
    (
        cd "$dir" || exit 1
        local t0 t1 t2 old new
        t0=$(date +%s)
        if [ "${AG_SANITIZER:-0}" = "1" ]; then
            # sanitizer mode exercises the new code only: the flag-off run is the
            # default path, which the sanitizer job does not need and which is slow
            # under ASan; the flag-on run must finish without a report in new files
            # shellcheck disable=SC2086
            "$BIN" $args -nt 1 -seed $SEED --prefix on -redo --analytical-gradients --ag-stats > on.stdout 2>&1 || { echo "  flag-on run failed"; exit 1; }
            local asan ubsan_new
            asan=$(grep -c "ERROR: AddressSanitizer" on.stdout)
            ubsan_new=$(grep -E "runtime error" on.stdout | grep -c -E "$SAN_NEW_FILES")
            echo "  sanitizer: asan_reports=$asan ubsan_in_new_files=$ubsan_new"
            [ "$asan" = "0" ] && [ "$ubsan_new" = "0" ] || exit 1
            grep -c "AG stats" on.stdout | sed 's/^/  AG stats lines: /'
            exit 0
        fi
        # shellcheck disable=SC2086
        "$BIN" $args -nt 1 -seed $SEED --prefix off -redo > off.stdout 2>&1 || { echo "  flag-off run failed"; exit 1; }
        t1=$(date +%s)
        # shellcheck disable=SC2086
        "$BIN" $args -nt 1 -seed $SEED --prefix on -redo --analytical-gradients --ag-stats > on.stdout 2>&1 || { echo "  flag-on run failed"; exit 1; }
        t2=$(date +%s)
        old=$(grep -m1 "^Log-likelihood of the tree" off.iqtree | grep -Eo '[-]?[0-9]+\.[0-9]+' | head -1)
        new=$(grep -m1 "^Log-likelihood of the tree" on.iqtree | grep -Eo '[-]?[0-9]+\.[0-9]+' | head -1)
        grep -c "AG stats" on.stdout | sed 's/^/  AG stats lines: /'
        echo "  logl flag-off=$old flag-on=$new  time flag-off=$((t1-t0))s flag-on=$((t2-t1))s"
        awk -v a="$old" -v b="$new" -v tol="$tol" 'BEGIN { if (b >= a - tol) { print "  quality: PASS"; exit 0 } else { print "  quality: FAIL (flag-on worse than flag-off by more than " tol ")"; exit 1 } }'
    ) > "$REP/$id.txt" 2>&1
    rc=$?
    (echo "== $id (exit $rc)"; cat "$REP/$id.txt") > "$REP/$id.tmp" && mv "$REP/$id.tmp" "$REP/$id.txt"
    [ "$rc" = "0" ] || touch "$REP/$id.FAIL"
}

run_recovery() {   # simulated two-profile mixture: weights and profiles must be recovered (absolute tolerance)
    local id="q_recovery"
    local dir="$OUT_DIR/$id" rc=0
    mkdir -p "$dir"
    (
        cd "$dir" || exit 1
        if [ "${AG_SANITIZER:-0}" = "1" ]; then
            # AliSim has a pre-existing AddressSanitizer report of its own
            # (new-delete-type-mismatch), so the simulation cannot run under the
            # sanitizer binary; the optimiser itself is covered by the other cases
            echo "  skipped in sanitizer mode (AliSim is not sanitizer-clean)"; exit 0
        fi
        P1="0.18/0.1/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.02/0.02/0.02/0.02/0.02/0.02"
        P2="0.02/0.02/0.02/0.02/0.02/0.02/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.1/0.18"
        echo "(((A:0.08,B:0.12):0.05,(C:0.1,D:0.07):0.06):0.04,((E:0.09,F:0.11):0.05,(G:0.06,H:0.13):0.07):0.03);" > t.nwk
        "$BIN" --alisim sim -t t.nwk -m "MIX{LG+F{$P1}:1:0.3,LG+F{$P2}:1:0.7}+G4{0.8}" --length 1500 -seed 11 -redo -quiet > alisim.stdout 2>&1 || { echo "  alisim failed"; exit 1; }
        "$BIN" -s sim.phy -m "MIX{LG+FO,LG+FO}+G4" -te t.nwk -nt 1 -seed $SEED --prefix on -redo --analytical-gradients --ag-force --ag-stats > on.stdout 2>&1 || { echo "  flag-on run failed"; exit 1; }
        grep "AG stats" on.stdout | head -1
        python3 - "$P1" "$P2" <<'PY' || exit 1
import sys, re, itertools, math
p1 = [float(x) for x in sys.argv[1].split("/")]; p2 = [float(x) for x in sys.argv[2].split("/")]
truth = [(0.3, p1), (0.7, p2)]
rows = []
for l in open("on.iqtree"):
    m = re.match(r"\s*\d+\s+\S+\s+([0-9.]+)\s+([0-9.]+)\s+\S+FO\{([^}]*)\}", l)
    if m:
        rows.append((float(m.group(2)), [float(x) for x in m.group(3).split(",")]))
if len(rows) != 2:
    print("  could not parse two components from on.iqtree"); sys.exit(1)
best = None
for perm in itertools.permutations(range(2)):
    dw = max(abs(rows[perm[i]][0] - truth[i][0]) for i in range(2))
    rmse = max(math.sqrt(sum((a - b) ** 2 for a, b in zip(rows[perm[i]][1], truth[i][1])) / 20) for i in range(2))
    if best is None or dw + rmse < best[0] + best[1]:
        best = (dw, rmse)
print("  recovery: max |weight error| = %.3f (tol 0.08), max profile RMSE = %.4f (tol 0.02)" % best)
sys.exit(0 if best[0] <= 0.08 and best[1] <= 0.02 else 1)
PY
        if [ "$FULL" = "1" ]; then
            "$BIN" -s sim.phy -m "MIX{LG+FO,LG+FO}+G4" -te t.nwk -nt 1 -seed $SEED --prefix off -redo > off.stdout 2>&1 || { echo "  flag-off run failed"; exit 1; }
            old=$(grep -m1 "^Log-likelihood of the tree" off.iqtree | grep -Eo '[-]?[0-9]+\.[0-9]+' | head -1)
            new=$(grep -m1 "^Log-likelihood of the tree" on.iqtree | grep -Eo '[-]?[0-9]+\.[0-9]+' | head -1)
            echo "  logl flag-off=$old flag-on=$new"
            awk -v a="$old" -v b="$new" 'BEGIN { exit (b >= a - 0.1) ? 0 : 1 }' || { echo "  flag-on worse than flag-off"; exit 1; }
        fi
        echo "  recovery: PASS"
    ) > "$REP/$id.txt" 2>&1
    rc=$?
    (echo "== $id (exit $rc)"; cat "$REP/$id.txt") > "$REP/$id.tmp" && mv "$REP/$id.tmp" "$REP/$id.txt"
    [ "$rc" = "0" ] || touch "$REP/$id.FAIL"
}

run_starts() {   # warm (default), cold and multi-start must reach the same optimum on the simulated mixture
    local id="q_starts"
    local dir="$OUT_DIR/$id" rc=0
    mkdir -p "$dir"
    (
        cd "$dir" || exit 1
        if [ "${AG_SANITIZER:-0}" = "1" ]; then echo "  skipped in sanitizer mode (AliSim is not sanitizer-clean)"; exit 0; fi
        P1="0.18/0.1/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.02/0.02/0.02/0.02/0.02/0.02"
        P2="0.02/0.02/0.02/0.02/0.02/0.02/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.05/0.1/0.18"
        echo "(((A:0.08,B:0.12):0.05,(C:0.1,D:0.07):0.06):0.04,((E:0.09,F:0.11):0.05,(G:0.06,H:0.13):0.07):0.03);" > t.nwk
        "$BIN" --alisim sim -t t.nwk -m "MIX{LG+F{$P1}:1:0.3,LG+F{$P2}:1:0.7}+G4{0.8}" --length 1500 -seed 11 -redo -quiet > alisim.stdout 2>&1 || { echo "  alisim failed"; exit 1; }
        for v in "warm|--ag-start warm --ag-multistart 0" "cold|--ag-start cold --ag-multistart 0" "multi|--ag-multistart 20"; do
            name="${v%%|*}"; opts="${v#*|}"
            # shellcheck disable=SC2086
            "$BIN" -s sim.phy -m "MIX{LG+FO,LG+FO}+G4" -te t.nwk -nt 1 -seed $SEED --prefix $name -redo --analytical-gradients --ag-force --ag-stats $opts > $name.stdout 2>&1 || { echo "  $name run failed"; exit 1; }
            grep -E "AG: (cold|multi)" $name.stdout | head -2
            echo "  $name: $(grep -m1 "^Log-likelihood of the tree" $name.iqtree | cut -d" " -f1-5) $(grep "AG stats" $name.stdout | head -1 | grep -o "likelihood_evaluations=[0-9]*")"
        done
        python3 - <<'PY' || exit 1
import re
v = {}
for n in ("warm", "cold", "multi"):
    v[n] = float(re.search(r"tree: (-?[0-9.]+)", open(n + ".iqtree").read()).group(1))
best = max(v.values())
worst = min(v.values())
print("  starts: spread %.3f (tol 0.5)" % (best - worst))
raise SystemExit(0 if best - worst <= 0.5 else 1)
PY
        echo "  starts: PASS"
    ) > "$REP/$id.txt" 2>&1
    rc=$?
    (echo "== $id (exit $rc)"; cat "$REP/$id.txt") > "$REP/$id.tmp" && mv "$REP/$id.tmp" "$REP/$id.txt"
    [ "$rc" = "0" ] || touch "$REP/$id.FAIL"
}

run_partitions() {   # -Q takes the new path per partition; -p falls back with the warning
    local id="q_partitions"
    local dir="$OUT_DIR/$id" rc=0
    mkdir -p "$dir"
    (
        cd "$dir" || exit 1
        "$BIN" -s "$WD/turtle_aa.fasta" -Q "$WD/turtle_aa.nex" -m LG+F+G4 -nt 1 -seed $SEED --prefix q -redo --analytical-gradients --ag-stats > q.stdout 2>&1 || { echo "  -Q run failed"; exit 1; }
        n=$(grep -c "AG stats" q.stdout)
        echo "  -Q: $n per-partition AG stats lines"
        [ "$n" -ge 2 ] || { echo "  -Q did not use the analytic path"; exit 1; }
        "$BIN" -s "$WD/turtle_aa.fasta" -p "$WD/turtle_aa.nex" -m LG+F+G4 -nt 1 -seed $SEED --prefix p -redo --analytical-gradients > p.stdout 2>&1 || { echo "  -p run failed"; exit 1; }
        grep -q "not applicable to edge-linked partition models" p.stdout || { echo "  -p did not print the fallback warning"; exit 1; }
        [ "$(grep -c "AG stats" p.stdout)" = "0" ] || { echo "  -p used the analytic path"; exit 1; }
        echo "  -p: fell back with the warning"
    ) > "$REP/$id.txt" 2>&1
    rc=$?
    (echo "== $id (exit $rc)"; cat "$REP/$id.txt") > "$REP/$id.tmp" && mv "$REP/$id.tmp" "$REP/$id.txt"
    [ "$rc" = "0" ] || touch "$REP/$id.FAIL"
}

# ---- robust suite: integration behaviour of the live optimiser ----
robust_case() {   # id: runs the named check inside its own directory, marker FAIL on non-zero exit
    local id="$1"
    local dir="$OUT_DIR/$id" rc=0
    mkdir -p "$dir"
    ( cd "$dir" && robust_$id ) > "$REP/$id.txt" 2>&1
    rc=$?
    (echo "== $id (exit $rc)"; cat "$REP/$id.txt") > "$REP/$id.tmp" && mv "$REP/$id.tmp" "$REP/$id.txt"
    [ "$rc" = "0" ] || touch "$REP/$id.FAIL"
}
logl_of() { grep -m1 "^Log-likelihood of the tree" "$1" | grep -Eo '[-]?[0-9]+\.[0-9]+' | head -1; }
within() { awk -v a="$1" -v b="$2" -v tol="$3" 'BEGIN { d = a - b; if (d < 0) d = -d; exit (d <= tol) ? 0 : 1 }'; }

robust_x_abort_resume() {   # X1: abort after the start-point phase, resume from the checkpoint, same optimum
    local ARGS="-s $EX/aa_example.phy -m LG+F2+G4 -te $HERE/data/aa_example_lg.nwk -nt 1 -seed $SEED --analytical-gradients --ag-force"
    # shellcheck disable=SC2086
    "$BIN" $ARGS --prefix full -redo > full.stdout 2>&1 || { echo "  uninterrupted run failed"; return 1; }
    # shellcheck disable=SC2086
    "$BIN" $ARGS --prefix part -redo --ag-abort-after init > part1.stdout 2>&1
    grep -q "checkpoint written, exiting" part1.stdout || { echo "  abort did not happen"; return 1; }
    [ -f part.ckp.gz ] || { echo "  no checkpoint written"; return 1; }
    # shellcheck disable=SC2086
    "$BIN" $ARGS --prefix part > part2.stdout 2>&1 || { echo "  resumed run failed"; return 1; }
    local a b; a=$(logl_of full.iqtree); b=$(logl_of part.iqtree)
    echo "  uninterrupted=$a resumed=$b"
    # the report prints four decimals, so 1e-3 is the finest honest tolerance here
    within "$a" "$b" 1e-3 || { echo "  resumed optimum differs by more than 1e-3"; return 1; }
    echo "  abort/resume: PASS"
}

robust_x_ckp_crosscompat() {   # X11: a checkpoint written by either path is read by the other
    local ARGS="-s $EX/example.phy -m GTR+FO+G4 -te $HERE/data/example_gtr_g.nwk -nt 1 -seed $SEED"
    # shellcheck disable=SC2086
    "$BIN" $ARGS --prefix a -redo --analytical-gradients > a1.stdout 2>&1 || return 1
    local a1; a1=$(logl_of a.iqtree)
    # shellcheck disable=SC2086
    "$BIN" $ARGS --prefix a --undo > a2.stdout 2>&1 || { echo "  default path could not continue from the flagged checkpoint"; return 1; }
    grep -q "CHECKPOINT: Model parameters restored" a2.stdout || { echo "  model not restored from the flagged checkpoint"; return 1; }
    within "$a1" "$(logl_of a.iqtree)" 1e-3 || { echo "  logl changed after reading the flagged checkpoint"; return 1; }
    # shellcheck disable=SC2086
    "$BIN" $ARGS --prefix b -redo > b1.stdout 2>&1 || return 1
    local b1; b1=$(logl_of b.iqtree)
    # shellcheck disable=SC2086
    "$BIN" $ARGS --prefix b --undo --analytical-gradients > b2.stdout 2>&1 || { echo "  flagged path could not continue from the default checkpoint"; return 1; }
    grep -q "CHECKPOINT: Model parameters restored" b2.stdout || { echo "  model not restored from the default checkpoint"; return 1; }
    within "$b1" "$(logl_of b.iqtree)" 1e-3 || { echo "  logl changed after reading the default checkpoint"; return 1; }
    echo "  flagged->default $a1, default->flagged $b1: PASS"
}

robust_x_modelfinder() {   # X7: ModelFinder with the flag exits 0 and uses the analytic path for its candidates
    "$BIN" -s "$EX/example.phy" -m MF -mset GTR -mrate G,I+G -nt 1 -seed $SEED --prefix mf -redo --analytical-gradients --ag-stats > mf.stdout 2>&1 || { echo "  ModelFinder run failed"; return 1; }
    local n; n=$(grep -c "AG stats" mf.stdout)
    echo "  ModelFinder: $n analytic optimisations, best model $(grep -m1 "Best-fit model" mf.iqtree | cut -c1-60)"
    [ "$n" -ge 1 ] || { echo "  analytic path not used"; return 1; }
    echo "  modelfinder: PASS"
}

robust_x_pmsf() {   # X8: PMSF guide-tree fit uses the flag, the site-specific second pass falls back with a NOTE
    "$BIN" -s "$EX/aa_example.phy" -m LG+C10+G4 -ft "$HERE/data/aa_example_lg.nwk" -nt 1 -seed $SEED --prefix pmsf -redo --analytical-gradients --ag-stats > pmsf.stdout 2>&1 || { echo "  PMSF run failed"; return 1; }
    grep -q "not applicable (site-specific model)" pmsf.stdout || { echo "  no NOTE for the site-specific pass"; return 1; }
    [ "$(grep -c "AG stats" pmsf.stdout)" -ge 1 ] || { echo "  guide-tree pass did not use the analytic path"; return 1; }
    echo "  pmsf: PASS (guide-tree pass analytic, site-specific pass fell back)"
}

robust_x_gammai_restart() {   # X12: +I+G restarts (--opt-gamma-inv) with the flag exit 0
    "$BIN" -s "$EX/example.phy" -m GTR+F+I+G4 --opt-gamma-inv -te "$HERE/data/example_gtr_g.nwk" -nt 1 -seed $SEED --prefix gi -redo --analytical-gradients --ag-stats > gi.stdout 2>&1 || { echo "  run failed"; return 1; }
    "$BIN" -s "$EX/example.phy" -m GTR+F+I+G4 --opt-gamma-inv -te "$HERE/data/example_gtr_g.nwk" -nt 1 -seed $SEED --prefix gi0 -redo > gi0.stdout 2>&1 || return 1
    echo "  logl flag-off=$(logl_of gi0.iqtree) flag-on=$(logl_of gi.iqtree)"
    awk -v a="$(logl_of gi0.iqtree)" -v b="$(logl_of gi.iqtree)" 'BEGIN { exit (b >= a - 0.1) ? 0 : 1 }' || { echo "  flag-on worse"; return 1; }
    echo "  gamma-invar restarts: PASS"
}

robust_x_fault() {   # X14: an exception inside the outside pass leaves the tree reusable (fallback for that step)
    AG_TEST_FAULT_EDGE=3 "$BIN" -s "$EX/aa_example.phy" -m LG+F2+G4 -te "$HERE/data/aa_example_lg.nwk" -nt 1 -seed $SEED --prefix fault -redo --analytical-gradients --ag-stats > fault.stdout 2>&1 || { echo "  faulted run failed"; return 1; }
    grep -q "analytic gradient failed (AG_TEST_FAULT_EDGE" fault.stdout || { echo "  fault was not reported"; return 1; }
    "$BIN" -s "$EX/aa_example.phy" -m LG+F2+G4 -te "$HERE/data/aa_example_lg.nwk" -nt 1 -seed $SEED --prefix clean -redo --analytical-gradients > clean.stdout 2>&1 || return 1
    echo "  logl faulted=$(logl_of fault.iqtree) clean=$(logl_of clean.iqtree) $(grep -o "fd_fallbacks=[0-9]*" fault.stdout | head -1)"
    grep -q "fd_fallbacks=1 " fault.stdout || { echo "  expected exactly one fallback (the injected fault is one-shot)"; return 1; }
    within "$(logl_of fault.iqtree)" "$(logl_of clean.iqtree)" 0.5 || { echo "  faulted run ended far from the clean run"; return 1; }
    echo "  fault injection: PASS"
}

robust_x_threads_partitions() {   # X5: -Q at 4 threads (partitions enter the hook concurrently) is deterministic
    local i
    for i in 1 2; do
        "$BIN" -s "$WD/turtle_aa.fasta" -Q "$WD/turtle_aa.nex" -m LG+F+G4 -nt 4 -seed $SEED --prefix r$i -redo --analytical-gradients > r$i.stdout 2>&1 || { echo "  run $i failed"; return 1; }
    done
    echo "  logl run1=$(logl_of r1.iqtree) run2=$(logl_of r2.iqtree)"
    within "$(logl_of r1.iqtree)" "$(logl_of r2.iqtree)" 1e-3 || { echo "  runs differ"; return 1; }
    echo "  concurrent partitions: PASS"
}

ROBUST=(x_abort_resume x_ckp_crosscompat x_modelfinder x_pmsf x_gammai_restart x_fault x_threads_partitions)

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
if [ "$SUITE" = "quality" ] || [ "$SUITE" = "all" ]; then
    for e in "${QUALITY[@]}"; do
        IFS='|' read -r id tol args <<< "$e"
        ORDER+=("$id"); throttle; run_quality "$id" "$tol" "$args" &
    done
    ORDER+=("q_partitions"); throttle; run_partitions &
    ORDER+=("q_recovery"); throttle; run_recovery &
    ORDER+=("q_starts"); throttle; run_starts &
fi
if [ "$SUITE" = "robust" ] || [ "$SUITE" = "all" ]; then
    for id in "${ROBUST[@]}"; do ORDER+=("$id"); throttle; robust_case "$id" & done
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
