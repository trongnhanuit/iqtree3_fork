#!/bin/bash
# Default-path identity gate for --analytical-gradients work.
#
# Runs a fixed command list with a reference binary (built from the merge-base
# with master, untouched by this branch) and a new binary (this branch), at
# -nt 1 with a fixed seed, and requires the outputs to be identical apart from
# volatile lines (timestamps, timings, host, binary path).
#
# Usage:
#   test_scripts/ag/identity.sh <ref_binary> <new_binary> [out_dir] [--full] [--expect-flag-differs] [-j N]
#
# Modes:
#   A  ref (flag off)  vs new (flag off): strict identity  (always)
#   B  new (flag off)  vs new (flag on):  identity until the flag changes
#      behaviour; pass --expect-flag-differs (from Stage 3 on) to require that
#      at least one flagged run DOES differ (comparator self-test).
#
# Commands run in parallel, at most N at a time (default: half the cores, so a
# shared machine is never more than 50% loaded). Each run is single-threaded.
# The default list holds the fast cases; --full adds the slow ones.
#
# Exit code 0 = gate passed, 1 = a difference/self-test failure, 2 = usage.

set -u

REF_BIN="${1:-}"
NEW_BIN="${2:-}"
OUT_DIR="${3:-identity_out}"
FULL=0
EXPECT_FLAG_DIFFERS=0
MAXJOBS=0
argv=("$@")
for ((k=0; k<${#argv[@]}; k++)); do
    case "${argv[$k]}" in
        --full) FULL=1 ;;
        --expect-flag-differs) EXPECT_FLAG_DIFFERS=1 ;;
        -j) MAXJOBS="${argv[$((k+1))]:-0}" ;;
    esac
done
if [ -z "$REF_BIN" ] || [ -z "$NEW_BIN" ] || [ ! -x "$REF_BIN" ] || [ ! -x "$NEW_BIN" ]; then
    echo "usage: $0 <ref_binary> <new_binary> [out_dir] [--full] [--expect-flag-differs] [-j N]" >&2
    exit 2
fi
if [ "$MAXJOBS" -le 0 ]; then
    NCPU=$( (nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2) )
    MAXJOBS=$(( NCPU / 2 )); [ "$MAXJOBS" -lt 1 ] && MAXJOBS=1
fi

# runs happen inside per-run output directories, so binaries must be absolute
REF_BIN=$(cd "$(dirname "$REF_BIN")" && pwd)/$(basename "$REF_BIN")
NEW_BIN=$(cd "$(dirname "$NEW_BIN")" && pwd)/$(basename "$NEW_BIN")

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
SED_FILTER="$HERE/filter_volatile.sed"
SEED=73073
WD="$ROOT/test_scripts/test_data"
EX="$ROOT/example"

mkdir -p "$OUT_DIR"/ref "$OUT_DIR"/new "$OUT_DIR"/flag "$OUT_DIR"/report
OUT_DIR=$(cd "$OUT_DIR" && pwd)
REF_DIR="$OUT_DIR/ref"; NEW_DIR="$OUT_DIR/new"; FLAG_DIR="$OUT_DIR/flag"; REP_DIR="$OUT_DIR/report"
rm -f "$REP_DIR"/*

# Command list: "<id>|<arguments>" (prefix is added per run).
CMDS=(
  "ex_gtr_g|-s $EX/example.phy -m GTR+F+G4 -nt 1 -seed $SEED"
  "ex_gtr_fo_r4|-s $EX/example.phy -m GTR+FO+R4 -nt 1 -seed $SEED"
  "ex_mix_link|-s $EX/example.phy -m MIX{GTR+FO,GTR+FO} --link-exchange-rates -nt 1 -seed $SEED"
  "aa_lg_f2_g4_te|-s $EX/aa_example.phy -m LG+F2+G4 -te $HERE/data/aa_example_lg.nwk -nt 1 -seed $SEED"
  "ta_lg_c10_g4|-s $WD/turtle_aa.fasta -m LG+C10+G4 -nt 1 -seed $SEED"
  "ex_te_blfix|-s $EX/example.phy -m GTR+FO+G4 -te $HERE/data/example_gtr_g.nwk -blfix -nt 1 -seed $SEED"
  "ex_part_Q|-s $EX/example.phy -Q $EX/example.nex -m GTR+G -nt 1 -seed $SEED"
  "ex_part_p|-s $EX/example.phy -p $EX/example.nex -m GTR+G -nt 1 -seed $SEED"
)
if [ "$FULL" = "1" ]; then CMDS+=(
  "aa_lg_f4_r4|-s $EX/aa_example.phy -m LG+F4+R4 -nt 1 -seed $SEED"
  "ex_ufboot|-s $EX/example.phy -m GTR+F+I+G4 -B 1000 -nt 1 -seed $SEED"
  "ex_mfp|-s $EX/example.phy -m MFP -nt 1 -seed $SEED"
  "ta_mixfinder|-s $WD/turtle_aa.fasta -m MIX+MF -nt 1 -seed $SEED"
  "p_lg_f10_r4|-s $WD/prot_M126_27_269.phy -m LG+F10+R4 -nt 1 -seed $SEED"
); fi

run_one() {   # bin outdir id args...
    local bin="$1" dir="$2" id="$3"; shift 3
    ( cd "$dir" && "$bin" "$@" --prefix "$id" -redo -quiet > "$id.stdout" 2>&1 )
}

filtered() { sed -f "$SED_FILTER" "$1"; }

compare_pair() {   # dirA dirB id -> prints DIFF lines; returns 0 identical, 1 different
    local A="$1" B="$2" id="$3" rc=0 f n
    for f in treefile mldist contree ufboot splits.nex; do
        if [ -f "$A/$id.$f" ] || [ -f "$B/$id.$f" ]; then
            cmp -s "$A/$id.$f" "$B/$id.$f" || { echo "  DIFF $id.$f"; rc=1; }
        fi
    done
    for f in iqtree log; do
        if ! diff -q <(filtered "$A/$id.$f") <(filtered "$B/$id.$f") > /dev/null 2>&1; then
            echo "  DIFF $id.$f (after volatile-line filter)"
            diff <(filtered "$A/$id.$f") <(filtered "$B/$id.$f") | head -8; rc=1
        fi
    done
    if [ -f "$A/$id.ckp.gz" ] && [ -f "$B/$id.ckp.gz" ]; then
        if ! diff -q <(gunzip -c "$A/$id.ckp.gz" | grep -v -i -E 'time|command') \
                     <(gunzip -c "$B/$id.ckp.gz" | grep -v -i -E 'time|command') > /dev/null 2>&1; then
            echo "  DIFF $id.ckp.gz (after time/command filter)"; rc=1
        fi
    fi
    n=$(filtered "$A/$id.iqtree" 2>/dev/null | wc -l | tr -d ' ')
    if [ "${n:-0}" -lt 50 ]; then echo "  SELFTEST: filtered $id.iqtree has only ${n:-0} lines (filter too aggressive)"; rc=1; fi
    return $rc
}

# One job = the three runs of one command plus its comparisons. Writes
# report/<id>.txt and marker files A_FAIL/B_DIFF for the summary.
job() {
    local id="$1" args="$2" rep="$REP_DIR/$1.txt"
    {
        echo "== $id"
        # ref and new can run concurrently; the flagged run too.
        # shellcheck disable=SC2086
        run_one "$REF_BIN" "$REF_DIR" "$id" $args & local p1=$!
        # shellcheck disable=SC2086
        run_one "$NEW_BIN" "$NEW_DIR" "$id" $args & local p2=$!
        # shellcheck disable=SC2086
        run_one "$NEW_BIN" "$FLAG_DIR" "$id" $args --analytical-gradients & local p3=$!
        local ok=1
        wait $p1 || { echo "  ref run failed (see $REF_DIR/$id.stdout)"; ok=0; }
        wait $p2 || { echo "  new run failed (see $NEW_DIR/$id.stdout)"; ok=0; }
        wait $p3 || { echo "  flagged run failed (see $FLAG_DIR/$id.stdout)"; ok=0; }
        if [ "$ok" = "0" ]; then touch "$REP_DIR/$id.A_FAIL"; exit 0; fi
        if compare_pair "$REF_DIR" "$NEW_DIR" "$id"; then echo "  A: identical"; else echo "  A: FAILED"; touch "$REP_DIR/$id.A_FAIL"; fi
        if compare_pair "$NEW_DIR" "$FLAG_DIR" "$id" > /dev/null 2>&1; then
            echo "  B: flag on == flag off"
        else
            echo "  B: flag on != flag off"; touch "$REP_DIR/$id.B_DIFF"
        fi
    } > "$rep" 2>&1
}

# Throttle: each job launches 3 single-threaded runs, so allow MAXJOBS/3 jobs.
SLOTS=$(( MAXJOBS / 3 )); [ "$SLOTS" -lt 1 ] && SLOTS=1
echo "identity gate: ${#CMDS[@]} commands, up to $SLOTS concurrent (≤ $MAXJOBS processes)"
running=0
for entry in "${CMDS[@]}"; do
    id="${entry%%|*}"; args="${entry#*|}"
    job "$id" "$args" &
    running=$((running+1))
    if [ "$running" -ge "$SLOTS" ]; then wait -n 2>/dev/null || wait; running=$((running-1)); fi
done
wait

fail=0; flag_diff=0
for entry in "${CMDS[@]}"; do
    id="${entry%%|*}"
    cat "$REP_DIR/$id.txt"
    [ -f "$REP_DIR/$id.A_FAIL" ] && fail=1
    if [ -f "$REP_DIR/$id.B_DIFF" ]; then flag_diff=1; [ "$EXPECT_FLAG_DIFFERS" = "1" ] || fail=1; fi
done
if [ "$EXPECT_FLAG_DIFFERS" = "1" ] && [ "$flag_diff" = "0" ]; then
    echo "SELFTEST FAILED: --expect-flag-differs was given but no flagged run differed"; fail=1
fi
if [ "$fail" = "0" ]; then echo "IDENTITY GATE PASSED"; else echo "IDENTITY GATE FAILED"; fi
exit $fail
