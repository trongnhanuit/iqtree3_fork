#!/bin/bash
# Scalar-vs-SIMD equivalence gate for the per-chunk gradient accumulation
# (accumulateChunk vs accumulateChunkSSE, model/phylogradient.cpp).
#
# Runs a fixed command list twice with the same binary -- AG_SCALAR_ACCUM=1
# (forces the scalar path) vs unset (dispatches to SSE2 when vsize==2) -- at
# the initial point (--ag-gradient-check-only exits right after, so these
# stay fast even though it also runs its own, separate finite-difference
# check) and compares the dumped branch-length and theta gradients with
# compare_aggrad.py. `-lk SSE` forces vsize==2 regardless of the host's native
# kernel (an AVX/FMA/AVX512 machine -- most CI and dev hardware -- otherwise
# picks vsize 4/8 and accumulateChunkSSE would never actually run).
#
# Usage:
#   test_scripts/ag/simd_check.sh <binary> [out_dir]
#
# Exit code 0 = gate passed, 1 = a case failed, 2 = usage.

set -u

BIN="${1:-}"
OUT_DIR="${2:-simd_check_out}"
if [ -z "$BIN" ] || [ ! -x "$BIN" ]; then
    echo "usage: $0 <binary> [out_dir]" >&2
    exit 2
fi
BIN=$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
WD="$ROOT/test_scripts/test_data"
EX="$ROOT/example"
AGD="$HERE/data"
SEED=73073

mkdir -p "$OUT_DIR"
cd "$OUT_DIR" || exit 1

# label|args  -- covers estimated/fixed profiles, +I, -safe, and (turtle_aa,
# 407 patterns) an odd pattern count that exercises the last SIMD block's
# single real lane.
CMDS=(
    "dna_fo_g4|-s $EX/example.phy -m GTR+FO+G4 -te $AGD/example_gtr_g.nwk"
    "aa_fo_r4|-s $EX/aa_example.phy -m LG+FO+R4"
    "aa_c10_g4|-s $WD/turtle_aa.fasta -m LG+C10+G4"
    "dna_fo_i_g4|-s $EX/example.phy -m GTR+FO+I+G4"
    "dna_fo_g4_safe|-s $EX/example.phy -m GTR+FO+G4 -safe"
)

FAIL=0
for spec in "${CMDS[@]}"; do
    label="${spec%%|*}"
    args="${spec#*|}"
    rm -f "${label}"_scalar.* "${label}"_simd.*
    AG_SCALAR_ACCUM=1 "$BIN" -pre "${label}_scalar" $args -lk SSE \
        --ag-gradient-check-only --ag-dump-gradient -nt 1 -seed "$SEED" -redo > "${label}_scalar.stdout" 2>&1
    "$BIN" -pre "${label}_simd" $args -lk SSE \
        --ag-gradient-check-only --ag-dump-gradient -nt 1 -seed "$SEED" -redo > "${label}_simd.stdout" 2>&1
    if [ ! -s "${label}_scalar.aggrad.tsv" ] || [ ! -s "${label}_simd.aggrad.tsv" ]; then
        echo "== $label: MISSING aggrad.tsv"
        tail -3 "${label}_scalar.stdout" "${label}_simd.stdout"
        FAIL=1
        continue
    fi
    result=$(python3 "$HERE/compare_aggrad.py" "${label}_scalar.aggrad.tsv" "${label}_simd.aggrad.tsv")
    echo "== $label: $result"
    case "$result" in PASS*) ;; *) FAIL=1 ;; esac
done

if [ "$FAIL" -ne 0 ]; then
    echo "SIMD CHECK FAILED"
    exit 1
fi
echo "SIMD CHECK PASSED"
