#!/usr/bin/env bash
#
# Run a sequence of independent seeds in parallel and merge the results.
#
#   tools/run_seeds.sh <nseeds> [njobs]
#
#   nseeds  how many seeds to run
#   njobs   concurrent jobs (default: all cores)
#
# Each job writes <outdir>/LO_<seed>.dat, taking outdir from the run card, and
# the merged result lands beside them on the untagged <outdir>/LO.dat -- a name
# no per-seed file can ever have, so a merged report cannot be fed back into its
# own merge.
#
# A seed that fails does not abort the campaign: the failures are listed and
# whatever succeeded is still merged. The exit status is non-zero if any seed
# failed, so a wrapper can notice.
#
# SEEDS ARE 1001, 1002, 1003, ... -- the run index, as in the KALEU tree.
#
# That convention was worth checking rather than assuming, because this build's
# generator is BRM48 (RANMAR), not the RANLUX the KALEU tree seeds. BRM48I
# splits its seed as
#
#     IJ = seed/30082 ;  KL = seed - 30082*IJ
#     I = MOD(IJ/177,177)+2 ;  J = MOD(IJ,177)+2
#     K = MOD(KL/169,178)+1 ;  L = MOD(KL,169)
#
# and fills its 97-word table from M = MOD(MOD(I*J,179)*K,179), a recurrence in
# I, J and K alone. Every seed below 30082 gives IJ = 0, so seeds 1001, 1002
# and 1003 all get I = J = 2 and K = 6, differing ONLY in L. On paper that is a
# weak corner of RANMAR's seeding: its independent-sequence guarantee is
# claimed across the whole IJKL range, not inside one IJ cell.
#
# MEASURED, AND IT DOES NOT MATTER. Cross-correlating the raw streams of two
# seeds over 1e7 draws at lag 0, and over 2e6 draws scanning lags 0..400, puts
# every pair at the noise floor:
#
#     1001 vs 1002    max |r| = 2.26e-3      (consecutive, shared I,J,K)
#     1001 vs 1003    max |r| = 2.30e-3      (consecutive)
#     1001 vs 31083   max |r| = 2.17e-3      (one IJ cell apart)
#     40   vs 1002    max |r| = 2.21e-3      (unrelated)
#
# against a 1/sqrt(N) floor of 7.1e-4 -- and 2.2e-3 is just what the LARGEST of
# 401 noise samples looks like. Consecutive seeds are indistinguishable from
# unrelated ones. The 48-bit table fill evidently washes the shared recurrence
# out. So the readable convention stands, and a campaign is reproducible from
# the run index alone.
#
# SEED_STRIDE remains an escape hatch: set it to 30082 to put every seed in its
# own IJ cell if a future test ever finds something this one could not see.
#
set -uo pipefail

SEED_BASE=${SEED_BASE:-1001}
SEED_STRIDE=${SEED_STRIDE:-1}

die() { echo "ERROR: $*" >&2; exit 1; }

[[ $# -ge 1 ]] || die "usage: $0 <nseeds> [njobs]"

NSEEDS=$1
NJOBS=${2:-$(nproc)}

[[ "$NSEEDS" =~ ^[0-9]+$ && "$NSEEDS" -ge 1 ]] || die "nseeds must be a positive integer"
[[ "$NJOBS"  =~ ^[0-9]+$ && "$NJOBS"  -ge 1 ]] || die "njobs must be a positive integer"

# The script lives in tools/, but runLO reads its card as '../input.input', so
# the binary only runs with LO_check_DY as the working directory -- and the
# card's [output] dir is written relative to that same place.
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
RUNDIR="$ROOT/LO_check_DY"
CARD="$ROOT/input.input"

[[ -x "$RUNDIR/runLO" ]] || die "$RUNDIR/runLO not found -- run make first"
[[ -f "$CARD" ]]         || die "$CARD not found"

cd "$RUNDIR"

# Read the card rather than duplicating its values, so the launcher and a bare
# run always agree on where results go and what they are called.
card_value() {
    sed -n "s/^[[:space:]]*$1[[:space:]]*=[[:space:]]*\([^#!]*\).*/\1/p" "$CARD" \
        | head -1 | xargs
}

OUTDIR=$(card_value dir)
OUTFILE=$(card_value file)
OUTDIR=${OUTDIR:-../summary/run_00}
OUTFILE=${OUTFILE:-LO.dat}

STEM=${OUTFILE%.*}
EXT=${OUTFILE##*.}
[[ "$STEM" == "$OUTFILE" ]] && EXT=""          # a name with no extension
SUFFIX=${EXT:+.$EXT}

LOGDIR="$OUTDIR/logs"
mkdir -p "$OUTDIR" "$LOGDIR" || die "cannot create $OUTDIR"

# Build the seed list once and reuse it everywhere: for launching, for the
# survivor check, and for the merge input list. Nothing globs, so a stray file
# in the output directory cannot get swept into the merge.
SEEDS=()
for ((k = 0; k < NSEEDS; k++)); do
    SEEDS+=( $(( SEED_BASE + k * SEED_STRIDE )) )
done

echo "seeds   : ${SEEDS[*]}"
echo "stride  : $SEED_STRIDE"
echo "jobs    : $NJOBS"
echo "output  : $OUTDIR/${STEM}_<seed>${SUFFIX}  ->  $OUTDIR/${STEM}${SUFFIX}"
echo

# Each job records its own outcome; the fan-out must not stop at the first
# failure.
run_one() {
    local seed=$1
    if ./runLO "$seed" "$OUTDIR" > "$LOGDIR/${STEM}_${seed}.log" 2>&1; then
        echo "  seed $seed  ok"
    else
        echo "  seed $seed  FAILED (see $LOGDIR/${STEM}_${seed}.log)"
    fi
}
export -f run_one
export OUTDIR LOGDIR STEM

if command -v parallel >/dev/null 2>&1; then
    printf '%s\n' "${SEEDS[@]}" | parallel -j "$NJOBS" --will-cite run_one {}
else
    echo "(GNU parallel not found, using xargs -P)"
    printf '%s\n' "${SEEDS[@]}" | xargs -P "$NJOBS" -I{} bash -c 'run_one "$@"' _ {}
fi

# Judge success by the reports actually on disk rather than by exit codes alone
# -- a job killed by the scheduler may leave no report at all.
ok=(); bad=()
for seed in "${SEEDS[@]}"; do
    f="$OUTDIR/${STEM}_${seed}${SUFFIX}"
    # -f as well as -s: a directory also satisfies -s, so -s alone would count
    # a blocked output path as a successful run.
    if [[ -f "$f" && -s "$f" ]]; then
        ok+=("$f")
    else
        bad+=("$seed")
    fi
done

echo
echo "completed: ${#ok[@]}/$NSEEDS"
if (( ${#bad[@]} )); then
    echo "failed   : ${bad[*]}"
    echo "           logs in $LOGDIR"
fi

if (( ${#ok[@]} == 0 )); then
    die "every seed failed -- nothing to merge"
fi

echo
echo "merging ${#ok[@]} report(s) ..."
"$ROOT/tools/merge_runs.py" "${ok[@]}" -o "$OUTDIR/${STEM}${SUFFIX}" || die "merge failed"

echo
echo "result: $OUTDIR/${STEM}${SUFFIX}"
echo "  tools/extract_data.py $OUTDIR/${STEM}${SUFFIX} tot"

(( ${#bad[@]} == 0 ))
