#!/usr/bin/env bash
#
# Run a sequence of independent seeds in parallel and merge the results.
#
#   tools/run_seeds.sh <target> <nseeds> [njobs] [--merge-only] [--force]   # new campaign
#   tools/run_seeds.sh <target>                                            # reuse existing campaign
#
#   target  lo | kp | virt | rs
#   nseeds  how many seeds to run -- REQUIRED only when no NLO campaign
#           exists yet at the card's [output] dir. Once
#           tools/run_nlo_seeds.sh has created
#           <output_dir>/.nlo_campaign/campaign.info there, that manifest's
#           own seed list is authoritative: nseeds is then read from it, not
#           from the command line, and passing a conflicting nseeds is
#           rejected rather than silently overridden (pass none, or the
#           exact same count, to reuse the campaign).
#   njobs   concurrent jobs (default: all cores)
#
#   --merge-only   merge the reports already on disk, run nothing
#   --force        overwrite per-seed reports that already exist
#   --partial      run nothing; rebuild each seed's report from its
#                  checkpoint_<seed>.dat (tools/checkpoint_to_report.py) into
#                  <outdir>/partial/ and merge to <outdir>/partial/<stem>.dat.
#                  Works on a running or killed campaign; finished reports and
#                  campaign markers are never touched. A seed with no
#                  checkpoint but a finished report uses the report.
#
# When an existing campaign is reused, this script also updates that
# campaign's own per-seed marker files under
# <output_dir>/.nlo_campaign/status/<seed>/<target>.<state> (queued ->
# running -> done|failed, same convention tools/run_nlo_seeds.sh and
# tools/monitor_nlo.py use) -- but ONLY for the target selected here.
# Markers for the other three contributions are never touched, so this can
# resume/relaunch one target (e.g. lo) while the others keep running under
# the original campaign. Validation (target/dir resolution, manifest
# parsing, seed-list resolution, report-clash and executable checks) always
# runs to completion BEFORE any marker is written, so a rejected invocation
# never leaves a job falsely marked queued or running.
#
# ONE CARD, the root input.input, drives every target: each program reads
# ../input.input from its own directory, and appends its own fixed subdirectory
# (born/, kp/, virt/, rs/) to the card's [output] dir. Every program takes
# <binary> <seed> <outdir>, with outdir the shared run directory.
#
# A seed that fails does not abort the campaign: the failures are listed and
# whatever succeeded is still merged, with a non-zero exit status.
#
# SEEDS COUNT UP FROM THE CARD'S [run] seed -- 1001, 1002, 1003, ... -- the
# ONE seed every target starts from, since there is only one [run] seed in
# the whole file now (no more per-target numbering).
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

SEED_STRIDE=${SEED_STRIDE:-1}

die() { echo "ERROR: $*" >&2; exit 1; }

usage="usage: $0 <target> <nseeds> [njobs] [--merge-only|--partial] [--force]   # new campaign
       $0 <target>                                    # reuse existing campaign
       target: lo | kp | virt | rs"

MERGE_ONLY=0
PARTIAL=0
FORCE=0
POS=()
for arg in "$@"; do
    case "$arg" in
        --merge-only) MERGE_ONLY=1 ;;
        --partial)    PARTIAL=1; MERGE_ONLY=1 ;;
        --force)      FORCE=1 ;;
        -h|--help)    echo "$usage"; exit 0 ;;
        -*)           die "unknown option $arg"$'\n'"$usage" ;;
        *)            POS+=("$arg") ;;
    esac
done

[[ ${#POS[@]} -ge 1 ]] || die "$usage"

TARGET=${POS[0]}
POS=("${POS[@]:1}")
# POS now holds whatever is left of [nseeds] [njobs] -- resolved below, once
# we know whether an existing NLO campaign makes nseeds optional.

# tools/ sits at the repository root, so the merge script is found under
# REPO regardless of target. The binary is not: every one of these programs
# reads its card as '../input.input' relative to the working directory, so
# each only runs from its own study directory, RUNDIR -- and the card
# itself is the ONE file at input.input, always one level above
# RUNDIR.
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

CHANNEL=""

# SUBDIR/FILESTEM are each program's OWN fixed choice (its main.f90), not a
# card key any more -- see input.input's header. This table has to
# stay in step with those four main.f90 files by hand; there is nothing
# that checks it for you.
case "$TARGET" in
    lo)    RUNDIR="$REPO/lo_sa";   BIN=runLO   ; SUBDIR=born ; FILESTEM=born ;;
    kp)    RUNDIR="$REPO/kp_sa";   BIN=runKP   ; SUBDIR=kp   ; FILESTEM=kp   ;;
    virt)  RUNDIR="$REPO/virt_sa"; BIN=runVirt ; SUBDIR=virt ; FILESTEM=virt ;;
    rs)    RUNDIR="$REPO/rs_sa";   BIN=runRS   ; SUBDIR=rs   ; FILESTEM=rs   ;;
    *)     die "unknown target '$TARGET' -- expected lo, kp, virt or rs" ;;
esac
# TARGET doubles as the campaign marker/target name from here on -- it
# already matches tools/monitor_nlo.py's TARGETS tuple (lo, kp, virt, rsgg,
# rsqg, rsqq) exactly, so no separate name mapping is needed.

CARD="$RUNDIR/../input.input"

[[ -d "$RUNDIR" ]] || die "$RUNDIR not found"
[[ -f "$CARD" ]]   || die "$CARD not found"
if (( ! MERGE_ONLY )); then
    [[ -x "$RUNDIR/$BIN" ]] || \
        die "$BIN not found -- build it first:  cd $RUNDIR && make"
fi

cd "$RUNDIR"

# Read the card rather than duplicating its values, so the launcher and a bare
# run always agree on where results go and what they are called.
card_value() {
    sed -n "s/^[[:space:]]*$1[[:space:]]*=[[:space:]]*\([^#!]*\).*/\1/p" "$CARD" \
        | head -1 | xargs
}

# The card's [output] dir is the SHARED run directory now, resolved relative
# to RUNDIR (which is where the binary itself resolves it too) -- not
# already this target's leaf. BASE_OUTDIR is what gets passed as the
# outdir argument; the binary appends SUBDIR (and, for rs, CHANNEL) itself.
BASE_OUTDIR=$(card_value dir)
BASE_OUTDIR=${BASE_OUTDIR:-../summary/run_00}

# An NLO campaign (tools/run_nlo_seeds.sh) resolves this same [output] dir
# relative to lo_sa, one of the sibling study directories at the same
# depth as RUNDIR -- so resolving BASE_OUTDIR from here, in RUNDIR, lands on
# the identical absolute directory regardless of which target we are. That is
# what lets this script find the SAME campaign.info no matter which target
# is being resumed.
BASE_OUT_ABS=$(realpath -m "$BASE_OUTDIR")
CAMPAIGN_DIR="$BASE_OUT_ABS/.nlo_campaign"
MANIFEST="$CAMPAIGN_DIR/campaign.info"
STATUS="$CAMPAIGN_DIR/status"

IN_CAMPAIGN=0
declare -A CINFO
if [[ -f "$MANIFEST" ]]; then
    IN_CAMPAIGN=1
    while IFS='=' read -r k v; do
        [[ -z "$k" || "$k" == \#* ]] && continue
        CINFO["$k"]="$v"
    done < "$MANIFEST"
    [[ -n "${CINFO[seeds]:-}" && -n "${CINFO[nseeds]:-}" ]] || \
        die "malformed campaign manifest: $MANIFEST"
    [[ " ${CINFO[targets]:-} " == *" $TARGET "* ]] || \
        die "target '$TARGET' is not part of the campaign recorded in $MANIFEST (targets: ${CINFO[targets]:-none})"
fi

if (( IN_CAMPAIGN )); then
    # The campaign manifest is authoritative for the seed list: nseeds is
    # optional here, and if given it must match what the campaign was
    # created with -- it never silently overrides it.
    if [[ ${#POS[@]} -ge 1 ]]; then
        claimed=${POS[0]}
        [[ "$claimed" =~ ^[0-9]+$ ]] || die "nseeds must be a positive integer"
        [[ "$claimed" == "${CINFO[nseeds]}" ]] || \
            die "existing campaign at $MANIFEST was created with ${CINFO[nseeds]} seed(s); refusing to run $claimed instead. Omit nseeds to reuse the campaign's own seed list (${CINFO[seeds]})."
    fi
    NJOBS=${POS[1]:-$(nproc)}
    [[ "$NJOBS" =~ ^[0-9]+$ && "$NJOBS" -ge 1 ]] || die "njobs must be a positive integer"
    read -r -a SEEDS <<< "${CINFO[seeds]}"
    NSEEDS=${#SEEDS[@]}
    SEED_STRIDE=${CINFO[seed_stride]:-$SEED_STRIDE}
else
    [[ ${#POS[@]} -ge 1 ]] || \
        die "$usage"$'\n'"nseeds is required -- no existing campaign found at $MANIFEST"
    NSEEDS=${POS[0]}
    NJOBS=${POS[1]:-$(nproc)}
    [[ "$NSEEDS" =~ ^[0-9]+$ && "$NSEEDS" -ge 1 ]] || die "nseeds must be a positive integer"
    [[ "$NJOBS"  =~ ^[0-9]+$ && "$NJOBS"  -ge 1 ]] || die "njobs must be a positive integer"

    # The campaign starts at the card's own seed, so every target's numbering is
    # read from the same place a bare run would read it and none of them can
    # disagree -- there is only the one [run] seed in the file now. An explicit
    # SEED_BASE in the environment still wins.
    SEED_BASE=${SEED_BASE:-$(card_value seed)}
    [[ "$SEED_BASE" =~ ^[0-9]+$ && "$SEED_BASE" -ge 1 ]] || \
        die "[run] seed in $CARD is '$SEED_BASE', not a positive integer"

    # Build the seed list once and reuse it everywhere: for launching, for the
    # survivor check, and for the merge input list. Nothing globs, so a stray file
    # in the output directory cannot get swept into the merge.
    SEEDS=()
    for ((k = 0; k < NSEEDS; k++)); do
        SEEDS+=( $(( SEED_BASE + k * SEED_STRIDE )) )
    done
fi

SUFFIX=.dat
STEM="$FILESTEM"
OUTDIR="$BASE_OUTDIR/$SUBDIR"
if [[ -n "$CHANNEL" ]]; then
    OUTDIR="$OUTDIR/$CHANNEL"
    STEM="${STEM}_${CHANNEL}"
fi

LOGDIR="$OUTDIR/logs"
mkdir -p "$OUTDIR" "$LOGDIR" || die "cannot create $OUTDIR"

echo "target  : $TARGET   ($RUNDIR/$BIN)"
if (( IN_CAMPAIGN )); then
    echo "campaign: $MANIFEST (reusing its seed list)"
fi
echo "seeds   : ${SEEDS[*]}"
echo "stride  : $SEED_STRIDE"
echo "jobs    : $NJOBS"
echo "output  : $OUTDIR/${STEM}_<seed>${SUFFIX}  ->  $OUTDIR/${STEM}${SUFFIX}"
echo

# Only touches THIS target's own marker file in the seed's status directory
# -- never the other three contributions', which may still be running under
# the same campaign. Same marker names/precedence as tools/run_nlo_seeds.sh
# and tools/monitor_nlo.py: <target>.queued|running|done|failed, exactly one
# present at a time.
set_state() {
    local seed=$1 state=$2 dir="$STATUS/$seed"
    mkdir -p "$dir"
    rm -f "$dir/$TARGET.queued" "$dir/$TARGET.running" "$dir/$TARGET.done" "$dir/$TARGET.failed"
    : > "$dir/$TARGET.$state"
}
export -f set_state
export STATUS TARGET IN_CAMPAIGN

# Each job records its own outcome; the fan-out must not stop at the first
# failure.
run_one() {
    local seed=$1
    local args=("$seed" "$BASE_OUTDIR")
    [[ -n "$CHANNEL" ]] && args=("$CHANNEL" "$seed" "$BASE_OUTDIR")
    (( IN_CAMPAIGN )) && set_state "$seed" running
    if ./"$BIN" "${args[@]}" >> "$LOGDIR/${STEM}_${seed}.log" 2>&1 && \
       [[ -s "$OUTDIR/${STEM}_${seed}${SUFFIX}" ]]; then
        (( IN_CAMPAIGN )) && set_state "$seed" done
        echo "  seed $seed  ok"
    else
        (( IN_CAMPAIGN )) && set_state "$seed" failed
        echo "  seed $seed  FAILED (see $LOGDIR/${STEM}_${seed}.log)"
    fi
}
export -f run_one
export OUTDIR LOGDIR BASE_OUTDIR CHANNEL STEM BIN SUFFIX

if (( MERGE_ONLY )); then
    if (( PARTIAL )); then
        echo "partial: running nothing, rebuilding reports from checkpoints"
    else
        echo "merge-only: running nothing, merging what is already on disk"
    fi
else
    # A seed's report is not always a replica of its neighbours -- a directory
    # can hold runs made at different settings, each recorded under the seed
    # that produced it. Overwriting one silently would destroy that record, so
    # say what is in the way and make the caller ask for it.
    clash=()
    for seed in "${SEEDS[@]}"; do
        [[ -f "$OUTDIR/${STEM}_${seed}${SUFFIX}" ]] && clash+=("${STEM}_${seed}${SUFFIX}")
    done
    if (( ${#clash[@]} && ! FORCE )); then
        die "these reports already exist in $OUTDIR:"$'\n'"  ${clash[*]}"$'\n'\
"pass --force to overwrite them, SEED_BASE=<n> to start elsewhere, or
--merge-only to combine what is there without rerunning"
    fi

    # All validation is done -- resolve, campaign lookup, seed list, clash
    # check -- so only now, right before launch, do we touch any marker.
    if (( IN_CAMPAIGN )); then
        for seed in "${SEEDS[@]}"; do
            set_state "$seed" queued
        done
    fi

    if command -v parallel >/dev/null 2>&1; then
        printf '%s\n' "${SEEDS[@]}" | parallel -j "$NJOBS" --will-cite run_one {}
    else
        echo "(GNU parallel not found, using xargs -P)"
        printf '%s\n' "${SEEDS[@]}" | xargs -P "$NJOBS" -I{} bash -c 'run_one "$@"' _ {}
    fi
fi

# Judge success by the reports actually on disk rather than by exit codes alone
# -- a job killed by the scheduler may leave no report at all.
ok=(); bad=()
for seed in "${SEEDS[@]}"; do
    f="$OUTDIR/${STEM}_${seed}${SUFFIX}"
    # -f as well as -s: a directory also satisfies -s, so -s alone would count
    # a blocked output path as a successful run.
    if (( PARTIAL )); then
        ck="$OUTDIR/checkpoint_${seed}${SUFFIX}"
        pf="$OUTDIR/partial/${STEM}_${seed}${SUFFIX}"
        if [[ -s "$ck" ]] && "$REPO/tools/checkpoint_to_report.py" "$ck" -o "$pf"; then
            ok+=("$pf")
        elif [[ -f "$f" && -s "$f" ]]; then
            ok+=("$f")
        else
            bad+=("$seed")
        fi
    elif [[ -f "$f" && -s "$f" ]]; then
        ok+=("$f")
    else
        bad+=("$seed")
    fi
done
MERGED_DIR="$OUTDIR"
(( PARTIAL )) && MERGED_DIR="$OUTDIR/partial"

echo
if (( MERGE_ONLY )); then
    echo "found    : ${#ok[@]}/$NSEEDS"
    (( ${#bad[@]} )) && echo "missing  : ${bad[*]}"
    (( ${#ok[@]} )) || die "none of those seeds has a report in $OUTDIR"
else
    echo "completed: ${#ok[@]}/$NSEEDS"
    if (( ${#bad[@]} )); then
        echo "failed   : ${bad[*]}"
        echo "           logs in $LOGDIR"
    fi
    (( ${#ok[@]} )) || die "every seed failed -- nothing to merge"
fi

echo
echo "merging ${#ok[@]} report(s) ..."
mkdir -p "$MERGED_DIR"
"$REPO/tools/merge_runs.py" "${ok[@]}" -o "$MERGED_DIR/${STEM}${SUFFIX}" --force || die "merge failed"

# Printed from the repository root, since that is where tools/ is reachable --
# the working directory here is the study's, not the root's.
MERGED=$(cd "$MERGED_DIR" && pwd)/${STEM}${SUFFIX}
echo
echo "result: $MERGED"
echo "  tools/extract_data.py ${MERGED#$REPO/} tot"

(( ${#bad[@]} == 0 ))
