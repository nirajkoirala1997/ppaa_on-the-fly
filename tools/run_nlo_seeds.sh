#!/usr/bin/env bash
# Launch the four NLO contributions for a common sequence of seeds.
#
#   tools/run_nlo_seeds.sh <nseeds> [njobs] [--force] [--merge-only|--partial]
#
# --partial runs nothing and touches no finished report or campaign marker:
# it rebuilds every seed's report from its checkpoint_<seed>.dat
# (tools/checkpoint_to_report.py -- the checkpoint is rewritten after every
# production iteration, so this works on a running or killed campaign), merges
# them into <contrib>/partial/, and sums partial/nlo.dat.
# A seed with no checkpoint but a finished report uses the report.
#
# Unlike invoking run_seeds.sh four times, NJOBS is a global process limit.
# The live state is recorded under <output>/.nlo_campaign/ and is consumed by
# tools/monitor_nlo.py.  The calculation logs and reports remain in their
# normal contribution directories.
set -uo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }
usage="usage: $0 <nseeds> [njobs] [--force] [--merge-only|--partial]"

FORCE=0
MERGE_ONLY=0
PARTIAL=0
pos=()
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1 ;;
        --merge-only) MERGE_ONLY=1 ;;
        --partial) PARTIAL=1; MERGE_ONLY=1 ;;
        -h|--help) echo "$usage"; exit 0 ;;
        -*) die "unknown option: $arg" ;;
        *) pos+=("$arg") ;;
    esac
done
[[ ${#pos[@]} -ge 1 ]] || die "$usage"
NSEEDS=${pos[0]}
NJOBS=${pos[1]:-$(nproc)}
[[ "$NSEEDS" =~ ^[0-9]+$ && "$NSEEDS" -gt 0 ]] || die "nseeds must be positive"
[[ "$NJOBS" =~ ^[0-9]+$ && "$NJOBS" -gt 0 ]] || die "njobs must be positive"

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CARD="$REPO/input.input"
[[ -f "$CARD" ]] || die "missing card: $CARD"
card_value() {
    sed -n "s/^[[:space:]]*$1[[:space:]]*=[[:space:]]*\([^#!]*\).*/\1/p" "$CARD" | head -1 | xargs
}
BASE_REL=$(card_value dir)
BASE_REL=${BASE_REL:-../summary/run_00}
SEED_BASE=${SEED_BASE:-$(card_value seed)}
[[ "$SEED_BASE" =~ ^[0-9]+$ && "$SEED_BASE" -gt 0 ]] || die "invalid [run] seed: $SEED_BASE"
SEED_STRIDE=${SEED_STRIDE:-1}
[[ "$SEED_STRIDE" =~ ^[0-9]+$ && "$SEED_STRIDE" -gt 0 ]] || die "SEED_STRIDE must be positive"

# Every study resolves the shared output directory relative to itself.
BASE_OUT=$(cd "$REPO/lo_sa" && realpath -m "$BASE_REL")
CAMPAIGN="$BASE_OUT/.nlo_campaign"
STATUS="$CAMPAIGN/status"
mkdir -p "$STATUS" || die "cannot create $STATUS"

SEEDS=()
for ((i=0; i<NSEEDS; ++i)); do SEEDS+=( $((SEED_BASE + i * SEED_STRIDE)) ); done
TARGETS=(lo virt kp rs)

declare -A RUNDIR BIN SUBDIR STEM CHANNEL
RUNDIR[lo]="$REPO/lo_sa";     BIN[lo]=runLO;     SUBDIR[lo]=born; STEM[lo]=born; CHANNEL[lo]=''
RUNDIR[virt]="$REPO/virt_sa"; BIN[virt]=runVirt; SUBDIR[virt]=virt; STEM[virt]=virt; CHANNEL[virt]=''
RUNDIR[kp]="$REPO/kp_sa";     BIN[kp]=runKP;     SUBDIR[kp]=kp; STEM[kp]=kp; CHANNEL[kp]=''
RUNDIR[rs]="$REPO/rs_sa";     BIN[rs]=runRS;     SUBDIR[rs]=rs; STEM[rs]=rs; CHANNEL[rs]=''

for target in "${TARGETS[@]}"; do
    [[ -x "${RUNDIR[$target]}/${BIN[$target]}" || $MERGE_ONLY -eq 1 ]] || \
        die "missing executable: ${RUNDIR[$target]}/${BIN[$target]}"
done

MANIFEST="$CAMPAIGN/campaign.info"
if (( ! PARTIAL )); then
{
    echo "# Written by run_nlo_seeds.sh; monitor_nlo.py reads this file."
    echo "base_out=$BASE_OUT"
    echo "seed_base=$SEED_BASE"
    echo "seed_stride=$SEED_STRIDE"
    echo "seeds=${SEEDS[*]}"
    echo "targets=${TARGETS[*]}"
    echo "nseeds=$NSEEDS"
    echo "njobs=$NJOBS"
    echo "card=$CARD"
    echo "started=$(date -Is)"
} > "$MANIFEST"
fi

set_state() {
    local target=$1 seed=$2 state=$3 dir="$STATUS/$seed"
    mkdir -p "$dir"
    rm -f "$dir/$target.queued" "$dir/$target.running" "$dir/$target.done" "$dir/$target.failed"
    : > "$dir/$target.$state"
}

(( PARTIAL )) || for seed in "${SEEDS[@]}"; do
    for target in "${TARGETS[@]}"; do
        report="$BASE_OUT/${SUBDIR[$target]}/${STEM[$target]}_${seed}.dat"
        if [[ -f "$report" && $FORCE -eq 0 && $MERGE_ONLY -eq 0 ]]; then
            die "report already exists: $report (use --force or --merge-only)"
        fi
        if [[ -s "$report" && $MERGE_ONLY -eq 1 ]]; then
            set_state "$target" "$seed" done
        else
            set_state "$target" "$seed" queued
        fi
    done
done

run_one() {
    local target=$1 seed=$2 outdir="$BASE_OUT/${SUBDIR[$target]}" stem="${STEM[$target]}"
    local logdir="$outdir/logs"
    local log="$logdir/${stem}_${seed}.log" report="$outdir/${stem}_${seed}.dat"
    set_state "$target" "$seed" running
    mkdir -p "$logdir"
    local args=("$seed" "$BASE_OUT")
    [[ -n "${CHANNEL[$target]}" ]] && args=("${CHANNEL[$target]}" "$seed" "$BASE_OUT")
    if ( cd "${RUNDIR[$target]}" && ./"${BIN[$target]}" "${args[@]}" ) >> "$log" 2>&1 && [[ -s "$report" ]]; then
        set_state "$target" "$seed" done
    else
        set_state "$target" "$seed" failed
        return 1
    fi
}

if (( ! MERGE_ONLY )); then
    echo "NLO campaign: ${#SEEDS[@]} seeds x ${#TARGETS[@]} contributions, global jobs=$NJOBS"
    echo "Live monitor: $REPO/tools/monitor_nlo.py $MANIFEST"
    echo "Interim NLO : $0 $NSEEDS --partial"
    failures=0
    declare -a pids=()
    for seed in "${SEEDS[@]}"; do
        for target in "${TARGETS[@]}"; do
            run_one "$target" "$seed" & pids+=("$!")
            if (( ${#pids[@]} >= NJOBS )); then
                wait -n || failures=1
                pids=("${pids[@]:1}")
            fi
        done
    done
    for pid in "${pids[@]}"; do wait "$pid" || failures=1; done
else
    failures=0
fi

# Merge, per contribution. Normally only FINISHED reports, so an
# incomplete/failed seed never enters a final result. With --partial, every
# seed's report is rebuilt from its checkpoint into <contrib>/partial/ and
# everything below lands under partial/ instead.
P=""; (( PARTIAL )) && P="partial/"
for target in "${TARGETS[@]}"; do
    dir="$BASE_OUT/${SUBDIR[$target]}"
    inputs=()
    for seed in "${SEEDS[@]}"; do
        f="$dir/${STEM[$target]}_${seed}.dat"
        if (( PARTIAL )); then
            ck="$dir/checkpoint_${seed}.dat"
            pf="$dir/partial/${STEM[$target]}_${seed}.dat"
            if [[ -s "$ck" ]] && "$REPO/tools/checkpoint_to_report.py" "$ck" -o "$pf" >/dev/null; then
                inputs+=("$pf")
            elif [[ -s "$f" ]]; then
                inputs+=("$f")
            elif [[ -e "$ck" ]]; then
                echo "  $target seed $seed: checkpoint unreadable (mid-write?), skipped" >&2
            fi
        else
            [[ -s "$f" ]] && inputs+=("$f")
        fi
    done
    if (( ${#inputs[@]} )); then
        mkdir -p "$dir/$P"
        "$REPO/tools/merge_runs.py" "${inputs[@]}" -o "$dir/$P${STEM[$target]}.dat" --force || failures=1
    else
        echo "  $target: no report${P:+ or checkpoint} for any seed" >&2
        failures=1
    fi
done
# --force on every merge: without it merge_runs.py keeps an existing output
# and only warns, so a re-merge after more seeds finished would silently
# leave the old merged file in place.
PARTS=("$BASE_OUT/born/${P}born.dat" "$BASE_OUT/virt/${P}virt.dat" "$BASE_OUT/kp/${P}kp.dat" "$BASE_OUT/rs/${P}rs.dat")
if [[ -s "${PARTS[0]}" && -s "${PARTS[1]}" && -s "${PARTS[2]}" && -s "${PARTS[3]}" ]]; then
    mkdir -p "$BASE_OUT/$P"
    "$REPO/tools/merge_runs.py" --sum-components "${PARTS[@]}" -o "$BASE_OUT/${P}nlo.dat" --force || failures=1
    echo "NLO: $BASE_OUT/${P}nlo.dat"
else
    echo "NLO not summed: one of born/virt/kp/rs has no ${P:-merged }report" >&2
fi
if (( ! PARTIAL )); then
    echo "finished=$(date -Is)" >> "$MANIFEST"
    echo "Campaign manifest: $MANIFEST"
fi
exit "$failures"
