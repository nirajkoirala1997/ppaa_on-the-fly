#!/usr/bin/env bash
# The full diphoton campaign for the MATRIX comparison:
#   1. born, virt, kp, rs with the NLO set, into summary/run_01;
#   2. born again with the LO set, into summary/run_01_lo_pdf.
#
#   tools/run_diphoton_campaign.sh <nseeds> [njobs]
#
# Step 2 runs the runLO binary from a directory holding its own card, because
# the programs always read ../input.input.
set -uo pipefail

NSEEDS=${1:?usage: run_diphoton_campaign.sh <nseeds> [njobs]}
NJOBS=${2:-$(nproc)}
LO_PDF_SET=NNPDF30_lo_as_0118

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
if [[ -z "${LO_PDF_ONLY:-}" ]]; then
    "$REPO/tools/run_nlo_seeds.sh" "$NSEEDS" "$NJOBS"
    echo "NLO_CAMPAIGN_EXIT_STATUS $?"
fi

LO_BASE="${LO_BASE:-$REPO/summary/run_01_lo_pdf}"
# The card must not sit in $LO_BASE, where the program stamps its own copy.
mkdir -p "$LO_BASE/card/work" "$LO_BASE/born/logs"
sed "s/^set *=.*/set = $LO_PDF_SET/" "$REPO/input.input" > "$LO_BASE/card/input.input"

SEED_BASE=$(sed -n 's/^seed *= *\([0-9]*\).*/\1/p' "$REPO/input.input" | head -1)
cd "$LO_BASE/card/work" || exit 1
for ((index=0; index<NSEEDS; ++index)); do
    seed=$((SEED_BASE + index))
    "$REPO/lo_sa/runLO" "$seed" "$LO_BASE" > "$LO_BASE/born/logs/born_$seed.log" 2>&1 &
    if (( (index + 1) % NJOBS == 0 )); then wait; fi
done
wait

"$REPO/tools/merge_runs.py" "$LO_BASE"/born/born_*[0-9].dat -o "$LO_BASE/born/born.dat" --force
echo "LO_PDF_CAMPAIGN_EXIT_STATUS $?"
