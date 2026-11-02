#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
PRIMARY="${CLIF_FOUR_COX_PRIMARY:-output/run_nearest_assignment_20261030/primary}"
REVISION="${CLIF_FOUR_COX_REVISION:-output/run_nearest_assignment_20261030/revision}"
OUTPUT="${CLIF_FOUR_COX_OUTPUT:-output/run_nearest_assignment_20261030/four_endpoint_cox}"
if [ -e "$OUTPUT" ]; then echo 'Choose a fresh four-endpoint output' >&2; exit 1; fi
mkdir -p "$OUTPUT/logs"
PIDS=()
for SITE in UCMC NU RUSH; do
 Rscript diagnostics/four_endpoint_cox/site_capture.R "$SITE" "$PRIMARY" "$REVISION" "$OUTPUT/sites/$SITE" > "$OUTPUT/logs/$SITE.log" 2>&1 &
 PIDS+=("$!")
done
for PID in "${PIDS[@]}"; do wait "$PID"; done
Rscript diagnostics/four_endpoint_cox/collect.R "$OUTPUT/reports" "$OUTPUT/sites/UCMC/share" "$OUTPUT/sites/NU/share" "$OUTPUT/sites/RUSH/share"
