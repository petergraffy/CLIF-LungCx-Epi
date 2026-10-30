#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
source diagnostics/refresh20261026/environment.sh
export CLIF_FEDERATED_OUTPUT=output/run_nearest_assignment_20261030/primary
export CLIF_REVISION_PRIMARY_INPUT="$CLIF_FEDERATED_OUTPUT"
export CLIF_REVISION_OUTPUT=output/run_nearest_assignment_20261030/revision
export CLIF_SITE_TABLES_NU=/Users/saborpete/Library/CloudStorage/Box-Box/CLIF-NU/CLIF_2_1_1
mkdir -p "$CLIF_REVISION_OUTPUT/logs"
NEAREST_PIDS=()
for NEAREST_SITE in UCMC NU RUSH; do
 (
  Rscript diagnostics/revision20261022/site_review.R "$NEAREST_SITE"
  Rscript diagnostics/revision20261022/site_descriptive_tables.R "$NEAREST_SITE"
  Rscript diagnostics/revision20261022/site_arf_effects.R "$NEAREST_SITE"
  Rscript diagnostics/revision20261022/site_severity_source_audit.R "$NEAREST_SITE"
  Rscript diagnostics/revision20261022/site_code_source_audit.R "$NEAREST_SITE"
  Rscript diagnostics/revision20261022/site_scores.R "$NEAREST_SITE"
 ) > "$CLIF_REVISION_OUTPUT/logs/review_$NEAREST_SITE.log" 2>&1 &
 NEAREST_PIDS+=("$!")
done
for NEAREST_PID in "${NEAREST_PIDS[@]}"; do wait "$NEAREST_PID"; done
Rscript diagnostics/revision20261022/central_scores.R
NEAREST_PIDS=()
for NEAREST_SITE in UCMC NU RUSH; do
 Rscript diagnostics/revision20261022/site_score_transport.R "$NEAREST_SITE" > "$CLIF_REVISION_OUTPUT/logs/score_transport_$NEAREST_SITE.log" 2>&1 &
 NEAREST_PIDS+=("$!")
done
for NEAREST_PID in "${NEAREST_PIDS[@]}"; do wait "$NEAREST_PID"; done
Rscript diagnostics/revision20261022/central_review.R
