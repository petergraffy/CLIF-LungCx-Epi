#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
export CLIF_APACHE_PRIMARY=output/run_nearest_assignment_20261030/primary
export CLIF_REVISION_PRIMARY_INPUT="$CLIF_APACHE_PRIMARY"
export CLIF_APACHE_OUTPUT=output/run_nearest_assignment_20261030/apache_adapted
Rscript diagnostics/nearest_assignment/prepare_apache.R
mkdir -p "$CLIF_APACHE_OUTPUT/logs"
NEAREST_PIDS=()
for NEAREST_SITE in UCMC NU RUSH; do
 Rscript diagnostics/apache_adapted/site_train.R "$NEAREST_SITE" > "$CLIF_APACHE_OUTPUT/logs/train_$NEAREST_SITE.log" 2>&1 &
 NEAREST_PIDS+=("$!")
done
for NEAREST_PID in "${NEAREST_PIDS[@]}"; do wait "$NEAREST_PID"; done
Rscript diagnostics/apache_adapted/central_train.R
NEAREST_PIDS=()
for NEAREST_SITE in UCMC NU RUSH; do
 (
  Rscript diagnostics/apache_adapted/site_evaluate.R "$NEAREST_SITE"
  Rscript diagnostics/apache_adapted/site_publish.R "$NEAREST_SITE"
  Rscript diagnostics/apache_adapted/site_audit.R "$NEAREST_SITE"
 ) > "$CLIF_APACHE_OUTPUT/logs/evaluate_$NEAREST_SITE.log" 2>&1 &
 NEAREST_PIDS+=("$!")
done
for NEAREST_PID in "${NEAREST_PIDS[@]}"; do wait "$NEAREST_PID"; done
Rscript diagnostics/apache_adapted/finalize.R > "$CLIF_APACHE_OUTPUT/logs/finalize.log" 2>&1
