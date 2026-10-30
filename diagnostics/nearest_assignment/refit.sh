#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
NEAREST_PARENT="$1"
NEAREST_OUTPUT="$2"
NEAREST_VARIANT="${3:-}"
export CLIF_FEDERATED_OUTPUT="$NEAREST_OUTPUT"
if [[ -n "$NEAREST_VARIANT" ]]; then
 export CLIF_REVISION_VARIANT="$NEAREST_VARIANT"
 Rscript diagnostics/nearest_assignment/prepare.R "$NEAREST_PARENT" "$NEAREST_OUTPUT" "$NEAREST_VARIANT"
 NEAREST_WORKER=diagnostics/revision20261022/site_outcomes.R
 NEAREST_POOL=diagnostics/revision20261022/central_pool.R
else
 Rscript diagnostics/nearest_assignment/prepare.R "$NEAREST_PARENT" "$NEAREST_OUTPUT"
 NEAREST_WORKER=diagnostics/federated/site_outcomes.R
 NEAREST_POOL=diagnostics/federated/central_pool.R
fi
mkdir -p "$NEAREST_OUTPUT/logs"
NEAREST_PIDS=()
for NEAREST_SITE in UCMC NU RUSH; do
 (
  for NEAREST_VERSION in full loo_UCMC loo_NU loo_RUSH; do
   Rscript "$NEAREST_WORKER" "$NEAREST_SITE" "$NEAREST_VERSION" > "$NEAREST_OUTPUT/logs/${NEAREST_SITE}_${NEAREST_VERSION}.log" 2>&1
  done
 ) &
 NEAREST_PIDS+=("$!")
done
for NEAREST_PID in "${NEAREST_PIDS[@]}"; do wait "$NEAREST_PID"; done
Rscript diagnostics/federated/stage_outcomes.R
for NEAREST_VERSION in full loo_UCMC loo_NU loo_RUSH; do Rscript "$NEAREST_POOL" "$NEAREST_VERSION"; done
NEAREST_PIDS=()
for NEAREST_SITE in UCMC NU RUSH; do
 (
  Rscript diagnostics/federated/site_transport.R "$NEAREST_SITE"
  if [[ -z "$NEAREST_VARIANT" ]]; then
   Rscript diagnostics/federated/site_comparators.R "$NEAREST_SITE"
   Rscript diagnostics/federated/site_risks.R "$NEAREST_SITE"
  fi
 ) > "$NEAREST_OUTPUT/logs/${NEAREST_SITE}_transport.log" 2>&1 &
 NEAREST_PIDS+=("$!")
done
for NEAREST_PID in "${NEAREST_PIDS[@]}"; do wait "$NEAREST_PID"; done
if [[ -z "$NEAREST_VARIANT" ]]; then
 Rscript diagnostics/federated/stage_discovery.R
 Rscript diagnostics/federated/central_comparators.R
 NEAREST_PIDS=()
 for NEAREST_SITE in UCMC NU RUSH; do
  (
   Rscript diagnostics/federated/site_transport_extended.R "$NEAREST_SITE"
   Rscript diagnostics/federated/site_quality.R "$NEAREST_SITE"
  ) > "$NEAREST_OUTPUT/logs/${NEAREST_SITE}_extended.log" 2>&1 &
  NEAREST_PIDS+=("$!")
 done
 for NEAREST_PID in "${NEAREST_PIDS[@]}"; do wait "$NEAREST_PID"; done
 Rscript diagnostics/federated/central_finalize.R > "$NEAREST_OUTPUT/logs/finalize.log" 2>&1
 Rscript tests/federated_rehearsal.R > "$NEAREST_OUTPUT/logs/federated_tests.log" 2>&1
 Rscript tests/federated_immutability.R > "$NEAREST_OUTPUT/logs/immutability_tests.log" 2>&1
 python3 diagnostics/federated/write_report.py
fi
