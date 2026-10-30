#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
POOL_PRIMARY="${CLIF_POOLED_PRIMARY:-output/run_federated_rehearsal_20261026}"
POOL_OUTPUT="${CLIF_POOLED_SURVIVAL_OUTPUT:-output/run_pooled_survival_20261030_v2}"
if [ -e "$POOL_OUTPUT" ]; then echo 'Choose a fresh CLIF_POOLED_SURVIVAL_OUTPUT' >&2; exit 1; fi
IFS=',' read -r -a POOL_SITES <<< "${CLIF_POOLED_SITES:-UCMC,NU,RUSH}"
POOL_TIME_RETURNS=()
POOL_RISK_RETURNS=()
for POOL_SITE in "${POOL_SITES[@]}"; do
  Rscript diagnostics/pooled_survival/site_times.R "$POOL_SITE" "$POOL_PRIMARY" "$POOL_OUTPUT/sites/$POOL_SITE/times"
  POOL_TIME_RETURNS+=("$POOL_OUTPUT/sites/$POOL_SITE/times")
done
Rscript diagnostics/pooled_survival/central_grid.R "$POOL_OUTPUT/coordinator/grid" "${POOL_TIME_RETURNS[@]}"
for POOL_SITE in "${POOL_SITES[@]}"; do
  Rscript diagnostics/pooled_survival/site_risksets.R "$POOL_SITE" "$POOL_PRIMARY" "$POOL_OUTPUT/coordinator/grid" "$POOL_OUTPUT/sites/$POOL_SITE/risksets"
  POOL_RISK_RETURNS+=("$POOL_OUTPUT/sites/$POOL_SITE/risksets")
done
Rscript diagnostics/pooled_survival/central_pool.R "$POOL_OUTPUT/reports" "$POOL_OUTPUT/coordinator/grid" "${POOL_RISK_RETURNS[@]}"
CLIF_POOLED_SURVIVAL_OUTPUT="$POOL_OUTPUT" python3 diagnostics/pooled_survival/plot.py
