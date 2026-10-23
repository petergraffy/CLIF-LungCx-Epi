#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
export CLIF_FEDERATED_OUTPUT="${CLIF_REVISION_OUTPUT:-output/run_revision_gaps_20261022}/sensitivities/absorbing_v2"
pids=()
for site in UCMC NU RUSH; do
 mkdir -p "$CLIF_FEDERATED_OUTPUT/sites/$site/local"
 Rscript diagnostics/revision20261022/site_discover_absorbing.R "$site" > "$CLIF_FEDERATED_OUTPUT/sites/$site/local/discovery.log" 2>&1 &
 pids+=("$!")
done
for pid in "${pids[@]}"; do wait "$pid"; done
Rscript diagnostics/federated/stage_discovery.R
Rscript diagnostics/revision20261022/central_discover_absorbing.R none
for site in UCMC NU RUSH; do Rscript diagnostics/revision20261022/site_absorbing_summary.R "$site"; done
