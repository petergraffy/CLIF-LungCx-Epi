#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
revision_root="${CLIF_REVISION_OUTPUT:-output/run_revision_gaps_20261022}"
mkdir -p "$revision_root"
for variant in h24 h48 onset support_only low_flow carry0 carry12 raw_support; do
 export CLIF_REVISION_VARIANT="$variant"
 export CLIF_FEDERATED_OUTPUT="$revision_root/sensitivities/$variant"
 mkdir -p "$CLIF_FEDERATED_OUTPUT"
 pids=()
 for site in UCMC NU RUSH; do
  mkdir -p "$CLIF_FEDERATED_OUTPUT/sites/$site/local"
  Rscript diagnostics/revision20261022/site_discover.R "$site" > "$CLIF_FEDERATED_OUTPUT/sites/$site/local/discovery.log" 2>&1 &
  pids+=("$!")
 done
 for pid in "${pids[@]}"; do wait "$pid"; done
 Rscript diagnostics/federated/stage_discovery.R
 for hold in none UCMC NU RUSH; do Rscript diagnostics/revision20261022/central_discover.R "$hold"; done
 pids=()
 for site in UCMC NU RUSH; do
  (
   for version in full loo_UCMC loo_NU loo_RUSH; do
    Rscript diagnostics/revision20261022/site_outcomes.R "$site" "$version" > "$CLIF_FEDERATED_OUTPUT/sites/$site/local/outcomes_$version.log" 2>&1
   done
  ) &
  pids+=("$!")
 done
 for pid in "${pids[@]}"; do wait "$pid"; done
 Rscript diagnostics/federated/stage_outcomes.R
 for version in full loo_UCMC loo_NU loo_RUSH; do Rscript diagnostics/revision20261022/central_pool.R "$version"; done
 for site in UCMC NU RUSH; do Rscript diagnostics/federated/site_transport.R "$site"; done
 echo "Completed sensitivity $variant"
done
