#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
run_root="${CLIF_FEDERATED_OUTPUT:-output/run_federated_rehearsal_20261022}"
if [[ -e "$run_root" ]]; then
  echo "Choose a fresh CLIF_FEDERATED_OUTPUT; the existing rehearsal is preserved." >&2
  exit 1
fi
export CLIF_FEDERATED_OUTPUT="$run_root"
mkdir -p "$run_root"
for site in UCMC NU RUSH; do mkdir -p "$run_root/sites/$site/local"; done
pids=()
for site in UCMC NU RUSH; do
  Rscript diagnostics/federated/site_discover.R "$site" > "$run_root/sites/$site/local/discovery.log" 2>&1 &
  pids+=("$!")
done
for pid in "${pids[@]}"; do wait "$pid"; done
Rscript diagnostics/federated/stage_discovery.R
for hold in none UCMC NU RUSH; do Rscript diagnostics/federated/central_discover.R "$hold"; done
Rscript tests/federated_rehearsal.R > "$run_root/federated_tests.log" 2>&1
pids=()
for site in UCMC NU RUSH; do
  (
    for version in full loo_UCMC loo_NU loo_RUSH; do
      Rscript diagnostics/federated/site_outcomes.R "$site" "$version" > "$run_root/sites/$site/local/outcomes_$version.log" 2>&1
    done
  ) &
  pids+=("$!")
done
for pid in "${pids[@]}"; do wait "$pid"; done
Rscript diagnostics/federated/stage_outcomes.R
for version in full loo_UCMC loo_NU loo_RUSH; do Rscript diagnostics/federated/central_pool.R "$version"; done
pids=()
for site in UCMC NU RUSH; do
  (
    Rscript diagnostics/federated/site_transport.R "$site" > "$run_root/sites/$site/local/transport.log" 2>&1
    Rscript diagnostics/federated/site_comparators.R "$site" > "$run_root/sites/$site/local/comparators.log" 2>&1
    Rscript diagnostics/federated/site_risks.R "$site" > "$run_root/sites/$site/local/risks.log" 2>&1
  ) &
  pids+=("$!")
done
for pid in "${pids[@]}"; do wait "$pid"; done
Rscript diagnostics/federated/central_comparators.R
pids=()
for site in UCMC NU RUSH; do
  (
    Rscript diagnostics/federated/site_transport_extended.R "$site" > "$run_root/sites/$site/local/transport_extended.log" 2>&1
    Rscript diagnostics/federated/site_quality.R "$site" > "$run_root/sites/$site/local/quality.log" 2>&1
  ) &
  pids+=("$!")
done
for pid in "${pids[@]}"; do wait "$pid"; done
Rscript diagnostics/federated/central_finalize.R > "$run_root/finalize.log" 2>&1
Rscript tests/federated_immutability.R > "$run_root/immutability_tests.log" 2>&1
python3 diagnostics/federated/write_report.py
