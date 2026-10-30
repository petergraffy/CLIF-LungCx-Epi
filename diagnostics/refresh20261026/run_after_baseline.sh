#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
source diagnostics/refresh20261026/environment.sh
logs=output/run_three_site_refresh_20261026/logs
if [[ ! -f output/run_UCMC_ccx_box_20261026_joint/local/derivation.rds ]]; then
 echo 'Finish the fresh UCMC derivation before starting the federated refresh.' >&2; exit 1
fi
bash diagnostics/02_compare_waterfall.sh > "$logs/ucmc_waterfall.log" 2>&1
bash diagnostics/federated/run_rehearsal.sh > "$logs/federated.log" 2>&1
bash diagnostics/revision20261022/run.sh > "$logs/reviewer_sensitivities.log" 2>&1
Rscript diagnostics/refresh20261026/compare.R > "$logs/comparison.log" 2>&1
python3 diagnostics/refresh20261026/report.py
Rscript diagnostics/refresh20261026/audit.R > "$logs/final_audit.log" 2>&1
