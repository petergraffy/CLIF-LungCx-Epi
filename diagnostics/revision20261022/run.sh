#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
revision_root="${CLIF_REVISION_OUTPUT:-output/run_revision_gaps_20261022}"
if [[ -e "$revision_root" ]]; then
 echo "Choose a fresh CLIF_REVISION_OUTPUT; existing revision outputs are preserved." >&2
 exit 1
fi
export CLIF_REVISION_OUTPUT="$revision_root"
mkdir -p "$revision_root/logs"
Rscript tests/revision_gaps.R > "$revision_root/logs/design_tests.log" 2>&1
pids=()
for site in UCMC NU RUSH; do
 (
  Rscript diagnostics/revision20261022/site_review.R "$site"
  Rscript diagnostics/revision20261022/site_descriptive_tables.R "$site"
  Rscript diagnostics/revision20261022/site_arf_effects.R "$site"
  Rscript diagnostics/revision20261022/site_severity_source_audit.R "$site"
  Rscript diagnostics/revision20261022/site_code_source_audit.R "$site"
  Rscript diagnostics/revision20261022/site_scores.R "$site"
 ) > "$revision_root/logs/site_${site}.log" 2>&1 &
 pids+=("$!")
done
for pid in "${pids[@]}"; do wait "$pid"; done
Rscript diagnostics/revision20261022/central_scores.R > "$revision_root/logs/central_scores.log" 2>&1
pids=()
for site in UCMC NU RUSH; do
 Rscript diagnostics/revision20261022/site_score_transport.R "$site" > "$revision_root/logs/score_transport_${site}.log" 2>&1 &
 pids+=("$!")
done
for pid in "${pids[@]}"; do wait "$pid"; done
Rscript diagnostics/revision20261022/central_review.R > "$revision_root/logs/central_review.log" 2>&1
bash diagnostics/revision20261022/run_sensitivities.sh > "$revision_root/logs/sensitivities.log" 2>&1
bash diagnostics/revision20261022/run_absorbing.sh > "$revision_root/logs/absorbing.log" 2>&1
for site in UCMC NU RUSH; do
 Rscript diagnostics/revision20261022/site_audit.R "$site" > "$revision_root/logs/audit_${site}.log" 2>&1
 CLIF_FEDERATED_OUTPUT="$revision_root/sensitivities/absorbing_v2" Rscript diagnostics/revision20261022/site_absorbing_audit.R "$site"
done
Rscript diagnostics/revision20261022/central_audit.R > "$revision_root/logs/central_audit.log" 2>&1
apache_root="${CLIF_APACHE_OUTPUT:-$revision_root/apache_adapted}"
export CLIF_APACHE_OUTPUT="$apache_root"
export CLIF_APACHE_PRIMARY="${CLIF_REVISION_PRIMARY_INPUT:-output/run_federated_rehearsal_20261022}"
bash diagnostics/apache_adapted/run.sh > "$revision_root/logs/apache_adapted.log" 2>&1
for site in UCMC NU RUSH; do
 Rscript diagnostics/site_return/capture.R "$site" "${CLIF_REVISION_PRIMARY_INPUT:-output/run_federated_rehearsal_20261022}" "$revision_root" "$revision_root/site_returns/$site" > "$revision_root/logs/supplement_${site}.log" 2>&1
done
Rscript diagnostics/site_return/collect.R "$revision_root/reports/site_supplement" "$revision_root/site_returns/UCMC/share" "$revision_root/site_returns/NU/share" "$revision_root/site_returns/RUSH/share"
Rscript diagnostics/revision20261022/stage_additional.R
Rscript diagnostics/revision20261022/finalize.R > "$revision_root/logs/finalize.log" 2>&1
if [[ "${CLIF_SKIP_HISTORICAL_REPORTS:-1}" != "1" ]]; then
 python3 diagnostics/revision20261022/write_report.py
 python3 diagnostics/revision20261022/write_submission_materials.py
fi
Rscript diagnostics/revision20261022/provenance.R
