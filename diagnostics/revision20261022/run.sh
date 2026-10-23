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
Rscript diagnostics/revision20261022/stage_additional.R
Rscript diagnostics/revision20261022/finalize.R > "$revision_root/logs/finalize.log" 2>&1
python3 diagnostics/revision20261022/write_report.py
python3 diagnostics/revision20261022/write_submission_materials.py
Rscript diagnostics/revision20261022/provenance.R
