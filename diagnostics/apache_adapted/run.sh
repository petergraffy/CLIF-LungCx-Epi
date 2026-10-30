#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
root="${CLIF_APACHE_OUTPUT:-output/run_APACHE_adapted_NEW}"
if [[ -e "$root" ]]; then echo 'Choose a fresh CLIF_APACHE_OUTPUT.' >&2; exit 1; fi
export CLIF_APACHE_OUTPUT="$root"
mkdir -p "$root/logs"
Rscript tests/apache_adapted.R > "$root/logs/design_tests.log" 2>&1
IFS=',' read -r -a sites <<< "${CLIF_APACHE_SITES:-UCMC,NU,RUSH}"
for site in "${sites[@]}"; do
 Rscript diagnostics/apache_adapted/site_prepare.R "$site" > "$root/logs/prepare_${site}.log" 2>&1
 Rscript diagnostics/apache_adapted/site_train.R "$site" > "$root/logs/train_${site}.log" 2>&1
done
Rscript diagnostics/apache_adapted/central_train.R > "$root/logs/central.log" 2>&1
for site in "${sites[@]}"; do
 Rscript diagnostics/apache_adapted/site_evaluate.R "$site" > "$root/logs/evaluate_${site}.log" 2>&1
 Rscript diagnostics/apache_adapted/site_publish.R "$site" > "$root/logs/publish_${site}.log" 2>&1
 Rscript diagnostics/apache_adapted/site_audit.R "$site" > "$root/logs/audit_${site}.log" 2>&1
done
Rscript diagnostics/apache_adapted/finalize.R > "$root/logs/finalize.log" 2>&1
