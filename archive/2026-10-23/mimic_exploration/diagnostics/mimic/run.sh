#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
: "${CLIF_MIMIC_ARCHIVE:?Set CLIF_MIMIC_ARCHIVE to the supplied MIMIC CLIF ZIP}"
root="${CLIF_MIMIC_OUTPUT:-output/run_MIMIC_external_20261023}"
if [[ -e "$root" ]]; then echo "Select a fresh CLIF_MIMIC_OUTPUT; current outputs are preserved." >&2; exit 1; fi
export CLIF_MIMIC_OUTPUT="$root"
mkdir -p "$root/local" "$root/logs"
python3 diagnostics/mimic/prepare_archive.py "$CLIF_MIMIC_ARCHIVE" "$root/local/source_tables"
Rscript diagnostics/mimic/prepare.R > "$root/logs/prepare.log" 2>&1
for site in UCMC NU RUSH; do Rscript diagnostics/mimic/site_train_portable.R "$site" > "$root/logs/train_$site.log" 2>&1; done
Rscript diagnostics/mimic/freeze_portable.R
Rscript diagnostics/mimic/prepare_waterfall.R
CLIF_WATERFALL_OUTPUT="$root" "${CLIF_WATERFALL_PYTHON:-output/run_UCMC_waterfall_20261022/venv/bin/python}" diagnostics/waterfall/run_waterfall.py > "$root/logs/waterfall.log" 2>&1
Rscript diagnostics/mimic/derive.R > "$root/logs/derive.log" 2>&1
Rscript diagnostics/mimic/evaluate.R > "$root/logs/evaluate.log" 2>&1
Rscript diagnostics/mimic/associate.R
Rscript diagnostics/mimic/characterize.R
Rscript diagnostics/mimic/local_discover.R MIMIC > "$root/logs/discovery.log" 2>&1
Rscript diagnostics/mimic/audit.R > "$root/logs/audit.log" 2>&1
Rscript diagnostics/mimic/report.R
python3 diagnostics/mimic/write_report.py
Rscript diagnostics/mimic/provenance.R
