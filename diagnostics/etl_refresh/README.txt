UCMC respiratory ETL refresh audit — 2026-10-26

Default output: output/run_UCMC_etl_refresh_audit_20261026/etl_refresh_report.txt
Newer source: Box-Box/03-CLIF-2.1/2.1.0, not the unchanged configured Desktop copy.

This audit freezes the previous 650 primary landmark hospitalizations and all
nonrespiratory clinical inputs. It compares the saved respiratory source with
updated Box respiratory rows, rebuilds the 72h raw panel, then repeats the same
bounded CLIFpy waterfall while preserving refreshed raw-evidence ARF states.
Anchor tables are independently checked for cohort/outcome/timing changes.
No clustering, prognostic refitting or pooled-result replacement occurs.

From the repository root, choose a FRESH audit destination:
  export CLIF_UCMC_REFRESH_PATH=/path/to/updated/UCMC/2.1.0
  export CLIF_UCMC_REFRESH_OUTPUT=output/run_UCMC_etl_refresh_NEW
  Rscript diagnostics/etl_refresh/check_ucmc.R
  CLIF_WATERFALL_OUTPUT="$CLIF_UCMC_REFRESH_OUTPUT" \
    output/run_UCMC_waterfall_20261022/venv/bin/python \
    diagnostics/waterfall/run_waterfall.py
  Rscript diagnostics/etl_refresh/finish_ucmc.R
  Rscript diagnostics/etl_refresh/compare_anchors.R
  python3 diagnostics/etl_refresh/write_report.py

Requires the original preserved UCMC prepared/derived caches and the pinned
CLIFpy 2.1.0 environment. Raw-source rows and patient-level panels stay in ignored
local/. Aggregates are local diagnostics, not disclosure-approved institutional
returns. No existing primary artifacts are overwritten. Updated source-table
labels may reflect producing-ETL inference; raw EHR documentation provenance is
not established merely by improved CLIF device coverage.
