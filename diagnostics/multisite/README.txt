Three-site respiratory-support/ARF comparison (2026-10-22)

Run from the repository root. Site configs and patient-level results remain
inside ignored output/run_*/local directories. The existing UCMC config,
source code fingerprint and completed UCMC runs are unchanged.

1. Create site configs from config/config.json, retaining every scientific
   parameter. Change only site_name, tables_path and analysis.run_id.
2. Rscript diagnostics/multisite/run_site.R /absolute/path/to/site_config.json
3. CLIF_RUN_CONFIG=/absolute/path/to/site_config.json \
   CLIF_WATERFALL_OUTPUT=output/run_SITE_waterfall_20261022 \
   CLIF_WATERFALL_PYTHON=output/run_UCMC_waterfall_20261022/venv/bin/python \
   bash diagnostics/02_compare_waterfall.sh
4. After all sites finish:
   Rscript diagnostics/multisite/compare_sites.R

Each site runs all nine baseline analyses, all nine full-waterfall analyses,
and an additional joint-state arm with waterfall support and original ARF
states. First detected ARF anchors stay fixed to observed baseline evidence.
100 patient bootstrap replicates assess both 72-hour waterfall joint arms.

Comparisons use each site's locally derived clusters. Equal cluster numbers
DO NOT denote matched phenotypes across sites. No pooled effect estimate or
external prediction validation is claimed. Patient-level source data remain
local; comparison aggregates are exploratory and have not passed disclosure
review for external sharing. Calendar-year and cohort-size differences,
missing physiology/severity/code status, and filled FiO2/device conflicts
must be assessed before interpreting differences as clinical effects.

Rush loader compatibility
Rush code_status has hospitalization_id as well as patient_id. The existing
loader projects out the former, then mistakenly tries filtering it. Instead
of changing the frozen analytical source or original source files, a local
schema view drops that unused column; all other tables symlink to originals.
The reader then correctly uses its existing patient-ID filter. Clinical code
status features already use patient-ID and admission-bounded time windows.
Run prepare_code_status_view.py with a private config path to reproduce this
projection. local/schema_adapter.json records unchanged source/output rows;
local/source_config.json preserves the original data location. This is a
loader compatibility workaround, not a change to code-status categories.

Date-boundary audit
The frozen utc() helper parses date-only config strings as NA. Use full
UTC timestamps (2018-01-01 00:00:00 and 2025-01-01 00:00:00) in new configs.
UCMC and NU's complete included cohorts were checked: every admission is
already within the intended window, so their inclusion sets are unaffected.
Rush v3b exposed out-of-window admissions and was aborted. Only Rush v3c,
with explicit timestamp boundaries, enters the comparison. The final
comparison asserts exact admission timestamps are in the intended interval
for every site's full eligible cohort. After completing and archiving these calculations, core utc() was fixed to
parse date-only strings; read_run_config now rejects invalid/reversed bounds.
The core projected-ID filter was also fixed, eliminating the need for the
Rush schema workaround on future runs. tests/input_contracts.R covers both.
Completed analyses preserve their original source fingerprints. A frozen
code/ and diagnostics/ snapshot is stored in the comparison local/ folder as
analytical_source_before_input_fixes.tar.gz. To reproduce these exact cached
runs use that source snapshot; future runs with fixed core require fresh run
IDs under the source/parameter cache contract.
