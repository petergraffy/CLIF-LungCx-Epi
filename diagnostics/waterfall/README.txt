CLIF waterfall comparison

This experiment reads the configured completed baseline derivation and reuses
its exact cohorts, time origins, terminal events, clinical covariates and
clustering settings. It does not overwrite the original derivation or switch
the production primary analysis automatically.

Tested environment: official clifpy 2.1.0, pandas 2.3.3, Python 3.12. The
installed waterfall function's SHA256 and package versions are recorded in
local/waterfall_manifest.json. The exact Python environment is recorded in
local/python_requirements.txt. The imported package falls back to DuckDB for
an unrelated Polars runtime check; the invoked waterfall itself uses pandas
and DuckDB. No upstream waterfall algorithm code is modified.

Run from repository root with a fresh output directory:
  CLIF_RUN_CONFIG=config/config.json \
  CLIF_WATERFALL_OUTPUT=output/run_UCMC_waterfall_NEW \
  CLIF_WATERFALL_PYTHON=output/run_UCMC_waterfall_20261022/venv/bin/python \
  bash diagnostics/02_compare_waterfall.sh

The baseline private config and prepared/derivation RDS files must be present.
The baseline source fingerprint is checked. All inputs, processed records,
patient-level assignments and review queues remain in ignored local output.
Aggregate files are exploratory local outputs, not disclosure-reviewed
federated share packages.

Processing boundaries
For each cohort/window separately, input records are limited to the same
hospitalization, a six-hour pre-origin lookback, and the earlier of the
completed landmark or terminal event. Empty boundary markers extend the
waterfall's hourly scaffold across the requested window; they do not provide
any device values. Actual boundary-time source records are never duplicated.
Output is trimmed again to the exact half-open interval. Hospitalization,
not patient, is the waterfall grouping key. Absorbing states are added only
after processing and remain confined to the descriptive absorbing sensitivity.

The official function runs with bfill=False. This disables numeric backward
fill; standard device/mode heuristics may still look ahead within the already
completed window. Device forward fill has no six-hour expiry in the standard
waterfall; it can extend to the next documented change or window boundary.
Raw FiO2 is normalized and bounded using the baseline row-wise rule before
waterfall. Inferred/filled rows are distinguished from matching raw device
observations, including synthetic scaffold rows.

Comparison arms
The nine existing analyses are rerun with waterfall device and numeric fields.
The extra h72_waterfall_support_with_raw_arf arm uses waterfall support while
holding original hourly ARF states fixed. This isolates device filling from
the additional influence of filled FiO2 on ARF detection. First-detected-ARF
anchors remain based on the original observed evidence for a paired cohort
comparison; onset is not redefined using imputed FiO2.

Baseline ICU severity, diagnosis features and code status are held fixed.
Their missingness is not solved by this experiment. Adjusted association
models and competing-risk curves are run locally for the full waterfall
72-hour cohort. Predictive validation is not repeated: waterfall would need
to be included in the fold-isolated learning procedure first. Cluster labels
are provisional local identifiers and use no outcomes in their derivation.

Quality checks
test_boundaries.py verifies landmark and terminal truncation, preservation
of actual source values at the start boundary, no numeric backward fill,
no inference of a specific device from oxygen flow alone, and retention of
unobserved encounters. Real-data checks verify one assignment per cohort
hospitalization, expected sequence lengths, unchanged primary event totals,
before/after entirely unknown encounter counts and correct absorbing terminal bins.

qa.R also records source/processed conflicts, including positive source flow
with waterfall room air and filled FiO2 0.21 on low-flow support. These are
review flags, not confirmed source errors. Record-level queues stay local.
The source input should not be declared defective solely from pre-waterfall
missingness or from the decrease achieved by imputation.
