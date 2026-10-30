MIMIC lung-cancer transportability assessment — 2026-10-23

Report: output/run_MIMIC_external_20261023/reports/mimic_lung_cancer_report.txt

Cohorts:
- Current hospitalization must contain primary-organ ICD10 C34.x or ICD9 162.x.
  A metastasis-to-lung or personal-history code alone does not qualify.
- Adult ICU, alive and hospitalized at first ICU arrival +72h.
- Stricter subgroup: same current cancer criterion AND a previous completed
  lung-cancer-coded hospitalization ending by current hospital admission.
- Source POA flags are all null. Neither group is called POA-confirmed.
- Shifted dates are used only for relative timing; all source periods are kept.
  No shifted calendar year, discharge-code Charlson or advanced-cancer variable
  is used as a baseline predictor.
- Missing respiratory support never excludes a hospitalization.

Critical limitation: no room-air device category or device name occurs in this
selected CLIF extract. There are zero room-air hours and zero P1 assignments.
Missing support is not converted to room air. MIMIC cannot validate all four
shared phenotypes or supply P2–P4-versus-P1 association estimates in this state.
The source catalogue, support/ADT audit and separate prediction results are
retained. Do not pool this different-eligibility cohort into the primary model.
Original source oxygen-device documentation or an improved ETL is needed to
recover explicitly documented room-air states.

Predictors use unchanged three-site phenotypes. New portable ridge models are
learned only at UCMC/NU/Rush, excluding year/POA diagnosis covariates; this is not
transport of the old full-covariate equation. MIMIC has no predictor refitting,
recalibration or template learning. Models compare baseline, phenotype, modality
summaries (including hour-72 support), and phenotype added to summaries. Separate
modified-SOFA complete-component models use raw same-time GCS component triplets
only. gcs_eyes maps to gcs_eye; imputed/default GCS totals are not treated as raw.
Prediction intervals use 200 patient bootstrap replicates conditional on frozen
learners. Local OM clustering (100 patient bootstraps) is descriptive only and
occurs after frozen external evaluation; it does not update the shared templates.

Fresh-run reproduction from this repository:
  CLIF_MIMIC_ARCHIVE=/path/to/mimic-clif.zip \
  CLIF_MIMIC_OUTPUT=output/run_MIMIC_external_NEW \
    bash diagnostics/mimic/run.sh
Requires the original three-site frozen rehearsal and the CLIFpy 2.1.0 isolated
Python environment; override CLIF_WATERFALL_PYTHON if necessary. The supplied
source ZIP is read only. Selected Parquet tables are extracted under ignored
output/.../local/source_tables, preserving raw values. The eligibility adapter
uses a marked working copy solely to select cancer-coded encounters; observed
POA fields remain unchanged and diagnosis covariates are not fabricated.

Audit verifies actual current cancer codes, genuinely earlier completed cancer
encounters, all-case retention, no terminal exposure states, unchanged primary
artifacts, independent prediction reproduction and no identifiers in aggregates.

local/ contains patient data, raw source tables and predictions; never exchange
it. aggregate/ contains local aggregate results; exchange/training contains only
three-site coefficient returns. No real institutional exchange or disclosure
approval is implied. All source-data governance requirements continue to apply.
