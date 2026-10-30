Portable SOFA completeness capture

Use the normal site workflow; it now captures these diagnostics automatically:
  ./run.sh derive
For an already completed derivation:
  Rscript code/04_supplement.R /path/to/site_config.json /path/to/baseline_run
See docs/site_return_workflow.md for the complete frozen-model/sensitivity return.

The worker stages may also be invoked independently for debugging:
  export CLIF_SOFA_CONFIG=/path/to/site_config.json
  export CLIF_SOFA_BASELINE=/path/to/baseline_run
  export CLIF_SOFA_OUTPUT=/path/to/NEW_local_audit_output
  Rscript diagnostics/sofa_completeness/audit.R SITE
  Rscript diagnostics/sofa_completeness/extra_source_audit.R SITE
  Rscript diagnostics/sofa_completeness/extend.R SITE
SITE must match config.site_name. There are no inferred UCMC/NU/Rush paths.
CSV, parquet and FST source tables are supported; optional source names/LOINCs
are audited when present. Parameters come from the prepared source run.

All exact patient-level scores/evidence stay local. The code/04 wrapper creates
suppressed aggregate share/supplement tables, a required-file inventory and
manifest. Returning the standalone audit's local/ or reports/ is not the
institutional exchange contract. New source/parameters require fresh run IDs.
Frozen source scores, clusters and model fits are not changed by capture.

Measurement recovery covers exact-timestamp unambiguous SBP/DBP calculated MAP,
simultaneous GCS components where total is absent, and explicit arterial PO2 or
serum total bilirubin source labels/LOINCs plus same-result-time arterial panel
FiO2. Conflicting, out-of-range and incomplete BP/GCS tuples are rejected.
No missing domain is assigned normal.

Additional CV review handles clearly named active agents categorized as dextrose,
and dobutamine's any-positive-dose criterion without requiring weight conversion.
Unknown drug semantics outside these explicit mappings are not inferred.
A separate sensitivity uses latest measured weight within24h before ICU,
bounded by hospital admission, only when first24h weight is unavailable.
Another sensitivity estimates PaO2 from measured SpO2 50<=x<97 only when the
respiratory score is missing, paired within1h to raw FiO2 or documented room air.
The equation matches installed CLIFpy2.1.0 utils/sofa_polars.py
_impute_pao2_from_spo2. It is not measured PaO2 or validated standard SOFA.
Reference: https://github.com/Common-Longitudinal-ICU-data-Format/clifpy

First24h is [first ICU arrival,first ICU arrival+24h). Patients belong to the
72h landmark cohort. ABG pairing uses source result timestamps. The pre-ICU lab
check changes the lab window and is availability-only, not a primary score.
Renal scoring is creatinine-only; dialysis codes do not recover urine volume.
Adopting candidate scores requires new model versions and held-out evaluation.

Tests:
  Rscript tests/sofa_recovery.R
  Rscript tests/site_supplement.R

Historical development snapshot:
output/run_sofa_completeness_20261026/ contains the original three-site diagnostic.
finalize.R reproduces only that fixed historical report/figure (including its
version-specific interpretation); it is NOT the portable site finalizer.
Current distributed finalization is code/lib/site_supplement.R.
