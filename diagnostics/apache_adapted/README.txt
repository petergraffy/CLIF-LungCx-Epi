CLIF-adapted APACHE-II APS+age extension — October30,2026

This is an independently implemented modified score and a new predictive
sensitivity. It is NOT standard APACHE II or APACHE IV. Original72h phenotypes,
assignments, outcome fits and SOFA comparisons remain frozen.

Implementation: code/lib/apache_adapted.R
Protocol: docs/apache_adapted_methods.md
Completed results: output/run_APACHE_adapted_20261030_v2/reports/
Current site packets: same root/site_packets/SITE/share/

Portable component capture is automatic in ./run.sh derive or supplement.
It uses configured source metadata, matches only exact frozen numeric keys,
and keeps component rows/patient identifiers local. APACHE tables are in
share/supplement/apache_adapted_*.csv. No source score is silently replaced.

The reviewer runner diagnostics/revision20261022/run.sh now runs this complete
model extension and includes its allowlisted exports in the site returns.
For standalone reproduction from existing frozen primary caches:
  export CLIF_APACHE_PRIMARY=/path/to/frozen_primary_root
  export CLIF_APACHE_OUTPUT=/path/to/NEW_output
  export CLIF_APACHE_SITES=UCMC,NU,RUSH
  bash diagnostics/apache_adapted/run.sh
The runner refuses existing output roots. If a source folder moved, explicitly
set CLIF_APACHE_TABLES_SITE=/path/to/tables; only exact-match metadata may enrich
the cached numeric observations. Never replace original cached measurements
with a newer release inside a score-only comparison.

Phases:
  site_prepare.R SITE   local first24h scoring, complete/missing/proxy metadata
  site_train.R SITE     original other-site template assignments, local ridge fits
  central_train.R       allowlisted coefficients only; freeze other-site predictors
  site_evaluate.R SITE  same-patient held-out evaluation and200 patient bootstraps
  site_publish.R SITE   safe exact aggregate counts/provenance
  site_audit.R SITE     score/range/case/held-out/prediction invariants
  finalize.R           aggregate-only report/figures and source snapshot

All files in sites/SITE/local/ and source inputs remain local. Only reviewed
sites/SITE/export/ files are exchange candidates. Original frozen artifacts
are hash-checked. Spec incompatibility, duplicate source versions and unsafe
patient-column exports stop processing. Entirely unknown support stays remain
in the all-case cohort; no new primary cohort exclusion is introduced.

Use add_to_return.R SITE OLD_SHARE NEW_SHARE to version this score extension
onto an already completed site packet without overwriting the parent. Parent
consensus must match, site audits must pass, and the collector rejects mixed
APACHE-extension versions or score specifications across sites.

Cohorts:
- complete_hb_proxy: all12 observed/substituted APS components plus age.
- all_case_missing_indicators: all primary cases, partial sum plus12 indicators.
- matched_complete_SOFA: complete adapted score and original modified SOFA;
  both scores compared on the same held-out patients.

Chronic-health points and acute-renal-failure doubling remain unknown/excluded.
Hematology uses estimated Hct=3*measured Hb(g/dL) at every site. Observed-Hct
scores are separate availability diagnostics, not mixed into the primary proxy.
Arterial pH is preferred; absent pH may use bicarbonate as a labeled substitution.
Missing reference units assume canonical CLIF units and are counted explicitly.
Absent GCS is not set to15. Primary incomplete totals remain missing.

Nested comparisons report AUC/Brier, calibration point estimates, paired
increments, calibration bins and exploratory decision curves for death and
composite. Event partitions below10 are withheld. Coefficientwise other-site
synthesis retains flagged zero priors for unsupported terms; required score/
intercept coefficients must be estimable before evaluation. Bootstrap CIs
condition on frozen learners; no claim of independent prospective validation.

Tests:
  Rscript tests/apache_adapted.R
  Rscript tests/site_supplement.R
  Rscript tests/site_return_collection.R
