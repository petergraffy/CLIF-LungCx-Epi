Federated clustering/outcome rehearsal -- 2026-10-22

Completed output:
  output/run_federated_rehearsal_20261022/reports/
Read federated_rehearsal_report.txt for results and exact methods/limitations.
This is a logically separated rehearsal on one machine, not OS isolation,
cryptographic secure aggregation, or differential privacy.

Reproduce with a FRESH destination from the repository root:
  CLIF_FEDERATED_OUTPUT=output/run_federated_rehearsal_NEW \
    bash diagnostics/federated/run_rehearsal.sh
The runner stops rather than overwriting an existing rehearsal. Scientific
settings and seeds are defined in contract.R; model settings in model_lib.R.
Changing output paths/source hashes changes artifact IDs, not the intended
numerical analysis. This script uses the three already-validated local
waterfall caches; it does not repeat source ETL or waterfall processing.

Default local inputs, one per site worker:
  output/run_UCMC_waterfall_20261022
  output/run_NU_waterfall_20261022
  output/run_RUSH_waterfall_20261022
Override separately with CLIF_FEDERATED_INPUT_UCMC, _NU, or _RUSH. Each input
requires local/input_manifest.json, local/waterfall_manifest.json and
local/waterfall_derivation.rds with h72_waterfall_support_with_raw_arf.
The site field and timestamp/cohort/sequence invariants are checked locally.
This rehearsal targets the current three-site cohort; final expected cohort
size and interpretation checks intentionally stop on unexpected changes.

Directory/role boundary:
  sites/SITE/local/         patient records, sequences, fits, predictions
  sites/SITE/export/        allowlisted aggregate returns only
  exchange/                staged CSV/JSON returns; no patient records/RDS
  coordinator/             templates, frozen predictors, meta-analysis
  reports/                 aggregate tables, six figures, audits, report
Do not send local/ or legacy waterfall cache RDS files to a coordinator.

Passes:
  site_discover.R SITE       local common-cost OM/Ward fit; local k=2..8;
                            100 patient bootstraps; banded fingerprints
  stage_discovery.R         allowlist discovery exports and reject ID fields
  central_discover.R HOLD    learn full or two-site templates, choose global
                            k, freeze definitions; no outcome inputs
  site_outcomes.R SITE VER   frozen per-patient assignment; local Cox and
                            ridge predictive models; suppressed aggregates
  stage_outcomes.R          validate and stage local model packages
  central_pool.R VER        meta-analyze training-site associations and
                            freeze coefficient-synthesized predictors
  site_transport.R SITE     evaluate the two-site predictor at held-out site
  site_comparators.R SITE    exploratory harder respiratory-summary models
  central_comparators.R      synthesize training-site comparator coefficients
  site_transport_extended.R evaluate all four predictors and paired gains
  site_risks.R SITE         local 28-day absolute risks, patient-clustered SE
  site_quality.R SITE       verify all cases retained and reproduce transport
  central_finalize.R        aggregate-only tables/figures and export audit
  write_report.py           plain-text report from central aggregates only
VER is full, loo_UCMC, loo_NU or loo_RUSH. HOLD is none, UCMC, NU or RUSH.

Scope: the 72-hour primary workflow. The previously completed local 24/48h,
absorbing-event, ARF-onset and alternative-physiology sensitivities have not
all been refederated under the new prototype classifier.

Important design choices:
- Original hourly support/ARF inputs remain joint-state sequences. Waterfall
  changes support; ARF states remain based on the original observed evidence.
- This new rehearsal uses common constant substitution=2 and indel=1, rather
  than the prior site's learned TRATE matrix. Global k is not forced to four.
- Six-hour modal fingerprint intervals span at least ten possible counts;
  cluster sizes and stability are banded. The coordinator reconstructs an
  interval-feasible representative profile, not observed patient counts.
  Small template components can reflect interval uncertainty.
- Full hourly sequences are used for local OM; six-hour summaries feed the
  coordinator and frozen classifier. This is not the exact pooled OM/Ward
  solution. Finer privacy-compatible summaries remain a future sensitivity.
- Sites are equally weighted in template learning. Local clusters can map
  many-to-one. Final patient classification uses frozen template distances.
- Below 12 known-support hours, distance >.75 or margin <.01 produces I.
  All cases are retained, including I; I is not a clinical phenotype.
- Cox models censor follow-up at day 28 and use separate baseline hazards
  for classified versus I cases. Baseline age/sex/race/year/comorbidity,
  advanced-cancer proxy and admission code status are adjusted in core fits.
  SOFA is a labeled complete-case sensitivity, not missing-as-normal primary.
- Model estimates for rare event partitions are withheld; unavailable
  coefficients are not assigned zero for association meta-analysis.
- Prediction is exploratory ridge logistic synthesis, including a flagged
  zero prior for unsupported predictive terms. Held-out patients/outcomes
  never contribute to profile selection or predictor fitting.
- The harder summary comparator was added after the initial baseline result;
  it is explicitly exploratory and does not modify the frozen phenotypes.
- Patient bootstrap intervals condition on frozen training; they do not
  propagate every stage of phenotype discovery uncertainty.
- Data were examined previously. The validation is a development rehearsal,
  not independent prospective confirmation or a causal-effect analysis.

Tests:
  Rscript tests/federated_rehearsal.R
  Rscript tests/federated_immutability.R
The first tests interval utility, assignment and privacy boundaries, local
prediction estimation and held-out input exclusion. The second recomputes
frozen content hashes, validates all outcome-version hashes and confirms
that adding comparators leaves original held-out estimates unchanged.
In the completed run, provenance/discovery_source and
provenance/completed_analysis_source preserve phase-specific code. An
assignment-only matrix-shape issue was repaired before successful outcome
passes; path-override and figure-layout changes were added after analysis.
All scientific definitions and completed frozen artifacts were preserved.

No real between-institution exchange has occurred. External sharing requires
local disclosure review of the complete linked aggregate package. These
safeguards are not a formal differential privacy guarantee.
