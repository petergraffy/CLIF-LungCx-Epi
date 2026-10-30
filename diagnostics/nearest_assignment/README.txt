Current classification policy: always nearest for support-sufficient stays

Every site uses fed_assign in diagnostics/federated/contract.R. No hospitalization
can become unassigned because of a small distance margin or a large distance.
The closest weighted Euclidean template is always selected; exact ties use the
fixed template order. Low-margin/large-distance booleans remain diagnostic only
and never enter the clinical adjustment or prediction models. Support below the
existing information threshold remains I; its ETL repair is a separate task.

The primary threshold is12 known-support hours of72; the earlier-landmark
sensitivities retain their existing proportionally scaled information thresholds.
No missing support is imputed as room air and no cohort member is removed.
The same rule applies to all full/held-out and sensitivity/absorbing mappings.

Current run: output/run_nearest_assignment_20261030/
- primary/: unchanged learned local clusters/shared templates; fresh frozen
  assignment contracts; refit associations, predictors and held-out evaluation.
- revision/: all8 inferential sensitivities, descriptive terminal-state mapping,
  ARF/code/cancer sensitivities, modified-SOFA comparisons and audits rerun.
- apache_adapted/: same cached component values; classification-dependent models
  and transported comparisons rerun, without rereading clinical source values.
- exact/: fresh Table1/Table2, source/count/profile returns, PH/piecewise checks,
  KM/AJ curves, tdAUC/ROC and manuscript figure layouts. Exact counts, no masking.
- pooled_survival/: pooled KM from exact risk/event counts,200 within-site
  patient-cluster bootstrap replicates and exact risk tables.
- site_packets/SITE/share/: complete refreshed aggregate return for each site.
- reports/: independent assignment audit, transitions and full-return collection.

Full primary:54 formerly rejected stays now assigned (UCMC10, NU34, Rush10).
All3310 hospitalizations retained;3180 classified;130 remain support-insufficient
(NU4, Rush126). New pooled P sizes:1059 /1250 /379 /492; I130.
The remaining information gaps are not the ambiguity/poor-match issue.

Reproduction: prepare.R versions the copied parent templates/local clusters and
records their parent manifests/hashes. refit.sh runs fresh site outcome passes
and other-site predictors against these frozen templates. review.sh/apache.sh
reuse cached source components and rerun phenotype-dependent models. ancillary.sh
captures the exact supplementary outputs. audit.R independently recomputes all
nearest-template distances across108 site/version passes and verifies templates,
cohorts, raw outcomes and already-classified stays are unchanged. package.py
assembles only current model/descriptor exports, reusing classification-independent
source diagnostics for the identical baseline cohorts. Original results are kept.

For future work, source diagnostics/nearest_assignment/environment.sh from the
repository root. config/active_analysis.json records the active output paths.
The original support/ARF preprocessing, clustering weights, shared K and score
values are not changed by this assignment-only update. NU's source catalogue
audit uses its resolved current path; cached clinical values remain frozen.

This is a transparent posthoc policy revision after initial outputs were examined;
it is not new untouched validation. Outcomes were not used to choose assignments.
Uncertainty continues to condition on the frozen phenotype learning/classifier.
Tests: Rscript tests/nearest_assignment.R; frozen-artifact/held-out audits also run.
