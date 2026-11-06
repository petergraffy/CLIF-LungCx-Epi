# Two-stage federated run

Stage 1 learns local clusters and returns aggregate fingerprints. The coordinator
learns the common phenotype templates, freezes all mappings and issues the tracked
`frozen/stage2_v1/` bundle.

Stage 2 is independent of Stage 1 folders. Each admitted site configures its CLIF
source tables and pinned environment, then runs `scripts/site.R stage2` once.
It rebuilds the eligible and landmark cohorts, covariates, bounded respiratory
waterfall and raw-evidence ARF features. It assigns all support-sufficient stays
to the closest frozen prototype and retains insufficient support as I. It runs
no local clustering, silhouette selection or clustering bootstrap. Frozen primary,
sensitivity and leave-one-site-out templates are unchanged.

One aggregate `stage2_return/` contains all local analyses and inputs for central
pooling. It includes a rebuilt-versus-Stage1-reference cohort-count comparison and
source timestamp-schema audit. Current CLIF tables/shared derivation can differ
from the original Stage1 run; counts alone cannot prove patient-level equality.
Source changes during a Stage2 run are rejected using local file-size/mtime checks.
No source rows, patient identifiers or individual predictions are returned.

Prediction uses patient-grouped local cross-validation; leave-one-site-out
phenotype transportability is distinct from externally validated outcome models.
The coordinator pools associations, ROC/tdAUC, and exact event/censor histories for
KM/Aalen–Johansen curves without another site exchange.

Follow [the operational site/coordinator instructions](../code/README.md).
