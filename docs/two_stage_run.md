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

Follow [the site instructions](../code/README.md).

## Coordinator workflow

Inspect the original Stage 1 returns, choose supported phenotype counts and record
explicit legacy-source acceptance in a private coordinator config. All sites
participate in Stage 2. Any missing discovery contribution must be named and
justified in `discovery_exceptions`; it does not remove that site from modeling.

```sh
Rscript --vanilla scripts/coordinator.R inspect config/coordinator.json
Rscript --vanilla scripts/coordinator.R freeze config/coordinator.json
```

`freeze` learns outcome-free full and leave-one-site-out templates for every
variant, copies the fixed analysis plan, and writes a checksummed
`stage2_parameters/` bundle. Publish those aggregate CSV/JSON templates under
`frozen/stage2_v1/` in Git; source tables, institutional uploads, private configs
and patient caches remain outside Git. The shipped release records exactly which
sites contributed to each discovery, including explicitly accepted omissions.

Sites run `stage2` once and upload `stage2_return/`. After every participating site
returns, collect all packets with the same private coordinator config:

```sh
Rscript --vanilla scripts/coordinator.R collect config/coordinator.json /path/SITE_A/stage2_return /path/SITE_B/stage2_return /path/SITE_C/stage2_return
```

If any site used the separate Python fallback, collect the complete roster
with the compatible collector instead (continue the command with every site):

```sh
Rscript --vanilla python_fallback/collect.R frozen/stage2_v1 output/run_pooled_stage2_v1 .venv/bin/python /path/SITE_A/stage2_return /path/SITE_B/stage2_return /path/SITE_C/stage2_return
```

It validates native and fallback packets against their distinct issued
contracts and the same frozen parent phenotypes, then uses the existing
pooling methods. Received manifests are preserved. See the
[fallback instructions](../python_fallback/README.md).

The coordinator checks the release identity and source hashes, pools association
estimates and prediction histograms, reconstructs pooled survival/competing-risk
curves, and retains all site tables for the supplement. No outcome-predictor
bundle is sent back to sites.

A late site completes the retained Stage 1 workflow first. Create an admission revision without relearning the full prototypes:

```sh
Rscript --vanilla scripts/admit_stage2_site.R frozen/stage2_v1 /path/LATE_SITE/stage1_return frozen/stage2_v1_admission2 "Reviewed Stage1 discovery; admit site to independent Stage2 reconstruction."
```

Issue that revision to all sites before their runs. The late site uses the full frozen templates for its
external phenotype transport assessment because it contributed no discovery data.
