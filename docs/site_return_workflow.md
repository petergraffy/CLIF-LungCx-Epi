# Site returns

Stage 1 returns `stage1_return/` and waits for the frozen release. Stage 2 is one
local run (`scripts/site.R stage2`) and one upload (`stage2_return/`). Follow
[the complete site instructions](../code/README.md).

Stage 2 rebuilds its cohort, clinical covariates, bounded waterfall and hourly
features from the configured CLIF tables. Stage 1 folders are not needed. No local
clustering or prototype learning is performed. Cohort-count comparisons and an
input-rebuild audit accompany the return.

The frozen bundle lives in `frozen/stage2_v1/` and is part of the repository.
The single Stage 2 packet contains associations and diagnostics, all sensitivity
models, descriptors and supplement audits, exact trajectory and Sankey counts,
phenotype transport agreement, locally cross-validated prediction histograms,
and exact event/censor histories for centrally pooled KM/Aalen–Johansen curves.

Patient identifiers, clinical rows, individual predictions and RDS caches remain
local. Local cross-validation keeps each patient's admissions together. This is
internal prediction validation; leave-one-site-out phenotype mapping is distinct
from external validation of an outcome predictor. No second site pass is needed.
