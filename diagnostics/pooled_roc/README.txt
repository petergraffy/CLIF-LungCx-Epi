Pooled held-out ROC: four models, days 7/14/28

Run locally at each site:
Rscript diagnostics/pooled_roc/site_capture.R SITE PRIMARY_ROOT APACHE_ROOT NEW_SHARE

Return only NEW_SHARE: score_histograms.csv, cohort_counts.csv,
local_validation.csv, manifest.json. No local clinical caches are returned.
The script uses cached matched_complete_SOFA evaluation predictions for clinical
covariates, covariates+modified SOFA, covariates+adapted APACHE and
covariates+phenotype. All four models use identical evaluation stays, and each
site's predictors were frozen from the other sites; there is no local refit.
Both hospital death and hospital death-or-hospice are exported. Days are measured
AFTER the 72-hour landmark. Controls are all noncases, including competing exits
(timeROC AUC2); the evaluated cohorts have no ordinary outcome censoring.

Coordinator, from aggregate returns only:
python3 diagnostics/pooled_roc/central_plot.py NEW_REPORT SITE_SHARE [SITE_SHARE ...]

Exact case/control histogram counts at common public probability bins are summed
across sites before computing the pooled ROC; finished site curves are not
averaged. The frozen 28-day model risks are evaluated at days7/14/28, not refitted
at each horizon. Public bins have width0.00001. AUC and ROC refer to these grouped
scores, not mathematically exact unrounded predictions. The report includes
worst-case unrounded-AUC ambiguity bounds and each site's exact-versus-binned
AUC validation. Tied bins receive half credit. 200 independent within-site
patient-cluster bootstrap count replicates, shared across models/days/endpoints,
provide paired pointwise AUC intervals conditional on the frozen learners.

This is an unadjusted combined-cohort ROC of the site-held-out prediction
mixture; site case mix and between-site calibration differences contribute.
It is not an equal-site standardized estimate or one single jointly fitted
model. Keep site-specific curves as companions. The complete-case restriction
can select a different population from the full primary landmark cohort.
Modified SOFA and CLIF-adapted APACHE are proxies, not standard validated scores.

Current run: output/run_nearest_assignment_20261030/pooled_roc/
Sites: UCMC96, NU198, Rush47; N341. Original models and prior outputs preserved.
Tests: Rscript tests/pooled_roc.R
Reference: https://github.com/cran/timeROC/blob/master/man/timeROC.Rd
