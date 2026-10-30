Main-paper figure layout examples

Run from the repository: python3 diagnostics/manuscript_figures/examples.py
Requires numpy and matplotlib. Reads aggregate October26 primary/reviewer
results, October30 adapted-APACHE comparisons and publication descriptors.
Writes PNG and vector PDF examples with source hashes to:
output/run_manuscript_figure_examples_20261030/

These are proposed layouts, not final submission figures. No models are refit.
Figure1: cohort and federated workflow; detailed exclusion reasons not shown.
Figure2: shared template support and ARF marginals, plus observed site membership.
Figure3: released absolute risks and adjusted site/pooled associations.
Figure4: paired held-out AUC increments in explicitly separate evaluation cohorts.
Exact sufficient-statistic counts and released count bands are retained as given.
Supplemental death/discharge state profiles remain in the revision-gap outputs.

October30 refinement:
Rscript diagnostics/manuscript_figures/site_curves.R SITE PRIMARY_ROOT NEW_SITE_OUTPUT
python3 diagnostics/manuscript_figures/refinements.py

Current output: output/run_manuscript_figures_20261030_v2/
The plotting script opens aggregate returns only. One-site workers open only
that site's frozen cache and return daily AJ probabilities/SE/risk-count bands.
Released original checkpoints are checked for exact agreement. Early/rare
partitions are jointly masked as before, leaving visible gaps. These are daily
checkpoints, not exported patient event times, and are never smoothed.

Template smoothing is display-only common convex cubic interpolation between
six-hour centers. It preserves the original values at all centers, bounds and
sum-to-one; it does not re-estimate hourly frequencies or change any assignments.

Absolute AUROC plots compare models in identical cases within a panel/site.
Separate SOFA and APACHE panels use different cohorts. The matched-score plot
uses identical complete-SOFA/APACHE cases for direct comparison. Paired delta
AUROC intervals remain necessary for inference; marginal interval overlap does
not test a model difference. The earlier delta plots are preserved.

Future complete site returns automatically include the daily curve CSV through
site_return/capture.R; collect.R rejects mixed curve-return phases. Frozen
packets can be extended with add_curves_to_return.R without changing parents.
Tests: Rscript tests/competing_event_curves.R

Survival and time-dependent discrimination (October30 correction)

Rscript diagnostics/manuscript_figures/site_survival_tdauc.R SITE PRIMARY_ROOT NEW_SITE_OUTPUT
Rscript diagnostics/manuscript_figures/add_survival_to_return.R SITE SURVIVAL_SHARE OLD_SHARE NEW_SHARE
python3 diagnostics/manuscript_figures/plot_survival_tdauc.py

Current verified calculations: output/run_survival_tdauc_20261030_v2/
UCMC/NU successful local exports were retained from the first calculation;
Rush was rerun to handle empty event-free risk sets after all hospital exits.
Original model fits, predictions and prior figures are preserved.

Survival: all4533 eligible hospitalizations enter at first ICU arrival. Death is
an event for mortality; death/hospice is an event for the composite. Other exits
are censored, not counted as deaths. Unknown terminal dispositions (63 overall)
are censored at observed discharge, which assumes no known target event before
that time and cannot eliminate missing-outcome bias. These cause-specific KM
functions are not marginal hospital-event probabilities and not all-cause
post-discharge survival. Keep the competing-event AJ curves as companions.
No72h phenotype is retrospectively assigned to people at admission: primary
phenotype-specific survival and discrimination start at the72h landmark.

Time-dependent AUC: timeROC cumulative/dynamic IPCW definitions. Controls2
(main) include competing terminal exits and patients still event-free; controls1
(supplement) include only patients free of any terminal hospital event at t.
All primary evaluation hospital events are observed with no ordinary censoring
before the evaluation horizon, so IPCW weights equal1. Once every hospital exit
has occurred, controls1 are empty (AUC1 unavailable); competing-exit controls2
can still identify AUC2. The code checks this distinction explicitly.
Frozen28-day prediction scores are evaluated at horizons1–28; no horizon-specific
or held-out-site model is fit. Pointwise95% intervals use200 patient-cluster
bootstrap samples shared across models/days/control definitions. Bootstrap CIs
condition on the frozen learner; they are not simultaneous confidence bands.
At28d, primary AUC2 exactly reproduces the earlier binary-horizon AUROC.
Time-dependent ROC curves at7/14/28d use a common fixed FPR grid; individual risk
scores, score cutoffs and patient ROC vertices are never returned. Exact AUC
values come from timeROC, not integration of the sampled display.

The new returns include survival_curves.csv, survival_denominators.csv, tdauc.csv,
tdroc.csv and a methods/hash manifest under federated/survival_tdauc/. All site
phases must match for coordinator collection. Existing full/sensitivity models
remain the original primary analyses; these curves are exploratory extensions.
Tests: Rscript tests/survival_tdauc.R and tests/site_return_collection.R.

Primary methods references:
https://github.com/cran/timeROC/blob/master/man/timeROC.Rd
https://www.stat.ethz.ch/CRAN/web/packages/survivalVignettes/vignettes/tutorial.html
