Federated pooled Kaplan-Meier curves

Verified current run: output/run_pooled_survival_20261030_v2/
Main plots and aggregate tables: reports/
Primary phenotype consensus remains6245d53b0469ab7d1a5b; no clinical model refits.

The pooled product-limit estimator is computed from SUM(events)/SUM(at risk)
at each exact common follow-up duration. This is a combined hospitalization-
cohort curve, with larger sites contributing more observations. It is not an
unweighted or sample-size-weighted average of the finished site KM curves.
Do not call it an adjusted or equal-site-standardized survival curve. Retain
site-specific curves and the original adjusted random-effects HR meta-analysis.

Clocks and endpoints:
- All4533 eligible hospitalizations at first ICU admission, without requiring
  survival through72h. Unknown terminal dispositions are censored at observed
  discharge as before; possible missing-outcome bias remains explicit.
- All3310 landmark hospitalizations from ICU+72h, by the frozen shared P1–P4
  phenotypes plus I. No future72h phenotype is assigned at admission.
- Death and death-or-hospice are separate endpoints. Other hospital exits are
  censored in cause-specific KM. These curves are not marginal event incidence
  and do not establish all-cause survival after discharge; retain the AJ plots.

Two local passes and aggregate-only central steps:
1. Each site runs site_times.R SITE PRIMARY_ROOT NEW_TIME_SHARE. Return the two
   aggregate CSVs plus manifest. They contain counts at relative event durations
   and group/endpoint cohort totals, not identifiers or calendar timestamps.
2. Coordinator runs central_grid.R NEW_GRID SITE_TIME_SHARE [...]. Send that
   shared duration grid and manifest to each site. Integer microsecond ticks
   prevent floating-point CSV round trips from losing event/censor ties; this
   preserves the source timestamp precision rather than daily/hourly binning.
3. Each site runs site_risksets.R SITE PRIMARY_ROOT CENTRAL_GRID NEW_RISK_SHARE.
   Return risksets.csv, bootstrap_risksets.csv and manifest only. Each count is
   exact, including singleton cells. Individual admissions/patient identities
   remain local. Local reconstruction must match ordinary local KM at the grid.
4. Coordinator runs central_pool.R NEW_REPORT CENTRAL_GRID SITE_RISK_SHARE [...].
   It checks site/grid/file hashes, sums counts, computes pooled KM and intervals.
   It opens only aggregate returns. plot.py renders these aggregate results.

Intervals use200 independent within-site patient-cluster bootstrap replicates.
Every admission for a sampled patient receives the same weight; the weights are
shared across that site's curves and endpoints. The coordinator combines each
replicate's risk/event counts before estimating its pooled curve. Intervals are
pointwise percentile intervals, conditional on these three sites and the frozen
phenotype definitions. They do not include phenotype-learning uncertainty, site
sampling uncertainty for a new institution, or dependence from the same patient
appearing at different institutions (no cross-site identity linkage exists).
Hospitalization counts, not globally deduplicated patient counts, are pooled.

Local three-site rehearsal runner:
CLIF_POOLED_SURVIVAL_OUTPUT=/path/to/fresh/run bash diagnostics/pooled_survival/run.sh
Set CLIF_POOLED_PRIMARY and CLIF_POOLED_SITES if needed. In a real federation,
steps1 and3 run on the institutions' machines; the runner here rehearses them
locally. Sharing the two-phase aggregate returns requires no patient identifiers.

Tests: Rscript tests/pooled_survival.R
They compare exact federated counts to centralized KM on synthetic data with
unequal site sizes, tied events/censoring, patient resampling weights, empty
sites, and decimal-CSV round-trip errors. All real site reconstructions also
passed. Original fitted clinical models and source artifacts remain unchanged.

Methods reference:
https://www.stat.ethz.ch/R-manual/R-devel/library/survival/html/survfit.html

October30 discharge-panel update
The current nearest-template analysis now exports live non-hospice discharge
incidence using Aalen-Johansen, with death and hospice as competing outcomes.
This replaces the death-or-hospice-free KM panel in the main pooled figures.
The former composite KM remains in aggregate tables for supplementary use.
All endpoint grids use the union of exit times for each cohort/phenotype so
competing events update the shared survival probability at their exact times.
Both local conventional AJ checks and aggregate synthetic tests must pass.
200 patient-cluster bootstrap replicates are combined at the count level.
Unknown dispositions are censored at observed exit; no postdischarge events
are inferred. Discharge-alive here explicitly excludes hospice.

Final main return packaging:
Rscript diagnostics/pooled_survival/add_to_return.R SITE TIME_SHARE RISK_SHARE OLD_RETURN NEW_RETURN
This preserves the earlier packet, adds both aggregate phases, updates the
inventory and hashes, and sets pooled_survival_discharge_complete=true.
Use CLIF_REQUIRE_POOLED_SURVIVAL=true for final site_return/collect.R collection.
Current complete packets: output/run_nearest_assignment_20261030/site_packets_discharge/SITE/share/.

Paired competing incidence curves
--------------------------------
central_pool.R also returns pooled_competing_incidence.csv and
site_competing_incidence.csv for death_or_hospice and live_discharge.
These are the two absorbing causes: unknown dispositions are censored at
observed exit. Both CIFs use all known terminal exits in the common
Aalen-Johansen survival product. Their sum plus the probability of remaining
in hospital without either event is checked to equal one at every time.
The same existing aggregate risksets and patient-cluster bootstrap counts
support both curves; no additional patient-level exports are required.
Plot the landmark phenotype curves with:
  python3 diagnostics/pooled_survival/plot_competing.py PATH_TO_POOLED_REPORTS
The figure omits embedded captions/subtitles, retains all phenotype groups,
and reports pointwise bootstrap uncertainty and numbers at risk. Curves are
unadjusted and describe the pooled observed cohort, not equal site averages.
