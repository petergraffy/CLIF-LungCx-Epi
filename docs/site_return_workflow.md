# Site supplement and return workflow

**Current aggregate policy (October30): return exact aggregate counts, including
small cells and exact numbers at risk. Count banding and disclosure masking are
disabled by default. Identifiers, individual sequences/predictions and raw event
timestamps remain local. Historical banded packets are preserved separately.**


The current study uses bounded waterfall support with source-evidence ARF and
frozen interval-template classification. Do not mix this with the earlier raw
baseline one-to-one consensus. The active methods are in `analysis_protocol.md`.
The three-site rehearsal runners are not an installer that automatically runs
a new ten-site federation; the coordinator must distribute the same frozen
contracts and phase inputs to every participating site.

## Automatic local capture

`./run.sh derive` (or `run.bat derive`) now runs source preparation/baseline
derivation and then captures the site supplement. It writes:

- `output/run_SITE_RUN_ID/share/supplement/`: aggregate return tables, inventory
  and source/method manifest; apply local disclosure review before return.
- `local/supplement_exact/`: exact diagnostic counts; never return.
- `local/supplement_sofa/`: source checks, candidate scores, logs; never return.

Existing completed derivations can capture the supplement separately:

```sh
./run.sh supplement
# Or use an explicit site config and completed baseline cache:
Rscript code/04_supplement.R /path/to/site_config.json /path/to/baseline_run
```

Set `CLIF_SUPPLEMENT_OUTPUT` to a fresh directory for a read-only diagnostic
extension of a frozen run. Completed supplement folders are never overwritten.
Capture failure stops before a completion manifest; inspect local logs. New core
code changes the source fingerprint, so new derivations require fresh run IDs.
The explicit supplement command may read a historical cache without relabeling
its original source/method version.

Copy `config/site_metadata_template.csv` to `config/site_metadata.csv` and fill
site statements/provenance, or set `CLIF_SITE_METADATA` to a separate CSV with
`field,value` columns. Blank fields remain `site_confirmation_required`.
Do not place patient identifiers or patient-level clinical narratives in this form.
Fields include IRB board, number/date/title, consent/waiver, source release/ETL,
timezone/provenance, DNR/DNI/comfort mapping, oncology timing and disclosures.

SOFA capture uses each site's configured source tables and prepared cache, not
hardcoded UCMC/NU/Rush paths. CSV, parquet and FST source tables are supported.
Missing optional source labels produce explicit unverified coverage rather than
invented mappings. Calculated MAP and recovered source lab data are diagnostic
candidate scores. Recent pre-ICU weight and SpO2-derived PaO2 remain separate
sensitivities. This step does not refit outcomes or replace frozen scores.

## Complete federated return after all frozen site passes

The current `diagnostics/revision20261022/run.sh` now automatically captures each
site's complete return after local models, validation and sensitivity audits.
For an already completed run, execute on that site's machine:

```sh
Rscript diagnostics/site_return/capture.R SITE PRIMARY_ROOT REVISION_ROOT NEW_RETURN_ROOT
```

`PRIMARY_ROOT/sites/SITE/local/` and the corresponding reviewer/sensitivity
site-local caches must be available locally. This worker opens only that site's
clinical caches. Other-site models/templates must already be supplied as the
frozen phase inputs required by the original analysis. No capture step contacts
another institution or sends files.

The only folder to return is **`NEW_RETURN_ROOT/share/`**. It contains:

| Folder | Contents |
|---|---|
| `supplement/` | SOFA completeness/recovery/source audits, missingness patterns, cohort/early-event reconciliation, oncology/code-status availability, annual evidence, reporting form |
| `federated/primary/` | Local k/silhouette/stability/fingerprints, mapping coverage, separate death/hospice/composite models, actual released model denominators/formulas, competing events, risks, transported prediction |
| `federated/scores/` | Modified-SOFA nested comparisons, paired increments, calibration/Brier, decision curves and model metadata |
| `federated/review/` | ARF types/timing/contribution, code transitions/outcomes, cancer proxy components, COVID/full-code/first-hospitalization subgroups and source feasibility |
| `federated/sensitivities/` | All eight independently refederated variants with their own frozen version references |
| `federated/absorbing/` | Descriptive terminal-state fingerprints and early-event counts; no outcome regression |

The package also captures PH diagnostics and waterfall conflict checks. Clinical
characteristics from baseline local clusters are explicitly separate from shared
phenotype tables. Recovered SOFA diagnostic scores do not replace the original
scores used by these model exports. A manifest records that distinction.
Do not return the surrounding `local/` folder or an entire analysis run.

## Coordinator-only collection

```sh
Rscript diagnostics/site_return/collect.R NEW_OUTPUT SITE_A_SHARE SITE_B_SHARE SITE_C_SHARE
```

Collection reads aggregate return folders only, checks distinct sites and a common
frozen primary consensus, and retains site rows and missing/suppressed statuses.
It does not infer zeros, pool suppressed counts or average site medians. The
existing meta-analysis and frozen-prediction scripts remain responsible for
inferential synthesis; this collector assembles supplement diagnostics.

Review `site_return_requirements.csv` against actual manifests. A capture-complete
package can still have non-estimable models, unavailable oncology/APACHE data or
pending author confirmations. These must be reported, not treated as successful
comparisons or proof that treatment limitations were absent.

The diagnostic return includes global and term-level Schoenfeld tests. For older
saved Cox fits without a design matrix, a diagnostic-only local refit must
reproduce the original coefficients before tests are reported; future site fits
retain the matrix locally. Existing source fits are never overwritten.
Because NU diagnostics suggested nonproportional phenotype associations, every
site also returns exploratory posthoc0–7 and7–28-day association fits and actual
risk-set metadata. These are time-window sensitivity estimates, not replacement
primary HRs or proof that proportionality holds within each interval.

For a versioned diagnostic extension of an already frozen return, use
`model_diagnostics.R` and `add_model_diagnostics.R` in `diagnostics/site_return/`.
The new packet retains the parent manifest hash. The current verified packets
are `output/run_site_returns_20261026/site_packets/SITE/share/`; their combined
supplement is `output/run_site_returns_20261026/reports/site_supplement/`.

The adapted APACHE extension is now automatic in new reviewer runs. Component
capture is also automatic in the standard supplement phase. The complete score
comparison is returned under `federated/apache_adapted/`, with separate complete,
all-case missing-indicator and same-patient SOFA cohorts. The collector rejects
mixed extension phases or incompatible score specifications. This is not standard
APACHE II/IV. Current packets including this extension are
`output/run_APACHE_adapted_20261030_v2/site_packets/SITE/share/`; the combined
supplement is `reports/site_supplement/apache_adapted/` under that run.

## Publication Table 1 and Table 2

New core derivations capture `share/descriptives/` with overall eligible-cohort
and landmark characteristics plus explicitly provisional local-cluster tables.
The complete frozen federation return additionally captures
`federated/descriptives/`: Table 1 for eligible and landmark cohorts separately,
and Table 2 for shared P1–P4, insufficient-information I, and Overall. Never use
provisional local clusters as shared phenotypes. `federated/cohort_audit/` retains
screening, annual denominators, patient counts and outcome/source audits.

Descriptors cover demographics, coded noncancer comorbidities, cancer proxies,
severity, support and ARF evidence/unknown hours, code status and outcomes.
Whole-stay outcomes and hospital LOS are descriptive; post-landmark28-day events
have separate labels. Optional stage, histology, ECOG and treatment fields have
explicit availability flags; timing must be confirmed before covariate use.

Continuous exports provide observed/missing N, sum, sum of squares, mean/SD and
site median/IQR. Exact pooled means/SD use sufficient statistics. Site medians
are never averaged. Released fixed-bin counts support exact integer quantiles
or continuous quantile bounds; otherwise retain site median/IQR. Categorical
counts have explicit total/observed/missing denominators. Linked small-cell
suppression can prevent pooling; missing or suppressed values are never zeros.
Site-distinct patient counts cannot be summed into a unique cross-site count.

The main collector automatically writes `publication_tables/` for consistent
new returns. Current verified three-site packets are
`output/run_descriptives_20261030_v2/site_packets/SITE/share/`; publication tables
are in `output/run_descriptives_20261030_v2/reports/`. Only share folders are
returned, after institutional disclosure review.

Daily competing-event curves are captured automatically from the frozen full
assignments into `federated/primary/outcomes/full/cumulative_incidence_daily.csv`.
These are Aalen–Johansen probabilities at days0–28 after landmark, with clustered
SE and released risk-set bands; linked small event/risk-set partitions remain
withheld. The coordinator retains site-specific curves rather than averaging
incomplete curves into a pooled estimate. Collection rejects mixed curve phases.
Current packets including these curves are under
`output/run_manuscript_figures_20261030_v2/site_packets/SITE/share/`.

The survival/time-dependent discrimination extension is a separate local pass:
`site_survival_tdauc.R SITE PRIMARY_ROOT NEW_SITE_OUTPUT` in
`diagnostics/manuscript_figures/`. It requires the frozen primary, reviewer-score
and adapted-APACHE evaluation caches for that site (`CLIF_REVISION_OUTPUT` and
`CLIF_APACHE_OUTPUT` can select their roots). Add the resulting aggregate share
folder with `add_survival_to_return.R SITE SURVIVAL_SHARE OLD_SHARE NEW_SHARE`.
Return only the new share folder. Collection rejects mixed survival/tdAUC phases.
Current packets are under `output/run_survival_tdauc_20261030_v2/site_packets/`.

Whole-cohort cause-specific KM curves start at ICU admission; shared-phenotype
curves start at72h. Exits are censored in KM; no post-discharge all-cause survival
is inferred. tdAUC exports contain cumulative/dynamic AUC(t), pointwise patient-
bootstrap intervals and fixed-grid ROC at7/14/28d under both competing-control
definitions. They evaluate frozen28-day risk scores over time; the earlier
single-horizon AUROCs remain separately labeled. See the figure README for
missing-outcome assumptions and the distinction from competing-event risks.

The current regenerated exact-count packets are
`output/run_exact_aggregates_20261030/site_packets/SITE/share/`.
They supersede the historical count/descriptor/curve tables and contain no old
count intervals. Counts are regenerated from local caches, never inferred from
old interval endpoints. They include Table1/Table2, daily KM/AJ/tdAUC summaries,
all primary/sensitivity cohort counts, saved model denominators, severity
completeness, hourly state counts and saved exact supplement diagnostics.
Clinical model fits and the frozen phenotype consensus remain unchanged.

`diagnostics/exact_aggregates/README.txt` documents capture and aggregate-only
collection. These count packets accompany the preserved frozen inferential
results; they do not claim previously unavailable estimates are estimable.

Pooled survival adds a two-pass aggregate return because existing daily curve
estimates alone cannot reconstruct the exact combined Kaplan–Meier curve.
See `diagnostics/pooled_survival/README.txt`: sites return an aggregate relative-
event-duration registry, receive the common grid, then return exact risk/event
counts and200 patient-cluster bootstrap count replicates. The coordinator pools
counts before computing survival. No identifiers or calendar timestamps move.
Current pooled figures are in `output/run_pooled_survival_20261030_v2/reports/`.
These are unadjusted combined-cohort curves; retain site-specific displays and
adjusted effect meta-analysis. They do not replace competing-event incidence.


Pooled survival and discharge: required two-pass site returns

Each institution must run `diagnostics/pooled_survival/site_times.R` first.
Return `event_duration_counts.csv`, `cohort_counts.csv`, and `manifest.json`.
After receiving the coordinator's shared grid, run `site_risksets.R` and return
`risksets.csv`, `bootstrap_risksets.csv`, and `manifest.json`. Both passes include
the admission-time whole cohort and landmark phenotype cohorts. Endpoints must
include death, death_or_hospice, live_discharge, and any_exit. The last two
support Aalen-Johansen discharge incidence with death/hospice competing.
The 200 patient-cluster bootstrap replicates use identical within-site patient
weights across groups/endpoints. Returns contain exact counts and relative
durations, with no patient identifiers, individual risk scores or calendar dates.

Bundle both passes into the final main packet using:

```sh
Rscript diagnostics/pooled_survival/add_to_return.R SITE TIME_SHARE RISK_SHARE OLD_RETURN NEW_RETURN
```

The final packet includes `federated/pooled_survival/times/` and
`federated/pooled_survival/risksets/`; its manifest must set
`pooled_survival_discharge_complete=true`. The collector validates local AJ/KM
checks, all four endpoints, file hashes, consensus, common grid and bootstrap
count. This extension requires a coordinator round trip; the initial local run
alone cannot produce the final risk-set return. Do not pool finished site curves.

For final publication collection, set `CLIF_REQUIRE_POOLED_SURVIVAL=true` when
running `diagnostics/site_return/collect.R`; a missing site extension then fails
collection. Earlier-phase packets remain collectable without this final-phase
requirement.


Pooled four-model ROC curves require the additional aggregate return generated
by `diagnostics/pooled_roc/site_capture.R SITE PRIMARY_ROOT APACHE_ROOT NEW_SHARE`.
Return score_histograms.csv, cohort_counts.csv, local_validation.csv and
manifest.json. These include fixed common score-bin counts for clinical
covariates alone and plus modified SOFA, adapted APACHE or phenotype, using
identical matched complete-case hospitalizations at days 7/14/28 after72h.
Both hospital-death and death/hospice endpoints and 200 paired within-site
patient-cluster bootstrap count replicates are included. The coordinator uses
`diagnostics/pooled_roc/central_plot.py` to pool counts before constructing ROC.
These are additional required returns for this figure: the earlier site ROC
curves at a common false-positive-rate grid cannot reconstruct pooled ROC.
See diagnostics/pooled_roc/README.txt for quantization and estimand details.
