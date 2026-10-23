# Reviewer concern implementation map

This tracks implementation, not new scientific results. All analyses require fresh source-table execution and clinical review. “Implemented” does not mean the revised manuscript has been updated or the corresponding result is established.

| Concern | New implementation or required follow-up |
|---|---|
| Early death/discharge in phenotype inputs | ICU landmark survivor derivation; explicit early-event counts; separate absorbing-state sensitivity; terminal sensitivity excluded from outcome models |
| Hospital-versus-ICU landmark | Primary origin is first ICU arrival; follow-up begins at ICU +72 hours; earlier landmarks are sensitivities |
| ARF-onset anchoring | First-detected positive physiology anchors a sensitivity with its own 72-hour landmark; missing/nonqualifying evidence audited; detection is not assumed biological onset |
| Event timestamp errors | Validated in-hospital death time or documented discharge proxy; hospice uses discharge time; no last-vital fallback |
| Missing physiology interpreted as no ARF | Separate positive, nonqualifying, and unknown evidence for hypoxemia and hypercapnia |
| HFNC/CPAP/NIPPV pooled | Distinct support states and cumulative hours |
| ARF detection depends on FiO₂ | Missingness audits, support-only sensitivity, low-flow FiO₂/SF sensitivity with stated assumptions |
| Cohort called ARF despite no ARF requirement | Broad lung cancer ICU cohort explicitly defined; ARF/high-support subgroup analysis |
| Outcome-informed consensus/reference | Hourly respiratory profiles only; deterministic respiratory-input reference; no mortality/SOFA/vasopressor/LOS matching |
| Temporal ordering discarded centrally | Hour-specific joint distributions, rather than whole-window proportions |
| Variable local k and matching bug | Local silhouette selection; rectangular assignment considers every prototype; unmatched profiles retained explicitly |
| Reproducibility | Local silhouette tables, patient bootstrap stability, mapping distances, consensus LOSO |
| Unclassified patients | Explicit landmark, support-missingness, and mapping audits; every eligible sequence, including entirely unknown support, has a local assignment |
| Prognostic superiority | Training-fold phenotype learning and nested A/B/C comparisons; AUC, calibration, Brier, decision curves and paired increments; remove superiority claims if unsupported |
| Baseline severity absent | Modified measured SOFA proxy in models; actual missingness and model exclusions; standard SOFA/APACHE claims require more data |
| Code status and treatment limitation | Admission and window-end status, new DNI/comfort status, full-code sensitivity; unknown and DNR-only remain separate |
| Cancer information | Explicit POA proxies, noncancer Charlson, optional descriptive oncology table; site availability and provenance still required |
| Composite described as mortality | Separate death, hospice, and composite outcomes and labels throughout new exports; manuscript wording remains to revise |
| Competing live discharge | Cause-specific estimand stated; Aalen–Johansen curves for mutually exclusive terminal hospital events |
| ICU LOS overlaps exposure and is truncated by death | Total LOS models removed; full ICU-free-day outcome requires additional validated follow-up |
| 24/48-hour utility | Separate derivation and model passes for each landmark; avoid interpreting different eligible cohorts as paired comparisons |
| Repeated patients | Hospitalization unit, patient bootstrap, patient-grouped validation, patient-clustered Cox and Aalen–Johansen uncertainty with separate hospitalization episodes, first-hospitalization sensitivity |
| Pollution covariates and double-counted cancer | No pollution/geography dependency; noncancer Charlson plus explicit cancer proxy only in retrospective associations |
| Model denominators differ | Model-specific actual N/events/patient count/formulas; missingness and exclusion audits; no cohort-N substitution |
| Table percentages and IMV/vasopressor labels | Distinct variables generated from source records; numeric n/denominator fields; no manual count repairs |
| Site count rescaling and implausible trends | Historical returns never loaded; no central rescaling, missing-year interpolation, or conversion of suppressed counts to zero |
| Meta-analysis uncertainty | REML with modified Hartung–Knapp, τ², I², prediction intervals, site-level effects |
| Missing phenotype main-effect pooling | Explicit phenotype-vs-shared-reference effects and separate pooling script |
| COVID inconsistency | Defined POA U07.1 exclusion with separate results; manuscript must use that exact definition |
| Ethics/reporting/literature | RECORD/STROBE, board/approval/waiver details, disclosures, data availability, and revised literature discussion remain manuscript work |

## Before site distribution

Pilot with real CLIF data; inspect device and disposition mappings, event proxy use, code-status coverage, six-hour carry, support missingness, SOFA proxy completeness, runtime, suppressed-profile overlap, and cluster interpretability. Missing support does not exclude anyone. Freeze the source commit and parameters across sites after the pilot. Withdrawing a site or changing a mapping requires a new consensus version and a new model pass.

The existing public/submitted code has been preserved in the archive. The new pipeline is a new analysis and may produce different cohort sizes, cluster counts, labels, and effect estimates. Historical phenotype names should not be imposed on its results.

## UCMC pilot: 2026-10-21

The adopted primary is the joint support–ARF ICU +72-hour landmark analysis. The fresh pilot contains 849 cohort hospitalizations: 650 enter the primary derivation (six local clusters), and all 849 enter the absorbing sensitivity. Post-primary-landmark outcomes are 152 hospital deaths, 72 hospice discharges, and 426 other live discharges. No support-coverage exclusion is applied. The first-detected-ARF sensitivity has 276 eligible hospitalizations and five local clusters. All nine derivations complete.

Only 99 primary hospitalizations have a complete measured-component modified SOFA proxy; the complete-case adjusted pilot models produce estimation warnings. Some local clusters also show limited bootstrap stability. Both measured-component and missing-domain-zero predictive validation are incomplete; no performance claim is exported. These are unresolved interpretation issues, not evidence of prognostic improvement.

The conservative small-cell suppression of hourly joint-state distributions leaves insufficient profile overlap for shared consensus matching in this pilot. Do not lower disclosure protection or interpret suppressed blanks as zero. Complete local profiles may support provisional within-site models, explicitly labeled as such; they do not establish multicenter consensus. A permissible aggregate profile representation remains necessary before site-wide harmonization.

## Unknown-support investigation: 2026-10-22

A source-level audit of the unchanged UCMC primary derivation found 14,525 unknown-support hours among 46,800 hours (31.0%). Every populated device category maps successfully. Of unknown hours, 11,272 (77.6%) have no respiratory-support row in that hour and 3,253 have source rows without an identifiable device. The missing-device rows have no recoverable device name or mode category. An older local respiratory extract supplies no additional populated device labels at the missing timestamps.

The problem occurs in both ICU and ward observations and increases markedly after 2019. ICU unknown hours rise from 2.8% in 2018 to 49.2% in 2024; ward unknown hours rise from 6.2% to 77.6%. This is consistent with a documentation or extraction coverage change, not proof of a specific ETL defect. The initial ICU facility is the same throughout; the pattern also occurs across ICU types. Raw EHR flowsheets and the producing ETL are unavailable, preventing source verification.

An exploratory, outcome-free repair distinguishes documented supplemental oxygen with unspecified device from completely unknown support. Valid positive oxygen flow (0 < lpm_set <= 60) appears directly in 3,001 unknown hours. A further 4,717 hours have eligible six-hour forward carry of positive flow, bounded by admission and the landmark, stopped by a newer valid zero-flow reading or conflicting named-device record. The current-hour rule represents any positive flow observed within that hour; carrying is used only when that hour has no valid flow reading. No backfill, low-flow/HFNC/NIPPV/IMV guess, estimated FiO2, physiology change, or patient exclusion is applied. Remaining fully unknown hours are 6,807 (14.5%); exact device information is still missing in the oxygen-unspecified hours.

This exploratory alphabet selects seven clusters (silhouette 0.314; ARI 0.692 versus primary). Of the original 212 hospitalizations in the unknown-dominant cluster, 136 move to an oxygen-unspecified-dominant group and 62 move to a residual unknown-dominant group; 14 move to other groups. These are descriptive diagnostic findings, not validated replacement phenotypes or evidence of better prognostication. Primary analysis definitions, assignments, and outcomes remain unchanged.

Reproduce with `Rscript diagnostics/01_audit_support.R` against the configured matching completed run. Source data, patient-level diagnostics, and exploratory assignments stay in the ignored output directory. The next available-data option is an explicit oxygen-unspecified sensitivity; recovering exact modality requires improved source device documentation/extraction.

## CLIF waterfall comparison 2026-10-22

The earlier missingness diagnostics used recorded devices plus bounded six-hour carry, without applying CLIF waterfall in this analysis. The subsequent comparison uses the official CLIFpy 2.1.0 waterfall with numeric bfill disabled. Input and output are separately bounded at each landmark or terminal event, with the same six-hour pre-origin lookback. Standard categorical heuristics can use within-window lookahead; device forward fill has no six-hour expiry. Empty boundary markers extend the hourly scaffold through the requested window. They add no observed device values.

All nine prior derivations complete with unchanged cohort membership. An additional 72-hour comparison holds original ARF states fixed and changes only respiratory support. Among the same 650 primary hospitalizations, UNKNOWN support falls from 14,525/46,800 hours (31.0%) to 6,166 (13.2%); the same 64 hospitalizations remain entirely unknown. Full device-and-numeric waterfall selects seven clusters (silhouette 0.391; ARI 0.385 versus baseline); waterfall support with original ARF states selects six (silhouette 0.403; ARI 0.405). This is a paired preprocessing comparison, not a new validated primary result.

Filled FiO2 changes 4,758 hourly physiology states and increases hospitalizations with positive ARF evidence from 215 to 255. Baseline ICU severity and code-status covariates remain frozen for comparison. All 152 post-landmark hospital deaths and 72 hospice discharges remain unchanged. Full-waterfall adjusted models continue to have estimation warnings with the 99 measured-component complete cases. Predictive validation has not been rerun with waterfall inside the training folds.

The full-waterfall unknown-dominant group contains 81 hospitalizations; the corresponding group in the support-only-change comparison contains 86. Some full-waterfall groups remain unstable, including a mean bootstrap Jaccard of 0.364 for cluster 7. Improved coverage includes inferred/filled exposure, not additional directly observed device documentation: 29,039 primary hours are filled or inferred support states.

Source-context QA flags 117 records across 107 hours in 23 hospitalizations with positive raw oxygen flow but waterfall-assigned room air. It also identifies 1,870 hours containing a filled FiO2 of 0.21 on nasal-cannula records. These are review flags, not proof of incorrect source documentation, and warrant review before treating filled FiO2 as observed physiology. The pre-waterfall audit alone should not be interpreted as evidence of an ETL defect.

Reproduction: `bash diagnostics/02_compare_waterfall.sh`, with a fresh output directory and an isolated Python environment containing CLIFpy 2.1.0. Full methods, boundaries, versions, and interpretation are described in `diagnostics/waterfall/README.txt`. Primary baseline files and config remain unchanged.

## NU and Rush comparison 2026-10-22

The same nine baseline derivations and ten waterfall comparisons now complete at NU and Rush. The 2018–2024 adult POA lung-cancer ICU cohorts contain 2,325 NU and 1,359 Rush hospitalizations; primary +72-hour landmark cohorts contain 1,691 and 969, respectively, compared with 650 at UCMC. No support-coverage criterion is applied. NU has three unclassifiable outcomes; Rush has 60 (59 source dispositions marked Needs Coding by HIM and one Other transfer). These are distinct from unknown support.

Unknown primary hours before/after waterfall are UCMC 31.04%/13.18%, NU 3.16%/0.33%, and Rush 17.69%/13.15%. Wholly unknown stays after processing number 64, 3, and 126. Among those stays, no respiratory-table rows occur in the window for 17 UCMC, all three NU, and 125 Rush stays. Positive source oxygen flow is present in 43 of the 64 UCMC wholly unknown stays. Consequently the residual source-information gaps differ substantially across sites; this comparison does not prove any particular ETL defect.

The conservative waterfall-support/original-ARF arm selects six UCMC, five NU and five Rush local clusters. Four broad support patterns recur (room-air recovery, low-flow oxygen, HFNC/noninvasive support and persistent IMV). Additional clusters include an unknown-dominant UCMC group (86), an entirely unknown Rush group (126), a small NIPPV-dominant UCMC group (25; mean Jaccard 0.403), and a small NU room-air group (32; mean Jaccard 0.181) with predominantly unknown hypoxemia states. Cluster numbering is local and unmatched. These findings require outcome-free phenotype harmonization and missingness sensitivity before pooling associations; bootstrap stability of an unknown group is not clinical validation.

Raw ARF detections number 215/650, 639/1,691 and 337/969; full numeric waterfall increases these to 255, 751 and 392. Post-landmark hospital death/hospice percentages are UCMC 23.4%/11.1%, NU 17.6%/13.3%, and Rush 9.4%/17.2%. These are descriptive risks over the remainder of hospitalization. Complete measured-component SOFA is available for only 99, 219 and 62 primary stays, respectively. No new pooled adjusted effect or external predictive validation is asserted.

Input auditing uncovered two reader bugs: date-only study boundaries became NA, and a projected-out optional hospitalization ID in Rush code_status was used for Arrow filtering. Rush was restarted with explicit timestamps; its aborted v3b run is excluded. UCMC and NU complete included cohorts were verified already within the intended dates. A local Rush code-status projection preserved all retained rows and values. After completing the calculations, core utc() and the selected-ID filter were fixed; invalid/reversed boundaries now stop configuration loading. New regression tests and the synthetic analysis suite pass. An additional audit verifies that fixed date-only parsing reproduces all three completed eligibility/landmark sets, including removal of the extra 2017/2025 stays in the aborted Rush preparation.

The exact analytical source before these reader fixes is preserved in output/run_multisite_comparison_20261022/local/analytical_source_before_input_fixes.tar.gz. Completed results retain their original source fingerprints; reproducing cached runs requires this snapshot, while new core runs require fresh IDs. Reproduction and schema-adapter details are in diagnostics/multisite/README.txt. Summary, all sensitivities, local profiles, source audits and verified figures are in output/run_multisite_comparison_20261022/. Those aggregates are exploratory local results, not disclosure-reviewed share packages.

## Federated primary-analysis rehearsal 2026-10-22

A separate rehearsal now completes the local-discovery → aggregate-fingerprint → central-consensus → frozen local-outcome → association-pooling → held-out-site workflow. Site workers reuse only their own validated waterfall-support/original-ARF hourly panels. All 3,310 primary +72-hour landmark hospitalizations remain in the main analysis (650 UCMC, 1,691 NU, 969 Rush), including an explicit insufficient/ambiguous-information assignment. Patient identifiers, hourly sequences, clinical rows and individual predictions stay in site-local directories. Coordinator scripts consume allowlisted aggregate exports and derived coordinator artifacts. This is logical process separation on one machine, with disclosure review still required before institutional exchange.

For this new rehearsal local OM/Ward clustering uses a common constant substitution cost of 2 and indel of 1, rather than site-estimated TRATE matrices. Silhouette selects local k=8, 4 and 5 at UCMC, NU and Rush. One hundred patient bootstraps characterize stability. Six-hour modal joint/channel distributions and support-change fingerprints are interval-valued; size and stability are banded. The coordinator reconstructs feasible probability templates and uses equal site weighting. Candidate global k=2..6 requires each group represented in at least two training sites; silhouette selects among supported candidates. Four shared phenotypes are selected in the three-site fit and independently in every two-site discovery: recovery toward room air, low-flow oxygen, HFNC/noninvasive support and invasive ventilation.

Frozen per-hospitalization template classification allows many local groups to converge on the same clinical phenotype, avoiding the prior one-to-one constraint. Fewer than twelve known-support hours, distance >0.75 or margin <0.01 yields I; these are assignment statuses rather than exclusions. Reported full-versus-held-out assignment agreement falls in the exported 95–100% interval at each site. The method is federated consensus/prototype classification, not an exact pooled OM/Ward solution. Six-hour aggregation and interval reconstruction are explicit information-loss/privacy tradeoffs.

Core Cox models adjust age, sex, race grouped White/Black/other-or-unknown, admission year, noncancer Charlson, a POA advanced-cancer proxy and admission code status. They use patient-clustered uncertainty, follow-up capped at 28 days after landmark, and separate baseline hazards for classified versus I cases. REML random-effects meta-analysis with modified Hartung–Knapp intervals yields composite death/hospice HRs versus room-air recovery of 1.76 (95% CI 1.02–3.03) for low-flow oxygen, 2.60 (1.54–4.40) for HFNC/noninvasive support and 2.61 (1.11–6.15) for invasive ventilation. The latter is heterogeneous (I² about 62%; prediction interval about 0.62–11.03). Small-event release rules leave only one exportable site for separate death contrasts, preventing meta-analysis; hospice uses two sites and has broad intervals. Sparse complete-case SOFA models do not supply enough usable site contrasts for severity-adjusted pooling.

Each held-out site receives profiles and ridge-logistic predictive coefficients learned exclusively at the other two sites, with no held-out refitting or recalibration. Baseline→baseline-plus-phenotype AUCs are UCMC 0.669→0.736, NU 0.613→0.755 and Rush 0.614→0.724. Patient-bootstrap intervals (200 replicates, conditional on frozen training) support improvement over baseline covariates. A harder exploratory comparator added during development uses support duration, observed ARF and missingness summaries. Adding phenotype beyond those summaries changes AUC by −0.0059 (95% CI −0.0124 to −0.0006) at UCMC, +0.0068 (0.0026–0.0107) at NU and −0.0078 (−0.0203 to 0.0071) at Rush. The results support reproducible descriptive phenotypes and outcome associations, while not establishing consistent predictive superiority over simple respiratory summaries.

Local cohort/prediction audits, export-header gates, held-out input checks, frozen consensus content hashes, source snapshot hashes and all twelve outcome-version references pass. No individual predictions or rows are in the exchange. The dataset had already been examined; this is an internal-external development rehearsal, not independent prospective validation. No causal effect is claimed. This federated run covers the 72-hour primary workflow; previous 24/48-hour, absorbing, ARF-onset and alternative-physiology pilot sensitivities have not all been refederated under the new prototype classifier.

Implementation and reproducible fresh-output runner: diagnostics/federated/README.txt and diagnostics/federated/run_rehearsal.sh. Aggregate report, six visually checked figures, all site/pooled effects, calibration/Brier estimates, suppressed competing-risk summaries, membership bands and audits are in output/run_federated_rehearsal_20261022/reports/. Completed source snapshots are under that run's provenance/. Existing primary manuscript and pilot outputs are preserved.

## Three-site reviewer-gap extension: 2026-10-22

The extension in `diagnostics/revision20261022/` completes eight independently
refederated inferential sensitivities (24/48-hour landmarks, first-detected ARF,
support-only, raw-source low-flow FiO2, 0/12-hour evidence carry, and no-waterfall
support), with 32 frozen consensus versions and 96 local outcome passes. A
separate absorbing-state federation uses descriptive D groups, with no outcome
regression. Its terminal-state transitions are explicitly unranked. The original
primary phenotype and model artifacts remain unchanged.

New analyses characterize ARF evidence type/timing/late positivity and admission
to landmark code-status documentation. They report separate advanced-cancer
proxy components, outcomes by code status, actual local model denominators and
export-safe continuous-variable observed/missing denominators. Strict full-code,
COVID, first-hospitalization and ARF/high-support sensitivities are run with the
frozen primary classifier. Full-code contrasts remain too sparse to pool.

Direct severity comparison now includes first-24-hour measured-component
**modified SOFA**, score alone, baseline plus score, simple summaries including
hour-72 support, and phenotype additions. Models use the same evaluation cases
within cohort/endpoint, with two-site training and held-out-site evaluation.
AUC/Brier, calibration, paired increments and exploratory decision curves are
provided for death and the composite. Only 380 primary hospitalizations have all
six measured domains. Creatinine-only renal scoring means this is not validated
standard SOFA. An explicitly labeled missing-domain-zero/missing-indicator
sensitivity retains all primary cases. Gains beyond severity and summaries are
small/uncertain at UCMC and NU and negative at Rush in complete cases.

**APACHE II/IV comparison remains unavailable**, rather than completed with a
proxy. No recorded scores exist; required physiology and chronic-health/admission
context are missing/unverified. UCMC has no hematocrit found in audited source
labels/LOINCs; NU has some arterial hematocrit mapped to `other`. Availability
and explicit missing comparison rows are exported. Do not claim superiority over
APACHE or inadequacy of established scores.

Support-only versus joint classification ARIs are approximately .859/.888/.909
at UCMC/NU/Rush, supporting substantial modality dominance. ARF-channel outcome
effects are imprecise. Primary four-group reproducibility does not imply an
invariant solution across all sensitivities: the 24-hour full solution selects
five groups, and some alternative-ARF held-out solutions select three.

The completed report, five figures, tables, reviewer-response/manuscript-text
drafts and reporting crosswalk are under
`output/run_revision_gaps_20261022/reports/`. These are three-site development
results, not a completed ten-site revision. Stage/ECOG/therapy intent, missing
care-preference documentation, unclassifiable source outcomes, scarce separate
death effects, author ethics/disclosure details and execution at remaining sites
remain unresolved. Submitted Word files were not overwritten. No institutional
exchange or disclosure approval is implied.

The final source code-status catalogue audit finds **no comfort-care category at
any site**, no explicit DNI category at NU, and Rush DNR-named orders mapped to
canonical DNR/DNI. These are unverified care-limit semantics, not evidence that
comfort care/DNI preferences were absent. The detailed status analyses require
this qualification. Catalogue exports are included in the new report package;
no source mapping was silently rewritten.
