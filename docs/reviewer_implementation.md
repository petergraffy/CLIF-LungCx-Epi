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

## UCMC source refresh verified: 2026-10-26

A newer respiratory table is present in Box under `03-CLIF-2.1/2.1.0` (October
23 update, about 225 MB). The configured Desktop table is still the prior copy;
its matched source records equal the saved snapshot after sorting. The source
refresh audit keeps the same 650 landmark hospitalizations and other clinical
inputs fixed. Eligible/landmark membership remains 849/650, with no changes in
matched event times/types, ICU origins, demographics or hospital admission/
discharge times.

Raw-source plus six-hour-carry UNKNOWN hours decrease from 14,525/46,800 (31.04%)
to 3,233 (6.91%). All 64 previously entirely unknown stays now have source device
information and none remains entirely unknown. With the same bounded waterfall,
UNKNOWN hours decrease from 6,166 (13.18%) to 30 (0.064%), across 11 stays. In this
paired respiratory-only comparison, raw-evidence ARF detection changes from 215
to 237 hospitalizations. The conservative waterfall comparison preserves updated
raw physiology rather than using filled FiO2.

This establishes major source coverage restoration, not a new pooled analysis.
Producing-ETL inference versus bedside-observed provenance still needs to be
represented accurately. The updated source has unannotated timestamps; overlap
with prior events is highest without a shift. Waterfall conflict flags are
recorded in the audit. Existing UCMC/NU/Rush/federated artifacts are unchanged.
A fresh UCMC run against the Box tables, regenerated fingerprints, versioned
consensus, and new outcome/pooling/transport passes are required before revised
results can be reported. The current private configuration still points to the
old Desktop copy; do not reuse its cached run as refreshed data.

Implementation: `diagnostics/etl_refresh/README.txt`. Paired report:
`output/run_UCMC_etl_refresh_audit_20261026/etl_refresh_report.txt`.

## Completed three-site refresh: 2026-10-26

The fresh UCMC preparation reloads all required updated Box tables with run ID
`ccx_box_20261026_joint`; its primary cohort remains 650 stays, with raw ARF
positive evidence in 237 and complete measured-component modified SOFA in 103.
The default private UCMC configuration now points to this Box source/run ID;
the previous configuration is backed up in the ignored orchestrator directory.
NU/Rush validated source caches are reused, and their panels are independently
verified unchanged. All three sites' local clusters, stability assessments,
shared definitions, assignments, outcome models and sensitivities are refitted.

The primary retains all 3,310 stays. Local k is 4/4/5 at UCMC/NU/Rush; the global
solution selects four groups in the full fit and in every held-out-site fit.
Its consensus ID is `6245d53b0469ab7d1a5b`. Predominant support labels are room
air, low-flow oxygen, HFNC and IMV. Core composite HRs versus room-air support
are 1.82 (95% CI 1.16–2.86), 2.61 (1.56–4.36), and 2.81 (1.41–5.59). The IMV
association remains heterogeneous (I² 46.4%; prediction interval 0.97–8.17).
Separate death effects now have two contributing sites, but very wide intervals.

Held-out baseline-plus-phenotype AUCs are .751/.756/.728. Increment beyond
simple summaries is inconsistent across sites; complete-component modified-SOFA
plus summary comparisons likewise do not establish universal added predictive
value. SOFA remains modified/incomplete; APACHE and oncology/care-preference
limitations remain. One small UCMC local group has a released mean-Jaccard band
of .5–.6, so reproducibility of global k does not imply every local group is
stable. Physiological unknown states remain explicitly represented.

All eight inferential sensitivities were refederated through full and three
held-out consensuses, local outcomes, pooling and prediction (32 sensitivity
consensuses; 96 sensitivity outcome passes). The separate absorbing analysis
remains descriptive. Primary plus sensitivity outcome passes total 108. The
24-hour full solution selects five groups; some physiology variants also differ
in held-out k. No historical k or clinical labels are imposed on separate fits.
All member-retention, training-split, frozen-prediction/hash, terminal-boundary,
original-artifact and export audits pass. Original results are preserved; MIMIC
remains excluded. Historical numerical manuscript drafts were not regenerated.

Current report:
`output/run_three_site_refresh_20261026/reports/three_site_refresh_report.txt`.
Primary results/figures: `output/run_federated_rehearsal_20261026/reports/`.
Sensitivity/score results: `output/run_revision_gaps_20261026/reports/`.
Reproduction/parameter details: `diagnostics/refresh20261026/README.txt`.

## SOFA component audit: 2026-10-26

First-24h source recovery increases complete measured-component modified SOFA
from 384 to 416 of 3,310 stays. UCMC calculated MAP from paired SBP/DBP raises
CV coverage 306->650 and complete totals 103->116. NU source arterial PO2
mappings and panel FiO2 raise respiratory coverage 444->489 and totals219->238;
additional evidence changes30 existing respiratory scores and18 complete totals.
Rush totals remain62; simultaneous GCS components and numeric strings provide
no missing CNS recovery. A separate recent measured pre-ICU weight sensitivity
raises complete totals to441. A further explicitly estimated-PaO2-from-SpO2
sensitivity reaches1066, not measured or validated standard SOFA. Bilirubin
and Rush GCS remain major gaps. Source drug-category flags are documented.
These are candidate audit scores; original primary outcomes, comparisons and
phenotype artifacts are preserved and require fresh model passes before updated
performance or severity-adjusted associations can be reported. See
`diagnostics/sofa_completeness/README.txt` and
`output/run_sofa_completeness_20261026/reports/sofa_completeness_report.txt`.

## Site supplement/return audit: 2026-10-26

The editor/reviewer letter was reread against the actual outputs. A distribution
gap was found: several diagnostic/reviewer tables were rehearsal-only and were
not returned by the default site command. The source/baseline `derive` entry point
now captures a portable supplement automatically; `supplement` captures one for
an existing derivation. Arbitrary configured site names and CSV/parquet/FST
sources are supported. Source recovery diagnostics are not silently substituted
into frozen outcome models. Exact counts/candidate scores remain local.

The reviewer runner now assembles a per-site return containing allowlisted
primary discovery/models, transported/nested score performance, ARF/code/cancer
characterization, all eight sensitivity versions, absorbing-event descriptions,
waterfall conflict checks and PH/denominator diagnostics. The coordinator-only
collector validates file hashes, site identity, common consensus and required
families; it reads aggregate return folders only. Suppression statuses remain
explicit and counts are not centrally repaired. Three current sites were captured
into `output/run_site_returns_20261026/`; candidate SOFA scores exactly reproduce
the previous audit. The verified collection is in that run's `verified_collection/`.

A concern-by-concern checklist is `docs/site_return_requirements.csv`. Administrative
and provenance fields are requested in `config/site_metadata_template.csv`;
blank IRB/waiver, order semantics or ETL fields remain pending site confirmation.
Molecular subtype and treatment type are now accepted in optional oncology
inputs/availability summaries. Model passes retain exact denominator metadata
locally and explicitly report reference-absent subsets.

The active protocol and site workflow now distinguish the adopted bounded-
waterfall/common-cost/template procedure from the earlier raw/TRATE/one-to-one
baseline procedure. The older protocol is archived. This does not turn the
three-site central orchestration into an already executed ten-site federation.
Remaining scientific/reporting limitations are unchanged: incomplete modified
SOFA, no validated APACHE, incomplete cancer/preferences, sparse separate death
and full-code effects, no validated ICU-free-day endpoint, source provenance
confirmation, and manuscript/RECORD/STROBE/literature work. No prognostic-
superiority, biological-endotype or causal intervention claim is supported.

The final Cox diagnostic review found that old federated objects lacked local
model matrices. Future runs retain them. Diagnostic refits on all three sites
reproduce original coefficients exactly and return global/term PH tests. NU
core composite PH global p=0.023, with phenotype-term flags; constant-HR
interpretation needs assessment. Every site now also returns exploratory
posthoc0–7/7–28-day effects and actual risk-set metadata. Original models remain
unchanged. Current complete packets are in
`output/run_site_returns_20261026/site_packets/SITE/share/`, with the verified
combined supplement at `reports/site_supplement/` under the same output root.

## External CLIF APACHE implementation reviewed: 2026-10-30

`parkegan/clif_hospice_trend` includes an APACHE II-named helper. At reviewed commit
`e71e6a36130e4e011c652fd93ce9433d22a86cb6`, its actual caller uses24h after
hospital admission, despite a48h helper default. It substitutes hemoglobin,
uses thresholds incompatible with our audited g/dL units without conversion,
defaults missing GCS/other components to0 points, fixes chronic points at0,
omits the renal doubling modifier, and labels>=10 point columns complete.
Synthetic probes reproduce decimal-bin gaps and forward-FiO2 matching under
its rolling-join end settings. Do not import it unchanged or call it APACHE IV.
A corrected explicitly modified APACHE-II-like sensitivity is feasible, while
full validated APACHE availability remains unresolved. No frozen analysis or
comparison was altered. Full review: `docs/apache_implementation_review.md`.

## Adapted APACHE sensitivity completed: 2026-10-30

The accepted recommendation is implemented independently in
`code/lib/apache_adapted.R` and `diagnostics/apache_adapted/`. New comparison
artifacts are under `output/run_APACHE_adapted_20261030_v2/`. Complete adapted
physiology is available for636 stays (171/391/74); the explicit partial-sum plus
missing-indicator sensitivity retains3,310. Matched APACHE-proxy/SOFA evaluation
uses341 stays (96/198/47), with identical patients across models. Original
phenotypes, outcomes, and SOFA model artifacts remain unchanged.

APACHE-proxy and SOFA score-only AUC differences are uncertain at every site.
In the complete adapted-component cohort, phenotype additions beyond the proxy
and respiratory summaries do not improve AUC at any site. All-case results show
small gains at NU/Rush, not UCMC. Modified-score availability does not establish
a validated APACHE II/IV comparison or justify superiority claims.

NU's original source folder moved. The new2.1.1 source supplies metadata only on
exact frozen numeric keys; newer measurements are ignored. Units absent from
matched metadata retain explicit CLIF-canonical assumptions. Source correspondence
counts are returned. Site supplement and reviewer runners now capture the new
component/model/validation outputs; collector checks score specification/consensus
compatibility and refuses mixed score-extension phases. See
`docs/apache_adapted_methods.md` and `diagnostics/apache_adapted/README.txt`.
