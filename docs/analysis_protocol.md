# Active three-site analysis and site-return protocol

This describes the October 26 three-site development workflow. The earlier
baseline protocol is preserved in `archive/2026-10-26/analysis_protocol_baseline.md`.
It used different support preprocessing, matching, and validation methods. Do not
combine its raw-profile consensus returns with the current interval-fingerprint
consensus. The scripts under `code/` prepare source data and baseline comparisons;
the adopted waterfall/template procedure is under `diagnostics/federated/` and
`diagnostics/revision20261022/`. `docs/site_return_workflow.md` specifies returns.

## Cohort and landmark

Include adult hospitalizations beginning in the configured study period with
primary lung cancer (ICD-10 C34 or ICD-9 162) present on admission and an ICU ADT
segment. History-only and secondary lung metastasis alone do not qualify. This
is a lung-cancer ICU cohort, not a cohort requiring measured ARF. The current
study interval is 2018–2024. No geography, pollution, complete demographics, or
minimum support coverage is required. All landmark-eligible stays remain,
including entirely unknown support and ambiguous shared assignments.

The unit is hospitalization. Retain repeat admissions; resample patients with
all their hospitalizations, cluster model uncertainty by patient, and perform a
first-eligible-hospitalization sensitivity. First ICU arrival is t0. Primary
sequences cover [t0,t0+72h); primary members are alive and in hospital through
that landmark. ICU-to-ward transfers remain included using hospital-wide
observations. Do not interpret unknown support as proof the patient is in ICU.

Final outcome and time must be classifiable. Unknown source dispositions and
invalid times are audited separately; no follow-up is fabricated. In-hospital
patient death time is used if within ICU arrival and hospital discharge;
otherwise discharge time is a documented death-time proxy. Hospice and other
live discharge use discharge time. No last-vital fallback or subsequent
outpatient death is used. Source event times, mappings, and timezone semantics
require site inspection.

## Support and ARF inputs

Support states are room air, low-flow oxygen, HFNC, CPAP, NIPPV, IMV, other,
and unknown. HFNC/CPAP/NIPPV are distinct. Original hourly source panels use
highest recorded support, with latest-record tie breaking and bounded six-hour
prior carry. A same-admission six-hour pre-origin lookback is allowed.

The adopted primary uses bounded CLIFpy2.1.0 waterfall support and preserves
ARF evidence from the CLIF source before the analysis-side waterfall. Numeric
bfill is disabled. Standard categorical waterfall rules can use within-window
lookahead and device forward fill has no six-hour expiry. Input/output are
bounded at each patient's landmark/terminal event. Empty scaffold markers add
no observed device values. These are explicit preprocessing assumptions, not
new bedside observations. Raw-support and other carry variants are sensitivities.
"Raw evidence" means source CLIF before our waterfall, not independently
verified raw EHR data or absence of source-ETL inference.

The physiology channel separately retains hypoxemia positive H1, observed
nonqualifying H0, unknown HU, and hypercapnia/acidemia positive C1,
nonqualifying C0, unknown CU. Missing is not negative. Primary positive hypoxemia
is paired arterial PaO2/FiO2<=300, or SpO2<90 on documented room air. FiO2 is
normalized from fraction/percent; only 0.21–1 is accepted. Explicit room air
permits0.21; other devices require source FiO2. Nearest pairing is within1h.
Hypercapnia/acidemia uses arterial PaCO2>=45 and pH<7.35 within2h. Recent
positive/nonqualifying evidence persists up to6h, with positive precedence;
expired evidence becomes unknown. Low-flow FiO2/SF estimation and0/12h carry
are separate sensitivities. Pulse oximetry bias remains a limitation.

## Discovery, fingerprints, and frozen assignment

Each site fits optimal-matching state sequences with common substitution cost2,
indel1, Ward.D2, and silhouette-selected k=2..8 (limited by sample size). This
replaces the earlier baseline TRATE-cost procedure. One hundred patient
bootstraps characterize local stability conditional on the chosen procedure.
No historical number or clinical name is imposed.

Sites return interval-valued six-hour modal joint/channel distributions and
support-change fingerprints, with banded cluster size/stability. The coordinator
reconstructs feasible templates with equal site weighting and size weighting
within site. Candidate global k=2..6 requires each shared group represented in
at least two training sites, then uses silhouette selection. Outcomes, SOFA,
vasopressors and LOS do not enter primary matching, reference selection or
ordering. Local-to-consensus mapping distances remain diagnostic outputs.

Frozen per-hospitalization classification allows many local clusters to share a
phenotype. Every hospitalization with sufficient known respiratory support is
assigned to the closest template by the weighted profile distance. Exact ties
use the frozen template order. Distance>0.75 and margin<0.01 are diagnostic flags
only and never reject an otherwise eligible assignment. Fewer than12 known-
support hours gives I in the72h primary (proportionally scaled in shorter-window
sensitivities). I remains in the main analysis under a separate baseline-hazard
stratum and represents insufficient support data, not a clinical phenotype.
This October30 policy revision replaces the prior ambiguity/distance rejection
rules; historical results remain preserved and classification-dependent models
are rerun. Templates/local clusters are unchanged in the assignment-only update. Reconstruction and six-hour
bins lose information; these are not exact pooled OM/Ward clusters or observed
patient fractions. Local and shared cluster-count reproducibility does not
establish stability of every group or biological endotypes.

Full and leave-one-site-out discovery are separately frozen. A held-out site
must not contribute to templates, classifier fitting or predictor coefficients.
This is internal-external development after the data have been examined, not
untouched prospective validation. Central run scripts currently orchestrate the
three named development sites; site-local capture itself accepts a configured
site name without hardcoded private paths.

## Outcomes and adjustment

Fit separate cause-specific models for hospital death, hospice and death-or-
hospice. Follow-up begins at ICU+72h and is capped at28 days afterward. Other
terminal events end cause-specific risk time. Export competing-event
Aalen–Johansen summaries and site absolute risks. The composite must never be
called mortality; subsequent post-discharge deaths are not observed.

Core covariates are age per10y, sex indicator, Black and other/unknown race
indicators, admission year, noncancer Charlson per5 points, POA advanced-cancer
proxy, and admission non-full/unknown code-status indicators. Patient-clustered
uncertainty and a separate I stratum retain all primary stays. Exact formulas,
actual fitted N/events, excluded cases, failed fits, and PH diagnostics are
captured locally; released model metadata are masked or banded. Measured SOFA
is a separate complete-case sensitivity, not in the core association estimates.

Advanced cancer means POA metastasis, malignant pleural effusion, or cachexia.
It is not stage. Cancer and metastatic-cancer components are removed from
Charlson. Report proxy components separately. These are retrospective billing
features, and the current transported predictor also uses the implemented
baseline feature matrix. Do not describe it as a fully prospective bedside
predictor with all covariates known at ICU admission. Associations are not
causal. Stage/ECOG/therapy intent and actual preferences may confound them.

Code status is reported at admission and landmark, with transitions/new orders,
source vocabulary, and full-code sensitivities. No recorded comfort/DNI category
is not absence of preference. Window-end code adjustment is explicitly a
sensitivity, not admission confounding control. Optional oncology fields are
descriptive until their timing/provenance is confirmed. No cancer field is
inferred to fill unavailable stage, molecular subtype, performance or therapy.

Site effects are pooled with REML and modified Hartung–Knapp intervals, tau2,
I2 and prediction intervals. Export contributor counts and withheld/failed
statuses. No missing coefficient becomes zero for association pooling.
Historical site counts are never rescaled; blanks are never interpreted as zero.
No total ICU LOS regression is fit. Validated ICU-free days remain unavailable.

## Direct severity comparison and completeness supplement

The frozen comparison uses first24h measured-component modified SOFA: measured
P/F, conventional measured MAP/standard vasoactive thresholds, platelets,
bilirubin, creatinine, and GCS. Renal scoring is creatinine-only; urine output is
not used. Missing domains remain missing. A separate zero-domain/missing-
indicator sensitivity retains all primary cases; it is not standard SOFA and
never changes clustering physiology. APACHE II/IV are not fabricated when
recorded scores/required context are unavailable or unverified.

Nested held-out comparisons include score alone, baseline, baseline+phenotype,
baseline+score, baseline+respiratory summaries, and their phenotype additions.
Harder summaries include the last completed hour's support, modality durations,
ARF detection, and support/physiology missingness. Evaluation patients are
identical across models within score cohort and endpoint. Export AUC/Brier,
calibration, paired increments and exploratory decision curves for death and
composite. Coefficient-wise two-stage ridge synthesis is exploratory;
unsupported predictive terms use a flagged zero prior. Patient-bootstrap
uncertainty conditions on frozen learners and does not propagate all discovery
uncertainty. No consistent superiority over severity/simple summaries is claimed.

Every site captures component/total completeness, joint missingness patterns,
source lab vocabulary, P/F evidence availability, GCS component/text audit,
vasoactive units/mapping flags, and recent-weight availability. Recovery candidates
include calculated MAP from unambiguous simultaneous SBP/DBP and explicit source
arterial labels/panel FiO2. Recent measured pre-ICU weight and SpO2-estimated PaO2
are separate sensitivities. SpO2 estimation is restricted to50<=SpO2<97 with
contemporaneous raw FiO2/documented room air; no missing domain becomes zero.
These are diagnostics only: frozen scores/models are NOT silently replaced.
Adoption requires new model versions and held-out evaluation.

## Sensitivities and terminal events

Eight inferential variants repeat discovery, stability, consensus, frozen local
outcomes, pooling and held-out evaluation:24h,48h, first-detected ARF, support-
only, low-flow FiO2/SF,0h and12h evidence persistence, and raw-support/no-waterfall.
Earlier windows use twelve relative-time bins and proportionally scaled known-
support thresholds. Different populations/labels/k are not assumed equivalent.
First detection is documentation-based, not biological onset. Its shifted
origin is followed by a full72h window and its own landmark.

A separate full-cohort absorbing analysis includes death, hospice and live
hospital discharge from the event-containing hour through72h. Terminal
transitions have no clinical escalation rank. D groups are descriptive and are
not regressed on the outcomes encoded in their exposures. Early events are
reported by those groups/local clusters, not imputed primary survivor P groups.
Other frozen-classifier sensitivities include first hospitalization, admission
and strict-through-landmark full code, POA U07.1 exclusion, and ARF/high support.

## Site returns and manuscript work

`share/supplement/` contains aggregate diagnostics with disclosure masking;
exact values and patient-level candidate scores remain in `local/`. The complete
federated return adds allowlisted discovery/model/validation/sensitivity files,
with frozen consensus IDs, source versions and an inventory. Sites provide or
explicitly leave pending IRB/waiver details, release/ETL/timestamp provenance,
order semantics, oncology timing, funding/conflicts and availability statements.
A completion manifest records capture, not institutional disclosure approval.

The lead team must still complete manuscript reframing, the point-by-point
response, clean/tracked versions, RECORD/STROBE, source limitations and literature
discussion. The decision letter specifies3000 words,250-word abstract, five main
figures/tables combined, and50 references; additional results belong in the
supplement. `docs/site_return_requirements.csv` maps concerns to site/coordinator
outputs, with non-computable analyses and author work identified explicitly.

### Additional model diagnostics in the site supplement

Global and term-level proportional-hazards tests are returned. Older fits are
reconstructed locally for diagnostics only, with coefficient equality required.
NU's core composite diagnostic flags possible time-varying phenotype effects.
Consequently all sites additionally return exploratory posthoc0–7 and7–28-day
cause-specific estimates and risk-set/event metadata. These do not replace the
frozen primary estimates and are not a formal demonstration that PH holds or a
prespecified confirmatory interaction analysis. Time-variation and the
interpretation of a single pooled HR remain matters for manuscript assessment.

### Implemented adapted APACHE comparison — October30

A separate corrected CLIF-adapted APACHE-II APS+age sensitivity now runs alongside
SOFA, with consistent Hb-derived hematocrit across sites, first-ICU-24h windows,
explicit unit/pairing/boundary rules and missing-component indicators. It omits
unverified chronic-health points/renal doubling and is not standard APACHE II
or APACHE IV. The original full-score-unavailability statements refer to those
standard scores. Complete-score and all-case analyses remain distinct; a matched
cohort supplies same-patient APACHE-proxy/SOFA comparisons. All original primary
artifacts remain frozen. Detailed methods/results: `apache_adapted_methods.md`.


The main descriptive time-to-event figure pairs cause-specific hospital-death
Kaplan–Meier survival with cumulative incidence of live non-hospice discharge.
The latter uses Aalen–Johansen with hospital death and hospice as competing
outcomes; it is not the complement of death-or-hospice-free KM. Both estimates
use exact pooled site risk/event counts, with 200 patient-cluster resamples
within each site for pointwise intervals. Unknown dispositions are censored at
observed exit. Cohort-wide curves start at ICU admission; phenotype curves
start at the 72-hour landmark. No postdischarge survival is inferred.


The adopted post-landmark association set includes four cause-specific Cox
endpoints: hospital death, hospice discharge, death-or-hospice, and live
non-hospice discharge. Each uses the same 28-day post-landmark clock and the
same cohort/covariates within each model specification. Competing exits are
censored at observed exit; patients with no terminal event by28days are
administratively censored. Discharge HR>1 means a faster live-discharge rate
among patients still hospitalized, not improved all-cause survival or an
absolute discharge probability. Site clustered SEs, nonestimability, robust
covariances and PH diagnostics are returned for all endpoints; pooled effects
use REML and modified Hartung-Knapp. Exploratory0-7/7-28day post-landmark effects
accompany PH assessment. The discharge P4 term shows time variation in all three
current sites, so single-HR interpretation requires the interval sensitivity.
The updated association-only extension preserves all frozen assignments and
original prediction models; prior death/hospice/composite coefficients reproduce.


## Portable two stage deployment

The external-site entry point is scripts/site.R, with the pinned environment
restored by scripts/setup.R from renv.lock. Stage 1 stops after local discovery,
stability, provisional descriptors and source/missingness diagnostics. The
coordinator alone reviews returns and issues frozen Stage 2 parameters using
scripts/coordinator.R. Reusable discovery/model workers accept a supplied site
list rather than only the three development-site names. These deployment
changes do not change the scientific exposure, assignment rule or estimands.
Stage 2 includes a coordinator exchange for frozen other-site prediction and
a common survival grid; these cannot be computed from Stage 1 fingerprints alone.
See the code README and docs/two_stage_run.md for exact commands and returns.
