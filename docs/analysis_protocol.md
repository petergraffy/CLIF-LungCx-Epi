# CCX revision analysis protocol

This protocol describes the new implementation. Parameters are proposed analysis choices to freeze after a clinical pilot, not findings from the new cohort. No historical cluster labels, counts, or site-specific repairs are reused.

## Analysis hierarchy

The sole primary analysis (`h72_primary`) clusters the corrected joint support × hypoxemic/hypercapnic ARF states in the ICU-admission +72-hour landmark cohort. Hospital death after that landmark is the primary inferential endpoint; hospice and death-or-hospice are separate secondary endpoints. Outcome events do not enter the primary clustering inputs.

Prespecified sensitivity derivations are full-cohort absorbing terminal states, support-only clustering, first-detected-ARF anchoring, ICU 24- and 48-hour landmarks, estimated low-flow FiO₂, and ARF persistence of zero or twelve hours. The absorbing analysis is descriptive only. All other sensitivity outcome models begin after their respective landmarks. `analysis_plan.csv` and run manifests record each analysis role, origin, inputs, and outcome-modeling status. Predictive validation and patient bootstrap stability target the primary analysis.

### First-detected-ARF sensitivity

Use the first qualifying positive physiology evidence detected during the first 72 ICU hours, with a six-hour pre-ICU lookback bounded by hospital admission. Detection occurs only when both members of a required pair are available; a later FiO₂ or pH measurement cannot assign an earlier detection time. Recent pre-ICU evidence anchors at ICU arrival. Evidence first detected after ICU arrival shifts the trajectory origin to detection, followed by a full 72-hour trajectory and an alive/in-hospital landmark at detection +72 hours. Patients without detected positive ARF have no onset anchor; they remain in the primary cohort when otherwise eligible. Report no evaluable measurements separately from observed nonqualifying evidence.

This is an observed detection-time sensitivity, not proven biological onset or incident ARF. Measurements and testing frequency can delay detection. The measured-FiO₂ ARF definition is used; estimated low-flow FiO₂ is a separate sensitivity. ICU-admission code status and the prepared ICU baseline severity remain baseline covariates even when the trajectory origin shifts. The landmark-flow audit reports early events relative to the shifted window. Missing support remains unknown and never excludes a hospitalization.

## Cohort and time axes

Include hospitalizations starting in `[start_date,end_date)` among adults with a primary lung cancer ICD-10-CM C34 or ICD-9-CM 162 diagnosis marked present on admission and an ICU ADT segment. Personal-history and secondary lung metastasis codes alone do not qualify. The unit is hospitalization; repeat admissions are retained, counted, grouped during resampling, and handled with patient-clustered Cox standard errors. A first-eligible-hospitalization sensitivity is included. This does not identify a patient's first lifetime hospitalization.

Geography, pollution linkage, and complete sex/race fields are not inclusion requirements. Covariate missingness is accounted for at model fitting. Cancer inclusion and POA covariates rely on retrospective coding and may misclassify active disease or baseline information availability.

First ICU arrival is `t0`. The primary derivation cohort is alive and in hospital at `t0 + 72 hours`. ICU-to-ward transfers remain eligible and use hospital-wide CLIF respiratory observations. The 24- and 48-hour analyses each have their own eligible cohort, phenotype derivation, and follow-up landmark. Their performance measures describe different populations and cannot be interpreted as a paired comparison without further restriction.

Known final hospital disposition and a valid event time are required in this implementation. Unknown disposition, missing discharge, and invalid timestamps are explicitly audited and excluded; there is no fabricated follow-up for still-admitted encounters. Recognized live dispositions are explicitly mapped in `cohort.R`; `Other` is deliberately unknown pending site adjudication. Review the disposition audit before release.

For hospital deaths, use the patient death time only if it falls between ICU arrival and hospital discharge; otherwise use discharge time as a documented proxy. Hospice and other live discharge always use hospital discharge time. Last vital signs and subsequent outpatient death never define these events. Sites must inspect discharge-proxy timing and ambiguous dispositions.

## Respiratory sequences

Each window contains exactly 24, 48, or 72 bins `[h,h+1)`. The highest recorded support tier in a bin is used; ties use the latest record. In empty bins, carry the latest prior known support for at most six hours. A pre-ICU observation may initialize the first bin only if it is inside the same hospitalization and carry limit. Never backward-fill from later observations. Unknown device values stay unknown. No hospitalization is excluded for missing support or insufficient support coverage; entirely unknown support sequences remain eligible. Support coverage and imputation are reported descriptively.

Support states are room air, low-flow oxygen, HFNC, CPAP, NIPPV, IMV, other, and unknown. “Other” has no validated severity interpretation and is excluded from reference-ranking calculations.

The physiology state preserves each component separately:

- Hypoxemia positive `H1`, nonqualifying observed evidence `H0`, or unknown `HU`.
- Hypercapnia/acidemia positive `C1`, nonqualifying paired evidence `C0`, or unknown `CU`.

The primary joint alphabet is all support-by-physiology combinations, plus terminal states only in the absorbing sensitivity. `H0_CU` does not establish normal blood gases or rule out respiratory failure. No physiology state is labeled “normal ABG.”

Primary hypoxemia evidence is paired arterial PaO₂/FiO₂ ≤300, or SpO₂ <90% on documented room air. Room air is assigned FiO₂ 0.21 only from an explicit room-air support record; other devices require measured FiO₂. FiO₂ fractions and percentages are normalized and values outside 0.21–1 are rejected. Nearest FiO₂ pairing allows at most one hour, including either direction within the completed trajectory window. Hypercapnic evidence requires arterial PaCO₂ ≥45 and paired arterial pH <7.35 within two hours. Results use CLIF reference units. A normal result for one component does not resolve missingness in the other.

Recent positive or nonqualifying evidence contributes for up to six hours; any positive evidence within the permitted lookback takes precedence. After expiry, the component becomes unknown. Sensitivities use zero- and twelve-hour persistence. The low-flow sensitivity estimates FiO₂ only for low-flow observations with flow 1–6 L/min using `min(0.45,0.21 + 0.04 × flow)` and applies SpO₂/FiO₂ ≤315 for saturation ≤97%. These approximations are sensitivity assumptions, not measured FiO₂. Pulse oximetry measurement bias, including differential bias by skin pigmentation, remains a limitation.

No support-coverage threshold is applied. The primary analysis retains every landmark-eligible hospitalization, and the absorbing sensitivity retains every cohort hospitalization with a valid terminal outcome/time, even if no pre-event hourly bins exist. Terminal states begin in the bin containing the event, discretizing timing by up to one hour. They are absorbing through the window. This analysis describes early outcomes; it is excluded from subsequent phenotype-outcome models. Unknown support may form its own clusters, which must be interpreted as documentation patterns rather than established clinical phenotypes.

## Derivation and harmonization

Local derivation uses optimal matching, transition-rate substitution costs, unit insertion/deletion costs, Ward D2 linkage, and the highest mean silhouette among k=2 through the configured maximum, capped at N−1. Ties select the smaller k. Four clusters are not enforced. Patient bootstrap stability resamples entire sets of eligible hospitalizations belonging to a patient, reclusters using the original dissimilarities, assigns original sequences to bootstrap medoids, and reports matched Jaccard and adjusted Rand indices. This assesses stability conditional on the original substitution-cost estimate and selected k.

Central harmonization compares hourly state probability profiles with root-mean-square distances over available common features. Whole hourly distributions are suppressed when necessary; missing features are not set to zero. At least 50% common features are required by default. All possible prototypes are considered in rectangular one-to-one assignment. Extra clusters or insufficiently overlapping profiles can remain unmatched. A seed site is chosen deterministically using profile availability, then site name. Prototypes are updated using unsuppressed cluster-size weights until assignments stabilize; failure to converge stops consensus generation.

The default consensus k is the integer median of local k; it may be explicitly specified. Phenotypes P1…Pk are ordered only by average respiratory support intensity among known ranked support states. This ordering chooses a consistent reference, not a validated ordinal disease scale. No outcome, vasopressor, SOFA, or LOS information enters matching or ordering. The absorbing sensitivity is an explicit exception because terminal outcomes are its sequence inputs. Mapping and matching distances are exported, and leaving out each site tests consensus sensitivity. This is a reproducibility check, not external validation of a prognostic model.

## Clinical characterization and outcomes

Code status uses documented orders in the current admission, with an optional explicit preadmission lookback. Unknown stays unknown. DNR without documented DNI is separate from DNI. Report admission status, status at window completion, and new DNI/comfort orders. Changes within the window are descriptive; they are not treated as baseline confounders. These orders do not fully capture goals-of-care discussions or actual treatment ceilings.

Noncancer Charlson removes the cancer and metastatic-cancer domains and applies severity hierarchy within other domains. Advanced cancer means recorded POA metastatic disease, malignant pleural effusion, or cachexia; its components are reported separately. This is a proxy, not stage. C80/199 unspecified malignancy and generic pleural effusion J91 are not automatically treated as metastatic disease or malignant effusion. Optional `clif_oncology`/`oncology` input may supply one row per hospitalization with histology, stage, ECOG, therapy line, and therapy intent. Availability does not establish baseline timing; these fields are descriptive until independently validated.

The implemented severity comparator is a **modified SOFA proxy**: measured PaO₂/FiO₂ paired in time, conventional thresholds for measured MAP/standard SOFA vasoactives, platelets, bilirubin, creatinine, and GCS. Renal scoring does not include urine output. Missing domains stay missing. No weight is invented for dose conversion and no saturation-derived PaO₂ is used. Unconvertible standard vasoactive doses make the cardiovascular domain unknown. This cannot be described as complete standard SOFA or APACHE IV. Sedation-related GCS and nonstandard vasoactive agents require clinical review. A separately labeled missing-domain-zero sensitivity is run for both associations and predictive validation; it never silently replaces the primary score. All-missing severity remains missing even in this sensitivity.

Primary inferential estimates are post-landmark cause-specific hazard ratios for hospital death. Hospice discharge and death-or-hospice are separate endpoints. Other terminal hospital events end cause-specific risk time. Aalen–Johansen curves treat death, hospice, and other live discharge as competing events, using hospitalization as the episode identifier and patient-clustered uncertainty for repeat admissions. Whole-stay outcomes are descriptive in phenotype tables. No total ICU LOS regression is fit. Actual fitted N, patient N, events, exclusions, formulas, and proportional-hazards diagnostics accompany each model.

Adjusted associations include measured severity, demographics, year, noncancer Charlson, cancer proxy, and admission code status. Constant covariates are omitted and the exact formula is exported. Models are retrospective associations, not causal treatment effects. Sensitivities use first hospitalization, admission full code, no POA U07.1 diagnosis, and measured ARF or high support. POA U07.1 exclusion does not exclude hospital-acquired COVID or all pandemic effects.

## Prognostic comparison

The endpoint is **in-hospital death within 28 days after the landmark**, not all-cause 28-day mortality. Earlier hospice/other live discharge is a competing outcome. All included encounters have known final hospital outcomes; therefore this fixed-horizon binary endpoint has no administratively censored observations. Discharged patients' subsequent deaths are not observed or imputed.

Models share the same complete-case sample:

1. Age, sex, race, admission year, measured modified SOFA, and admission code status.
2. Model 1 plus last completed hourly support state and cumulative hours in each nonreference modality, including unknown hours.
3. Model 2 plus phenotype.

Diagnosis-based comorbidity/cancer features are excluded from this prospective-information comparator because discharge billing diagnoses are not reliably known at the landmark. They remain in retrospective association models. Predictive validation groups all admissions of a patient in the same fold, and relearns transition costs, selected k, clusters, medoids, and test assignments using the training fold only. Performance is for the local phenotype-learning procedure, not a centrally trained deployable classifier. The full-data consensus mapping is not used to learn validation phenotypes.

Exports include out-of-fold AUC, Brier score, calibration intercept and slope, paired model increments, calibration bins, and decision curves. Patient bootstrap intervals are conditional on the out-of-fold predictions; they do not rerun the entire cross-validation procedure. Fold failures, separation, and unseen factor levels are reported; incomplete validation yields no performance claim. Decision curves describe hypothetical risk thresholds and do not validate treatment or goals-of-care decisions.

## Pooling and additional work

Pool identical phenotype-vs-P1 cause-specific log hazard ratios using REML heterogeneity estimation and modified Hartung–Knapp inference. Export site effects, contributing site count, τ², I², confidence intervals, and prediction intervals. One-site effects are not labeled pooled estimates. Never combine different consensus versions, references, or source versions. Site validation results remain site-specific; they are not presented as pooled external validation.

Complete 28-day ICU-free days require a validated observation scheme for subsequent ICU stays and deaths. The code deliberately does not invent post-discharge ICU-free time. A common-population paired comparison of early landmarks, a deployable external prediction model, and oncology-enhanced adjustment are further analyses to specify if data permit. Confirm IRB/waiver details, data availability, conflicts, and RECORD/STROBE reporting outside the analysis code.

Method references: [CLIF 2.1 dictionary](https://clif-icu.com/data-dictionary/data-dictionary-2.1.0), [survival competing-risk methods](https://stat.ethz.ch/R-manual/R-devel/RHOME/library/survival/doc/compete.pdf), [metafor inference](https://wviechtb.github.io/metafor/reference/rma.uni.html), [prediction intervals](https://wviechtb.github.io/metafor/reference/predict.rma.html).
