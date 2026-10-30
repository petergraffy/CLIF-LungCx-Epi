from pathlib import Path
import json,csv,os
r=Path(os.environ.get('CLIF_MIMIC_OUTPUT','output/run_MIMIC_external_20261023'));out=r/'reports';out.mkdir(exist_ok=True)
def read(name):
 with (r/'aggregate'/name).open() as f:return list(csv.DictReader(f))
c=read('cohort_counts.csv')[0];q=read('panel_quality.csv')[0];metrics=read('transport_metrics.csv');inc=read('paired_increments.csv')
lines=['MIMIC CLIF: restricted lung-cancer external transportability assessment','2026-10-23','',
'LUNG-CANCER RESTRICTION',
'All included hospitalizations have an ICD-10 C34.x or ICD-9 162.x code in that hospitalization. Secondary malignant neoplasm of lung codes (C78.0/197.0) and personal-history codes alone do not qualify. This matches the primary-organ code rule of the three-site analysis. It does not require lung cancer to be the principal reason for admission or confirm histology/stage.',
f"Adult ICU coded-lung-cancer cohort: {c['eligible_hospitalizations']} hospitalizations, {c['patients']} patients. Alive/in-hospital at first ICU arrival +72h: {c['landmark_hospitalizations']} hospitalizations ({q['patients']} patients).",
f"A prior completed lung-cancer-coded hospitalization is documented for {c['previous_lung_cancer_hospitalizations']} eligible hospitalizations; {q['prior_completed_lung_cancer_encounter_n']} remain at the +72h landmark. This stricter subgroup requires the earlier hospitalization to have ended by the current admission and still requires a current lung-cancer code. It establishes prior encounter documentation, not exact clinical diagnosis time or active stage.",
'POA flags are all null. No observed POA evidence was invented. A temporary eligibility-only adapter lets the existing code select current lung-cancer-coded encounters; the source diagnosis table remains unchanged, and discharge diagnoses are not used as baseline comorbidity/cancer predictors. The stricter documented-history subgroup is the strongest available option using these CLIF tables. Neither cohort is identical to the POA-confirmed primary population.',
'Shifted source dates preserve elapsed time and within-patient chronology. The real 2018–2024 filter is not applied to deidentified calendar years, and shifted year is not a predictor. All source periods are included; actual calendar-era matching cannot be recovered from the supplied CLIF patient table.', '',
'PRIMARY COMPATIBILITY FINDING',
'The supplied source has no room-air device category or device name in the selected cancer cohort. The processed panel also has zero room-air hours. Accordingly, zero hospitalizations map to the frozen P1 room-air recovery phenotype. This is not evidence that no lung-cancer ICU patient recovered to room air. It is a respiratory-documentation/ETL coverage limitation that prevents complete four-phenotype validation and prevents P2–P4 versus P1 association estimates.',
'All source device labels and the room-air audit are provided. Missing device/support is not converted to room air. Resolving this limitation requires original MIMIC oxygen-delivery documentation or a source conversion that preserves explicit room-air/off-oxygen observations, then a new versioned external analysis. No such observations can be recovered from all-null device rows by assumption.',
f"Waterfall changes unknown-support hours from {float(q['raw_unknown_pct']):.1f}% to {float(q['waterfall_unknown_pct']):.1f}%. Entirely unknown support remains in {q['entirely_unknown_stays']} stays. All {q['n']} landmark hospitalizations are retained; I is an information-status assignment, not an excluded patient or biological phenotype.",
'Unknown-support rates by hour-midpoint ADT location are reported separately. The window follows hospitalizations after ICU arrival, including subsequent ward/stepdown time; unknown support does not imply that the patient remained in the ICU.', '',
'FROZEN EXTERNAL TEST',
'The original UCMC/NU/Rush phenotype templates were unchanged. MIMIC did not participate in shared template learning, cluster-count selection or predictive coefficient fitting. Local MIMIC clustering was performed only as a separate descriptive diagnostic after frozen prediction evaluation; it did not modify the shared phenotypes.',
'Portable predictors were newly trained at the original three sites, without calendar year, POA Charlson or POA advanced-cancer covariates that cannot be measured comparably in MIMIC. Baseline predictors are age, sex, race grouping and admission code status. Respiratory-summary and phenotype additions follow the prior analyses, including hour-72 support and separate cumulative modality durations. Coefficientwise inverse-variance synthesis retains the earlier exploratory ridge-logistic design, with unsupported coefficients flagged as zero-prior terms. This evaluates a new portable predictor, not the earlier full-covariate risk equation.',
'There is no MIMIC refitting or recalibration. Endpoint: hospital death/hospice within 28 days after ICU +72h, with separate hospital-death performance also reported. Patient bootstrap intervals (200 replicates) condition on frozen templates and coefficients; they do not propagate every training/discovery uncertainty.',
'Composite AUC | baseline | +phenotype | +summaries | +summaries+phenotype']
for sub in ['coded_lung_cancer','prior_completed_lung_cancer_encounter']:
 vals=[]
 for model in ['baseline','phenotype','summaries','summaries_phenotype']:
  z=next(z for z in metrics if z['subgroup']==sub and z['score_cohort']=='all' and z['endpoint']=='death_or_hospice' and z['model']==model and z['metric']=='auc')
  vals.append(f"{float(z['estimate']):.3f}")
 lines.append(sub+' | '+' | '.join(vals))
for sub in ['coded_lung_cancer','prior_completed_lung_cancer_encounter']:
 z=next(z for z in inc if z['subgroup']==sub and z['score_cohort']=='all' and z['endpoint']=='death_or_hospice' and z['comparison']=='phenotype_beyond_summaries' and z['metric']=='auc')
 lines.append(f"Adding phenotype beyond summaries, {sub}: delta AUC {float(z['estimate']):+.4f} (95% CI {float(z['conf_low']):+.4f} to {float(z['conf_high']):+.4f}).")
lines+=['Discrimination is modest. Phenotypes improve on the limited baseline but do not show consistent incremental benefit beyond simple respiratory summaries. Results do not establish successful transport of all four clinical patterns. The absent room-air phenotype, different cohort ascertainment, differing source periods and sparse source observations limit interpretation.', '',
'CLINICAL AND SEVERITY SOURCES',
f"The source provides explicitly named DNR, DNI and comfort-measures-only orders. Admission/window-end status, newly documented changes and missingness are described by frozen phenotype; these remain charted documentation, not complete ascertainment of care preferences. First-24h complete modified SOFA is available for {q['measured_sofa_n']} landmark stays, using raw GCS components measured at the same timestamp. The gcs_eyes category is normalized to gcs_eye before forming complete triplets; incomplete/intubated verbal measurements are not defaulted to normal. Default/imputed GCS totals from the standard assessment table are excluded. Renal SOFA uses creatinine, not urine output, as in the earlier modified score.",
'Separate complete-component severity comparisons are included, with their actual denominators. They are small selected samples, not a standard-score benchmark. APACHE II/IV is not supplied as a validated score and has not been manufactured from incomplete CLIF fields.',
'Adjusted association attempts use only comparable demographic/code-status variables. All P2–P4 versus P1 contrasts are explicitly unavailable because P1 is absent. MIMIC estimates are not added to the POA-confirmed three-site meta-analysis. No causal effect, therapy threshold or superiority over established severity scores is claimed.', '',
'COHORT FLOW',
*['  '+z['landmark_reason']+': '+z['n'] for z in read('landmark_flow.csv')],
'Unclassifiable dispositions remain unknown; no death/live-discharge category is guessed. Early terminal events do not enter the primary landmark exposure/outcome association.', '',
'FILES AND NEXT DATA REQUIREMENT',
'aggregate/cohort_counts.csv and compatibility_manifest.json: exact cohort definition and source adaptations.',
'aggregate/device_catalogue.csv, room_air_compatibility.json, unknown_by_location.csv: respiratory-source limitation.',
'aggregate/transport_metrics.csv, paired_increments.csv, model_denominators.csv: frozen external performance, uncertainty and actual N.',
'aggregate/associations.csv: explicit absent-reference status; no primary meta-analysis pooling.',
'local_discovery/sites/MIMIC/export/discovery/: interval fingerprints, local silhouette/stability bands; descriptive only.',
'aggregate/final_audit.json: eligibility/history/retention/frozen-prediction checks.',
'Patient identifiers, sequences, timestamps and individual predictions remain under local/. Aggregates are local rehearsal outputs and require disclosure review before external exchange. The source archive and original three-site outputs were not modified.', '',
'PRIMARY REFERENCES',
'MIMIC CLIF v1.1.0 (CLIF 2.1): https://physionet.org/content/mimic-iv-ext-clif/1.1.0/',
'POA/GCS conversion notes: https://github.com/Common-Longitudinal-ICU-data-Format/CLIF-MIMIC/blob/main/ISSUESLOG.md',
'MIMIC date/anchor documentation: https://physionet.org/content/mimiciv/3.1/']
(out/'mimic_lung_cancer_report.txt').write_text('\n'.join(lines)+'\n')
print(out/'mimic_lung_cancer_report.txt')
