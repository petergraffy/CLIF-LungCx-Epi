from pathlib import Path
import csv, os, json
root=Path(os.environ.get('CLIF_REVISION_OUTPUT','output/run_revision_gaps_20261022'))
out=root/'reports'; tables=out/'tables';out.mkdir(parents=True,exist_ok=True)
def read(name):
 with (tables/name).open() as f:return list(csv.DictReader(f))
def num(x):
 try:return float(x)
 except:return float('nan')
def fmt(x):return f'{num(x):.3f}' if x else 'withheld'
metrics=read('score_transport_metrics.csv');increments=read('score_paired_increments.csv')
lines=['Reviewer-gap analyses: three-site federated rehearsal','Completed 2026-10-22; revised analyses, not a submission-ready ten-site manuscript.','',
'INTERPRETATION',
'The corrected primary remains the joint support–ARF ICU +72-hour landmark cohort (UCMC 650, NU 1,691, Rush 969; total 3,310). Its frozen phenotype definitions and original outputs are preserved. New work adds the requested severity comparison, harder summary comparator, ARF characterization, code-status analyses, separately learned sensitivity phenotypes, and reporting materials.',
'The results support descriptive respiratory trajectories and associations with subsequent hospital outcomes. They do not establish general prognostic superiority, biologic ARF subtypes, independence from treatment preferences, or causal treatment effects. All findings are exploratory development results from datasets already examined.', '',
'DIRECT SOFA COMPARISON',
'We compared eight transported models: baseline; baseline+phenotype; modified SOFA alone; baseline+modified SOFA; baseline+simple summaries; baseline+summaries+phenotype; baseline+modified SOFA+summaries; and that final model+phenotype. Summaries now explicitly include support in the last pre-landmark hour, separate CPAP/NIPPV/HFNC/IMV durations, low-flow duration, any ARF detection, and support/physiology unknown fractions.',
'Every comparison uses exactly the same evaluation hospitalizations within its cohort and endpoint. Complete measured-component SOFA is available for only 99 UCMC, 219 NU and 62 Rush hospitalizations (380/3,310; 11.5%). All six domains must be present. This is a modified SOFA: creatinine-only renal scoring, no urine-output criterion, and the existing measured first-24-hour respiratory/vasopressor/GCS implementation. It must not be called a fully validated standard SOFA score.',
'A separate missing-domain-zero sensitivity retains all 3,310 cases. Missing components contribute zero to this sensitivity index, and adjusted score models include six domain-missing indicators. The score-alone comparator has only the proxy index. This is not a missing-as-normal primary analysis, a standard SOFA score, or a substitute for recovering missing measurements.',
'Score-based risks are ridge-logistic transformations learned at the other two sites, not published SOFA mortality probabilities. Predictors and phenotype templates exclude the held-out site; no held-out fitting or recalibration occurs. Reported confidence intervals use 200 patient bootstraps conditional on the frozen learner. Calibration intercept/slope are point estimates in this added comparison; AUC, Brier and decision-curve estimates have bootstrap intervals. Coefficientwise inverse-variance synthesis remains exploratory, with unsupported coefficients explicitly flagged as zero-prior terms.',
'', 'Complete measured-component cohort: 28-day hospital death/hospice AUC',
'Site | modified SOFA alone | baseline+modified SOFA | baseline+phenotype | baseline+SOFA+summaries | +phenotype']
for site in ['UCMC','NU','RUSH']:
 vals=[]
 for model in ['sofa_only','sofa','phenotype','sofa_summaries','sofa_summaries_phenotype']:
  r=next(r for r in metrics if r['site']==site and r['cohort']=='measured_sofa' and r['endpoint']=='death_or_hospice' and r['metric']=='auc' and r['model']==model)
  vals.append(fmt(r['estimate']))
 lines.append(site+' | '+' | '.join(vals))
lines+=['','Adding phenotype beyond baseline + modified SOFA + simple summaries (paired AUC difference):']
for cohort in ['measured_sofa','missing_domain_zero_sensitivity']:
 for site in ['UCMC','NU','RUSH']:
  r=next(r for r in increments if r['site']==site and r['cohort']==cohort and r['endpoint']=='death_or_hospice' and r['metric']=='auc' and r['comparison']=='phenotype_beyond_sofa')
  lines.append(f"  {cohort}, {site}: {fmt(r['estimate'])} (95% CI {fmt(r['conf_low'])} to {fmt(r['conf_high'])})")
lines+=['The complete-case incremental gains are small and uncertain at UCMC and NU, and negative at Rush. The all-case sensitivity also shows no consistent incremental advantage. Higher AUC for one model at one site is not evidence of universal superiority. Decision curves are hypothetical risk-threshold assessments, not evidence that a treatment or goals-of-care intervention improves outcomes.',
'Death-only prediction is also evaluated separately. In the complete-case SOFA cohort, Rush performance is withheld because its endpoint partition is too small. All model denominators/statuses and suppression are explicit; no unavailable performance is assigned zero.', '',
'APACHE II AND IV: UNRESOLVED DATA REQUIREMENT',
'Neither a recorded APACHE score nor a validated APACHE input package exists in these folders. UCMC has no hematocrit under canonical categories, audited source labels, or the queried hematocrit LOINCs. NU has some arterial hematocrit recorded under lab_category=other; its availability is audited separately instead of being confused with absent source data. Rush has canonical hematocrit. The extracts do not establish the required chronic-health criteria, surgical admission context, or APACHE-II acute-renal-failure modifier; APACHE IV requires additional admission/history information. Oncology stage or an ICD metastasis proxy cannot stand in for those fields.',
'Therefore no valid direct APACHE II/IV comparison can be calculated from the current data. apache_comparison_status.csv contains explicit unavailable rows, not invented scores or AUCs. An incomplete APS or hemoglobin-derived hematocrit proxy has deliberately not been labeled APACHE. This requested comparison remains open unless additional validated score inputs become available. The narrative must use the editor’s descriptive-phenotyping option.', '',
'ARF TYPE, TIMING, PERSISTENCE AND CONTRIBUTION',
'New tables describe hypoxemic, hypercapnic/acidemic, simultaneous mixed, both-detected-at-different-times, and no-positive-evidence patterns within each frozen phenotype; first detected window hour and positivity in hours 48–71; and six-hour modal evidence-state profiles. H1_CU remains distinct from hypoxemia with negative hypercapnic evidence. No positive evidence does not mean no ARF. First detection is documentation-based and is not biological onset.',
'Hourly positive fractions represent observed-evidence rules plus the specified carry interval. They are not uninterrupted confirmed biological failure. The mixed-hour fraction overlaps the hypoxemic and hypercapnic fractions.',
'Adjusted models add hypoxemic, hypercapnic and simultaneous-mixed fractions jointly, scaled to 24 hours of positive evidence. Because channels overlap, the mixed coefficient is a conditional additional term, not the standalone prognosis of a mixed-ARF patient. Effects are imprecise. For the composite, pooled HRs are approximately 1.36 (0.82–2.26) for hypoxemic evidence, 0.99 (0.20–4.79) for hypercapnic evidence, and 0.83 (0.09–7.28) for the additional mixed term. None establishes additional independent prognostic information.',
'The support-only federated sensitivity independently learns four shared groups. Partition agreement with the original joint analysis is high (ARI approximately 0.859 UCMC, 0.888 NU and 0.909 Rush, including I as an assignment category). This indicates substantial support dominance; it does not prove physiology has zero value. Raw-support versus waterfall-support comparisons are provided separately.', '',
'TREATMENT LIMITATIONS, ONCOLOGY AND OUTCOMES',
'New site tables show admission and landmark code status, transitions, newly documented DNI/comfort status, outcomes by landmark code status, and missingness by phenotype. DNR alone remains separate from DNI. A newly documented order after ICU arrival can reflect delayed documentation, especially when baseline status is unknown. No record of a treatment limitation is not proof that no limitation existed.',
'Exploratory adjustment for landmark code status and newly documented DNI/comfort leaves composite HRs versus room-air recovery at approximately 1.64 (0.87–3.10), 2.29 (1.24–4.23), and 2.57 (1.06–6.26) for P2–P4. Adjustment for pre-landmark treatment decisions is descriptive and may condition on a mediator; it is not causal control of treatment-limit bias.',
'A source-order catalogue audit found no comfort-care category in any of the three extracts. NU has no explicit DNI category; Rush maps some DNR-named orders to canonical DNR/DNI, whose intubation implications are not independently verified. Consequently zero comfort detections cannot establish absence of comfort-care decisions, and canonical DNI flags are not uniformly verified clinical DNI orders. The baseline core covariate combines nonfull categories, but detailed DNI/change analyses require these qualifications.',
'Full-code-at-admission and documented-full-with-no-category-change-through-landmark sensitivities were run. Only one site supplies releasable contrasts, so these cannot be pooled. End status alone was not treated as proof of full code throughout the window. The core associations cannot be claimed independent of care preferences.',
'No-COVID, first-hospitalization and ARF-or-high-support subgroup models were run using the frozen primary classifier. The COVID definition is POA U07.1. Excluding COVID leaves the broad composite associations similar. Restricted subgroup analyses are less precise; all model-specific denominators are available locally, with banded exports.',
'POA metastatic disease, malignant pleural effusion and cachexia are reported separately, and define the explicit advanced-cancer proxy. Cancer stage, ECOG, molecular subtype, treatment line/intent and detailed goals-of-care discussions remain unavailable. These are substantive residual-confounding limitations.',
'Death, hospice and their composite remain separate endpoints. Separate phenotype-versus-reference death associations still have insufficient releasable site estimates for meta-analysis. Prediction performance is a different estimand and does not repair this association-pooling limitation. Total ICU LOS is not restored as an inferential outcome. Validated 28-day ICU-free days would require reliable follow-up and death ascertainment beyond hospital discharge.', '',
'FEDERATED SENSITIVITIES',
'All eight inferential sensitivities repeat local common-cost sequence clustering, 100 patient bootstraps, aggregate fingerprint exchange, central discovery, four frozen consensus versions, twelve local outcome passes per sensitivity, pooling and held-out prediction. Support is bounded waterfall with original/raw-source physiology, including raw low-flow/persistence variants. The raw-support sensitivity removes waterfall. Twenty-four/48-hour comparisons use twelve 2/4-hour bins instead of twelve 6-hour bins. Known-support requirements scale to one-sixth of the window. No support-coverage cohort exclusion is applied.',
'All three sites and all model versions retain every eligible hospitalization, including I. I means insufficient, ambiguous or poor-match information, not an additional biological phenotype. Different variants have independently learned labels; P2 in one sensitivity is not automatically P2 in another. Ordinal full-versus-held-out label agreement is not interpreted as biological accuracy when selected cluster counts differ.',
'Analysis | shared k with all sites | k holding out UCMC / NU / Rush']
cons=read('sensitivity_consensus.csv')
for variant in ['h24','h48','onset','support_only','low_flow','carry0','carry12','raw_support']:
 vals=[next(r['k'] for r in cons if r['analysis']==variant and r['version']==v) for v in ['full','loo_UCMC','loo_NU','loo_RUSH']]
 lines.append(variant+' | '+vals[0]+' | '+' / '.join(vals[1:]))
lines+=['Twenty-four-hour landmark cohorts contain 795 UCMC, 2,142 NU and 1,226 Rush hospitalizations; 48-hour cohorts contain 722, 1,908 and 1,115. First-detected-ARF cohorts contain 276, 774 and 363. These are different selected populations with different exposure windows and follow-up origins; their AUCs are not paired comparisons.',
'A separate federated absorbing-state analysis (absorbing_v2) includes death, hospice and live-discharge terminal states explicitly. Its groups are named D, not primary P, and early-event counts are exported by descriptive group. No outcome regression is run on these terminal-state-containing exposures. The initial absorbing diagnostic is preserved but superseded: v2 treats transitions into terminal states as unranked/unknown in the support-change channel. Primary landmark groups contain zero within-window terminal events by design.', '',
'COHORT AND REPORTING LIMITS',
'The three NU and sixty Rush encounters with unclassifiable source outcomes remain explicitly audited, separate from unknown support. They cannot be assigned a valid death/live-discharge outcome by assumption. Support missingness does not exclude anyone. Cohort flows, continuous-variable observed denominators and missing counts, categorical count bands, IMV versus vasopressor fields, and model-specific actual local N/events are generated independently.',
'This is a three-site development analysis, not a replacement return from the original ten-site consortium. Seven additional sites would need to execute and review the final frozen workflow for a ten-site revision. No patient identifiers or individual predictions are in exchanged exports. This remains one-machine logical federation; privacy/disclosure review and real institutional execution are not completed by these diagnostics.',
'IRB study number, approval date, board, study title and consent-waiver details cannot be recovered from the submitted manuscript. Author funding/conflict statements also need confirmation. The response and manuscript text are drafts with explicit unresolved items; the original submitted Word files have not been overwritten. RECORD/STROBE tracking is supplied as a reporting crosswalk, not a falsely completed page-number checklist.', '',
'FILES',
'reports/tables/score_transport_metrics.csv, score_paired_increments.csv and score_decision_curves.csv: direct modified-SOFA and exact nested-summary comparisons.',
'reports/tables/phenotype_characteristics.csv, code_transitions.csv, outcomes_by_code.csv, arf_type_profiles.csv, pooled_arf_associations.csv: physiology and care-limitation characterization.',
'reports/tables/sensitivity_*: independent federated sensitivity fits, cohorts, cluster counts, pooled effects, profiles and transported performance.',
'reports/tables/continuous_denominators.csv and sites/SITE/local/model_denominators.csv: export-safe descriptive denominators and exact local regression denominators.',
'reports/revision_audit.json: all local retention checks, 32 frozen consensus checks and 96 outcome-version checks. Absorbing checks are separate.',
'reports/figures/: direct SOFA comparison, incremental AUC, exploratory decision curves, ARF contribution and earlier landmarks.', '',
'PRIMARY SOURCES FOR THE REVISION',
'Original SOFA: https://pubmed.ncbi.nlm.nih.gov/8844239/',
'APACHE II: https://pubmed.ncbi.nlm.nih.gov/3928249/',
'APACHE IV: https://pubmed.ncbi.nlm.nih.gov/16540951/',
'Sathe respiratory subphenotypes: https://pubmed.ncbi.nlm.nih.gov/34526076/',
'Sathe trial-enrichment prediction: https://doi.org/10.1097/CCM.0000000000006181',
'Pulse-oximetry bias: https://doi.org/10.1056/NEJMc2029240',
'STROBE: https://www.strobe-statement.org/',
'RECORD: https://doi.org/10.1371/journal.pmed.1001885']
(out/'revision_gap_report.txt').write_text('\n'.join(lines)+'\n')
print(out/'revision_gap_report.txt')
