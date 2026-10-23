"""Create the aggregate-only plain-text report for a completed rehearsal."""
import csv,json,os
from pathlib import Path
root=Path(os.environ.get('CLIF_FEDERATED_OUTPUT','output/run_federated_rehearsal_20261022'))
report=root/'reports'
def read(name):
 with (report/'tables'/name).open() as f:return list(csv.DictReader(f))
def fnum(x):return float(x) if x else float('nan')
metrics=read('transport_metrics.csv');pool=read('pooled_associations.csv');coverage=read('transport_coverage.csv');counts=read('cohort_counts.csv')
labels={'P1':'Recovery toward room air','P2':'Low-flow oxygen','P3':'HFNC / noninvasive support','P4':'Invasive ventilation','I':'Insufficient or ambiguous information'}
lines=['Federated lung-cancer support / ARF rehearsal -- 2026-10-22','',
'Completed: three local discovery workers, aggregate fingerprint exchange, shared phenotype discovery, frozen assignments, local outcome models, meta-analysis, and three site-held-out prediction tests.',
'All 3,310 primary landmark hospitalizations were retained (UCMC 650; NU 1,691; Rush 969). Patient identifiers, individual sequences, covariates, predictions and outcomes remained in site-local directories.','',
'MAIN FINDINGS',
'Four shared phenotypes were selected from aggregate fingerprints. The same number was selected independently when each site was held out. Labels reflect support patterns within the joint support-ARF model: P1 recovery toward room air; P2 low-flow oxygen; P3 HFNC/noninvasive support; P4 invasive ventilation.',
'Core-adjusted 28-day death/hospice cause-specific hazard ratios relative to P1:']
for z in pool:
 if z['model']=='core' and z['endpoint']=='death_or_hospice' and z['status']=='ok':
  lines.append(f"  {z['phenotype']} {labels[z['phenotype']]}: HR {fnum(z['hazard_ratio']):.2f}, 95% CI {fnum(z['hr_low']):.2f} to {fnum(z['hr_high']):.2f}; contributing sites={z['k']}; I-squared={fnum(z['i2_percent']):.1f}%; prediction interval {fnum(z['hr_prediction_low']):.2f} to {fnum(z['hr_prediction_high']):.2f}.")
lines+=['','HELD-OUT PREDICTION','Outcome: hospital death or hospice within 28 days after the +72-hour landmark. Each site was tested using profiles and predictor coefficients learned only at the other two sites.',
'Site | Baseline AUC | Baseline + phenotype | Baseline + simple summaries | Baseline + summaries + phenotype']
for site in ['UCMC','NU','RUSH']:
 vals=[]
 for model in ['baseline','plus_phenotype','simple_summaries','summaries_plus_phenotype']:
  z=next(x for x in metrics if x['site']==site and x['model']==model and x['metric']=='auc');vals.append(f"{fnum(z['estimate']):.3f} ({fnum(z['conf_low']):.3f}-{fnum(z['conf_high']):.3f})")
 lines.append(site+' | '+' | '.join(vals))
lines+=['','ADDING PHENOTYPE BEYOND SIMPLE SUMMARIES']
for site in ['UCMC','NU','RUSH']:
 z=next(x for x in metrics if x['site']==site and x['model']=='phenotype_beyond_summaries' and x['metric']=='auc')
 b=next(x for x in metrics if x['site']==site and x['model']=='phenotype_beyond_summaries' and x['metric']=='brier')
 lines.append(f"{site}: delta AUC {fnum(z['estimate']):+.4f} (95% CI {fnum(z['conf_low']):+.4f} to {fnum(z['conf_high']):+.4f}); delta Brier {fnum(b['estimate']):+.4f} (95% CI {fnum(b['conf_low']):+.4f} to {fnum(b['conf_high']):+.4f}). Positive delta AUC and negative delta Brier favor adding phenotype.")
lines+=['Phenotype improves discrimination and Brier score over baseline covariates at all three sites. Its incremental value beyond simple respiratory summaries is small and inconsistent: a small improvement at NU, small worsening at UCMC, and no clear improvement at Rush. This supports reproducible descriptive phenotyping and outcome association; it does not establish superior prediction over straightforward exposure summaries. The harder comparator was added as an exploratory robustness check after the initial baseline comparison, without changing the frozen phenotypes.','',
'COVERAGE AND ASSIGNMENT REPRODUCIBILITY']
for c in coverage:
 n=int(next(z['n'] for z in counts if z['site']==c['site']))
 lines.append(f"{c['site']}: {n} evaluated; {c['classified_lower']}-{c['classified_upper']} assigned a clinical phenotype in the held-out run; full-versus-held-out assignment agreement falls in the exported {100*fnum(c['agreement_lower']):.0f}-{100*fnum(c['agreement_upper']):.0f}% interval. Remaining cases retain I and remain in outcome/prediction analyses.")
lines+=['','WORKFLOW AND METHODS',
'This federated rehearsal covers the 72-hour primary workflow. Prior local 24/48-hour, absorbing-event, ARF-onset and alternative-physiology sensitivities have not all been refederated using the new classifier.',
'1. Each local worker reused its own previously validated 72-hour hourly panel: bounded CLIFpy 2.1.0 waterfall support with original raw-evidence ARF states. No raw or individual sequence data were passed to the coordinator. This rehearsal rebuilds clusters, rather than copying earlier cluster assignments.',
'2. Local optimal-matching sequence clustering used a common prespecified constant substitution cost of 2, indel cost 1, and Ward.D2 linkage. This deliberately differs from earlier site-estimated TRATE costs. A learned shared TRATE matrix was not implemented. Candidate local k=2..8 was selected by silhouette; 100 patient bootstrap replicates assessed local stability. Selected local k: UCMC 8, NU 4, Rush 5.',
'3. Fingerprints summarize six-hour modal joint states, support, hypoxemia evidence, hypercapnia/acidemia evidence, and between-bin support changes. State probabilities are released as intervals, not exact sparse counts. Each interval spans at least ten possible patient counts. Cluster sizes and stability are banded. Below ten members, the fingerprint is uninformative. The coordinator reconstructs a feasible representative distribution within each interval.',
'4. Central discovery uses weighted profile distances: joint .4, support .3, hypoxemia .1, hypercapnia .1, support change .1. Sites receive equal weight, with size weighting within site. Candidate global k=2..6 is screened for representation in at least two training sites; silhouette selects among supported candidates. Missingness-dominant, small, or weakly stable local clusters do not define shared templates. The full solution and all three two-site solutions selected four.',
'5. Frozen probability templates are returned to sites. Patients are assigned using the same weighted profile distance, rather than requiring one-to-one local-cluster matches. Clinical assignment requires at least 12 hours of known support, distance at most .75 and nearest-versus-second distance margin at least .01. Others receive I. These are assignment statuses, not cohort exclusions. Six-hour summaries lose some within-bin detail; this is consensus/prototype classification rather than an exact centralized OM/Ward partition.',
'6. Sites fit unadjusted and core-adjusted cause-specific Cox models through day 28 for death, hospice, and their composite. Core covariates: age, sex, race grouped Black/White/other-or-unknown, admission year, noncancer Charlson score, a POA advanced-cancer proxy, and admission code status. These are proxy/available baseline covariates, not full oncology staging or acuity control. Information-insufficient cases have a separate baseline-hazard stratum. Patient IDs are used locally for clustered uncertainty. An explicitly labeled complete-case measured-SOFA sensitivity was also run.',
'7. Site log hazard ratios and robust SEs are combined using REML random effects with modified Hartung-Knapp intervals. Contrasts with small event/non-event partitions are withheld; fewer than two usable site estimates are not pooled. With only three sites, heterogeneity estimates and prediction intervals are imprecise. Separate death effects have only one exportable site and cannot be pooled; hospice has two contributing sites with broad intervals. Measured-SOFA contrasts lack enough exportable/estimable site effects for pooling.',
'8. Prediction uses local ridge logistic models (slope penalty=1; unpenalized intercept) for the observed 28-day hospital composite. Baseline models use the core covariates. Simple summaries add IMV, HFNC, NIV and low-flow duration fractions, observed ARF detection, and support/physiology missingness fractions. Site predictive coefficients are synthesized coefficientwise by inverse variance. Unsupported terms use an explicitly flagged zero prior, not a fabricated association estimate. This is an exploratory two-stage predictor, not a claim of equivalence to a centralized likelihood fit.',
'9. For each held-out site, profile learning, k selection and predictor fitting exclude that site. Frozen profiles and coefficients are then applied to all its hospitalizations. There is no held-out refitting or recalibration. AUC, Brier score and calibration have 200 patient-bootstrap intervals, conditional on the frozen learner. The entire federated discovery process is not refitted in each bootstrap replicate.','',
'PRIVACY AND AUDIT',
'This is a logical federation rehearsal on one machine. Site workers alone read local patient-level caches. Coordinator programs consume allowlisted site exports and coordinator artifacts. Exchange folders contain CSV/JSON aggregates only. No patient identifiers (including hashed identifiers), individual trajectories, individual predictions or raw timestamps were exported.',
'Frozen consensus hashes and source fingerprint hashes were independently checked after outcome fitting. All twelve site/consensus outcome passes used the unchanged frozen definitions. Local checks independently recomputed transported predictions using only two-site frozen coefficients. All original held-out estimates remained numerically unchanged after adding the harder comparator.',
'Interval fingerprints, small-cell rules and an export gate are safeguards, not a differential-privacy guarantee or proof against every reconstruction attack. Actual institutional exchange still requires local disclosure review of the linked output package. Coordinator privacy checks establish code-path/export separation, not operating-system isolation.',
'The datasets and exploratory outcomes had already been examined during earlier analyses. These site-held-out results are an internal-external development rehearsal, not untouched prospective validation. Neither causal treatment effects nor universal prognostic superiority are established.',
'Current inputs preserve the earlier eligibility/outcome rules. They do not recover the three NU or sixty Rush hospitalizations with unclassifiable outcomes outside the current landmark cohort. Oncology stage, ECOG, treatment intent, and extensive missing baseline severity remain unresolved.','',
'ARTIFACTS',
'figures/shared_phenotype_support_profiles.png -- support component of shared prototypes',
'figures/shared_phenotype_arf_profiles.png -- ARF evidence component of shared prototypes',
'figures/pooled_composite_associations.png -- primary core-adjusted composite HRs',
'figures/phenotype_absolute_risks.png -- site-specific observed 28-day risks',
'figures/held_out_auc.png -- four transported predictive models',
'figures/incremental_auc.png -- baseline and harder-comparator increments',
' tables/ contains all site effects, pooled effects, model denominators/status, censored/suppressed curves, composition bands, calibration, Brier and transport estimates.',
'federation_audit.json and frozen_artifact_ledger.json record the privacy/training split and immutable artifact checks.',
'Implementation: diagnostics/federated/README.txt and run_rehearsal.sh. Prior manuscript and pilot outputs are preserved.']
(report/'federated_rehearsal_report.txt').write_text('\n'.join(lines)+'\n')
print(report/'federated_rehearsal_report.txt')
