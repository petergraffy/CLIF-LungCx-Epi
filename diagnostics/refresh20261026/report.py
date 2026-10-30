"""Report refreshed results using coordinator aggregates only."""
from pathlib import Path
import csv,json,os
root=Path('output/run_three_site_refresh_20261026');out=root/'reports';fed=Path(os.environ['CLIF_REVISION_PRIMARY_INPUT']);rev=Path(os.environ['CLIF_REVISION_OUTPUT'])
def read(p):
 with Path(p).open() as f:return list(csv.DictReader(f))
def f(x):return float(x) if x else float('nan')
def ci(z):return f"{f(z['estimate']):.3f} ({f(z['conf_low']):.3f} to {f(z['conf_high']):.3f})"
sources=read(out/'tables/source_comparison.csv');metrics=read(fed/'reports/tables/transport_metrics.csv');pool=read(fed/'reports/tables/pooled_associations.csv');labels=read(fed/'reports/tables/phenotype_labels.csv');ph={r['phenotype']:r['label'] for r in labels};m=json.loads((fed/'coordinator/full/frozen_manifest.json').read_text());k={s:json.loads((fed/'coordinator'/('loo_'+s)/'frozen_manifest.json').read_text())['k'] for s in ['UCMC','NU','RUSH']}
lines=['Three-site full refresh after UCMC respiratory ETL update','2026-10-26','',
'SCOPE AND INPUTS',
'UCMC was rebuilt from all required updated Box CLIF tables in 03-CLIF-2.1/2.1.0 using a fresh source preparation and run ID. NU/Rush use their validated October 22 source caches; all three sites repeat local clustering, stability, central consensus, frozen local outcomes, pooling and held-out validation. MIMIC is excluded. Prior outputs are preserved.',
f"Fresh primary consensus: k={m['k']}, ID={m['consensus_id']}. Held-out discovery k: {k}.",
'No number of shared phenotypes or clinical label was forced. Labels below describe predominant support in the learned templates, not biological endotypes or an assumed ordering by mortality.','',
'PAIRED SOURCE / CLINICAL COHORT COMPARISON',
'Site | version | N | unknown support hours % | raw ARF positive | measured SOFA complete | local k']
for z in sources:lines.append(f"{z['site']} | {z['version']} | {z['n']} | {f(z['unknown_pct']):.3f}% | {z['arf_positive_n']} | {z['complete_sofa_n']} | {z['local_k']}")
lines+=['Raw-evidence means the refreshed CLIF table before our bounded analysis waterfall; it is not independently verified raw EHR observation provenance. Source-level ETL inference or lookahead cannot be assessed without the producing transformation code. Analysis-side input/output windows are enforced as specified.',
'The same primary landmark members are retained. No support-coverage exclusion is used. Information-insufficient/ambiguous cases remain I and contribute to the analyses under the stated models. NU/Rush source panels are verified unchanged, while their shared assignments can change because the consensus is relearned.', '',
'PRIMARY PHENOTYPES AND POST-LANDMARK ASSOCIATIONS']
for z in labels:lines.append('  '+z['label'])
for z in pool:
 if z['model']=='core' and z['endpoint']=='death_or_hospice':
  if z['status']=='ok':lines.append(f"{ph.get(z['phenotype'],z['phenotype'])} versus P1: HR {f(z['hazard_ratio']):.3f} (95% CI {f(z['hr_low']):.3f}–{f(z['hr_high']):.3f}); sites {z['k']}; I² {f(z['i2_percent']):.1f}%; prediction interval {f(z['hr_prediction_low']):.3f}–{f(z['hr_prediction_high']):.3f}.")
  else:lines.append(f"{z['phenotype']}: {z['status']} ({z['k']} usable sites).")
lines+=['Local bootstrap mean-Jaccard bands (released decile bands, not confidence intervals):']
for site in ['UCMC','NU','RUSH']:
 clusters=read(fed/'exchange/discovery'/site/'clusters.csv')
 lines.append(site+': '+', '.join('local '+z['local_cluster']+' = '+f"{f(z['jaccard_lower']):.1f}–{f(z['jaccard_upper']):.1f}" for z in clusters))
lines+=['Some local groups can remain weakly stable even when the shared solution selects the same k across held-out fits. Shared cluster-count reproducibility is not proof that every local group is stable.', 'Endpoint: hospital death or hospice within 28 days AFTER the ICU+72h landmark. Core adjustment: age, sex, race grouping, admission year, noncancer Charlson, POA advanced-cancer proxy and admission code status; patient-clustered uncertainty, separate baseline-hazard stratum for I. Site log hazards are pooled with REML modified Hartung–Knapp intervals. Associations are conditional, not causal.',
'Separate death and hospice fits and actual contributor counts/statuses are in the pooled tables. Rare contrasts are withheld; missing site effects are never assigned zero. Old versus new labels refer to independently learned partitions, so same P-number does not establish clinical equivalence.', '',
'HELD-OUT PREDICTION',
'Site | baseline | +phenotype | +simple summaries | +summaries+phenotype (AUC, 95% patient-bootstrap CI)']
for site in ['UCMC','NU','RUSH']:
 vals=[]
 for model in ['baseline','plus_phenotype','simple_summaries','summaries_plus_phenotype']:
  z=next(z for z in metrics if z['site']==site and z['model']==model and z['metric']=='auc');vals.append(ci(z))
 lines.append(site+' | '+' | '.join(vals))
for z in metrics:
 if z['model']=='phenotype_beyond_summaries' and z['metric']=='auc':lines.append(f"Phenotype beyond summaries, {z['site']}: ΔAUC {f(z['estimate']):+.4f} ({f(z['conf_low']):+.4f} to {f(z['conf_high']):+.4f}).")
lines+=['Predictors/templates exclude the held-out site; no held-out refit/recalibration. Coefficient synthesis is an exploratory two-stage predictor. Bootstrap uncertainty conditions on frozen training. Earlier data have been examined, so this is internal-external development, not untouched prospective validation. Harder severity comparisons below explicitly include hour-72 support in addition to modality durations.', '',
'MODIFIED SOFA AND HARDER NESTED COMPARISON']
score=read(rev/'reports/tables/score_transport_metrics.csv');inc=read(rev/'reports/tables/score_paired_increments.csv');meta=read(rev/'reports/tables/score_model_metadata.csv')
for cohort in ['measured_sofa','missing_domain_zero_sensitivity']:
 lines.append(cohort+': modified SOFA alone | baseline+SOFA | baseline+phenotype | baseline+SOFA+summaries | +phenotype')
 for site in ['UCMC','NU','RUSH']:
  vals=[]
  for model in ['sofa_only','sofa','phenotype','sofa_summaries','sofa_summaries_phenotype']:
   z=next(z for z in score if z['site']==site and z['cohort']==cohort and z['endpoint']=='death_or_hospice' and z['model']==model and z['metric']=='auc');vals.append(ci(z) if z['status']=='ok' else z['status'])
  lines.append(site+' | '+' | '.join(vals))
 for z in inc:
  if z['cohort']==cohort and z['endpoint']=='death_or_hospice' and z['comparison']=='phenotype_beyond_sofa' and z['metric']=='auc':lines.append(f"{z['site']} phenotype increment beyond SOFA+summaries: {f(z['estimate']):+.4f} ({f(z['conf_low']):+.4f} to {f(z['conf_high']):+.4f}).")
lines+=['Score comparisons use the same evaluation patients across models within each score cohort/endpoint. Modified SOFA uses measured first-24h components and creatinine-only renal scoring; it is not validated standard SOFA. The all-case zero-domain/missing-indicator sensitivity is labeled explicitly and does not treat missing physiology as normal in primary clustering. APACHE II/IV remain unavailable from current inputs; no incomplete proxy is called APACHE. Decision curves, calibration, Brier, actual denominators and separate death prediction are exported.', '',
'REFEDERATED SENSITIVITIES',
'Analysis | global k | k when UCMC / NU / Rush held out']
cons=read(rev/'reports/tables/sensitivity_consensus.csv')
for v in ['h24','h48','onset','support_only','low_flow','carry0','carry12','raw_support']:
 vals=[next(z['k'] for z in cons if z['analysis']==v and z['version']==ver) for ver in ['full','loo_UCMC','loo_NU','loo_RUSH']];lines.append(v+' | '+vals[0]+' | '+' / '.join(vals[1:]))
lines+=['Each sensitivity repeats local common-cost OM/Ward clustering, 100 patient bootstraps, approved-form aggregate fingerprints, full plus three held-out consensus versions, twelve local outcome passes, meta-analysis and held-out prediction. Earlier windows use twelve relative-time bins (2/4 hours), and the known-support assignment threshold scales with window length. Counts and selected labels are not assumed to match across sensitivity populations.',
'The separate absorbing-state federation includes early death, hospice and live discharge as terminal states; descriptive D groups are not primary P groups, and no prognostic regression uses those outcome-containing exposures. Code-status adjustments, strict full-code, first-hospitalization, COVID exclusion and ARF/high-support subgroup models use the frozen refreshed primary classifier.', '',
'ARF / TREATMENT-LIMITATION CHECKS',
'Updated ARF type/timing/persistence, independent channel terms, code-status transitions and source catalogues are in the reviewer-gap tables. Unknown physiology is not negative evidence; first detection is documentation-based, not biological onset. Mixed positive hours overlap the two individual channels. These models do not establish independent biological endotypes.',
'Code-status coverage and source semantics remain imperfect. A lack of a recorded comfort/DNI category is not absence of a clinical preference. Cancer stage, ECOG and treatment intent remain unavailable. Full-code and separate-death contrasts may lack enough releasable contributors; inspect statuses rather than interpreting absent effects as null.', '',
'PROVENANCE / REVIEW',
'Original three-site and retired diagnostic outputs are preserved. Primary and reviewer-gap frozen content hashes, site retention and transported predictions are audited. No patient identifiers, sequences or individual predictions are sent to the coordinator. This remains a logical federation rehearsal on one computer; aggregate disclosure review is needed before institutional exchange.',
'Author IRB/waiver details, funding/conflicts and completion of final manuscript/checklist pages still require author confirmation. Historical response/manuscript drafts were not regenerated with old numerical constants. This report is the current result summary.', '',
'OUTPUTS',str(fed/'reports'),str(rev/'reports'),str(out/'tables'),'Reports/figures include shared support–ARF profiles, pooled effects, absolute risks, held-out prediction, modified-SOFA comparisons, incremental performance, decision curves, ARF contribution and earlier landmarks.']
(out/'three_site_refresh_report.txt').write_text('\n'.join(lines)+'\n')
print(out/'three_site_refresh_report.txt')
