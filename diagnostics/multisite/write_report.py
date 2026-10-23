"""Write a local, reproducible plain-text summary of the three-site comparison."""
from pathlib import Path
import csv
root=Path('output/run_multisite_comparison_20261022')
def read(name):
 with (root/'aggregate'/name).open() as f:return list(csv.DictReader(f))
sites=read('site_comparison.csv');clusters=read('support_waterfall_raw_arf_clusters.csv')
lines=['UCMC / NU / RUSH comparison -- 2026-10-22','',
'Same adult POA lung-cancer ICU eligibility rules; admissions in [2018-01-01, 2025-01-01).',
'The primary exposure window starts at first ICU arrival and ends at +72 hours.',
'Landmark patients remain alive and hospitalized at +72 hours. Outcomes occur after the landmark.',
'No patient is excluded for missing respiratory support. Unknown remains an explicit state.','',
'The conservative comparison uses waterfall respiratory support with the original hourly ARF states.',
'Full waterfall (including filled numeric respiratory fields) is a separate sensitivity.',
'Waterfall inputs/outputs are bounded before landmark or terminal event; numeric backfill is disabled.',
'Baseline severity and clinical covariates are held fixed. First detected ARF anchors use raw evidence.','',
'SUMMARY','Site | Eligible ICU stays | Landmark N | Unknown hours: original -> waterfall | Wholly unknown: before -> after | Local clusters: original / support+raw ARF / full waterfall | Post-landmark death / hospice']
for s in sites:
 lines.append(f"{s['site']} | {s['eligible_hospitalizations']} | {s['primary_n']} | {float(s['original_unknown_pct']):.2f}% -> {float(s['waterfall_unknown_pct']):.2f}% | {s['entirely_unknown_before']} -> {s['entirely_unknown_after']} | {s['original_k']} / {s['support_waterfall_raw_arf_k']} / {s['full_waterfall_k']} | {s['deaths']} ({float(s['death_pct']):.1f}%) / {s['hospice']} ({float(s['hospice_pct']):.1f}%)")
lines+=['','SOURCE ROWS IN WHOLLY UNKNOWN STAYS']
for u in read('entirely_unknown_source_audit.csv'):
 lines.append(f"{u['site']}: {u['entirely_unknown_n']} wholly unknown stays; {u['no_respiratory_rows_in_72h_n']} have no respiratory-table rows in the 72 hours; {u['any_positive_source_flow_n']} have positive source oxygen flow; {u['any_source_fio2_n']} have source FiO2. Flow alone is not assigned a specific support device.")
lines+=['','MEASUREMENT AND STABILITY']
for s in sites:
 site=s['site'];n=int(s['primary_n']);j=[float(c['mean_jaccard']) for c in clusters if c['site']==site]
 lines.append(f"{site}: raw ARF evidence in {s['raw_arf_n']}/{n} ({100*int(s['raw_arf_n'])/n:.1f}%); full waterfall ARF in {s['full_waterfall_arf_n']}/{n} ({100*int(s['full_waterfall_arf_n'])/n:.1f}%). Complete measured-component SOFA available in {s['measured_sofa_n']}/{n}; admission code status unknown in {s['admission_code_unknown_n']}/{n}. Mean patient-bootstrap Jaccard across local clusters ranges {min(j):.2f} to {max(j):.2f} (100 replicates per cluster).")
lines+=['','LOCAL SUPPORT + ORIGINAL ARF CLUSTERS',
'Cluster numbers do not correspond across sites. Dominant support is a descriptive mean over 72 hours, not a matched phenotype.']
for c in clusters:
 lines.append(f"{c['site']} C{c['local_cluster']}: n={c['n']}; dominant={c['dominant_support']} ({float(c['mean_hours_pct']):.1f}% of hours); unknown={float(c['unknown_hours_pct']):.1f}%; ARF-positive stays={c['arf_n']}; death={c['death_n']} ({float(c['death_pct']):.1f}%); hospice={c['hospice_n']}; mean Jaccard={float(c['mean_jaccard']):.2f}.")
lines+=['','IMV-DOMINANT LOCAL GROUPS']
for c in clusters:
 if c['dominant_support']=='IMV':
  combined=int(c['death_n'])+int(c['hospice_n'])
  lines.append(f"{c['site']}: n={c['n']}, death/hospice={combined}/{c['n']} ({100*combined/int(c['n']):.1f}%). This is a descriptive comparison of local IMV-dominant groups, not a formal cross-site phenotype match.")
lines+=['','INTERPRETATION AND LIMITS',
'NU has much more complete respiratory-support coverage than UCMC or Rush after the same waterfall procedure. The remaining missingness is strongly site dependent. Reduction through imputation does not prove that a source extract was defective.',
'Full numeric waterfall increases detected ARF at every site. These additional labels depend on filled values; they are not additional observed measurements. Source conflict checks and review queues are retained locally.',
'Unadjusted hospital outcomes differ. They are descriptive and should not be interpreted as site effects. Hospital death and hospice are separate outcomes. Rates here cover the remainder of hospitalization, not a fixed 28-day risk.',
'Cluster numbers and dominant-support labels are local. Outcome-free phenotype matching is required before pooled phenotype associations. A stable unknown-support cluster can reflect stable missingness rather than a clinical phenotype.',
'Measured-component severity and admission code status remain incomplete. Cancer stage, functional status and treatment intent are not supplied as an oncology extension in these datasets.',
'No new cross-site or fold-contained predictive validation has been performed. The full waterfall association models are exploratory; they do not validate the conservative arm as a prognostic model.',
'NU has 3 and Rush has 60 hospitalizations with unclassifiable outcomes under current rules. Rush source names are Needs Coding by HIM (59) and a transfer categorized Other (1). These exclusions are unrelated to support missingness.',
'Calendar comparisons use supplied timestamps. Confirm whether Rush deidentification shifts dates before interpreting secular trends.','',
'INPUT AUDIT',
'After completing these analyses, the core date parser and projected-ID filter were fixed and regression-tested. Invalid/reversed study boundaries are now rejected. Completed calculations retain their original source fingerprints; code/ and diagnostics/ are preserved in local/analytical_source_before_input_fixes.tar.gz. New runs with the fixed core require fresh run IDs.',
'The frozen date parser did not parse date-only boundaries. UCMC and NU included cohorts were verified entirely within the intended dates; Rush was restarted with explicit timestamp boundaries. The aborted Rush v3b run (2017/2025 cases included) is excluded. Only Rush v3c enters these results.',
'Rush code_status contains an extra hospitalization_id. A local projection omitted that unused field so the frozen loader could apply its patient-ID filter. Every retained code-status row/value was verified identical to source. Original data files were not changed.',
'All nine baseline analyses and ten waterfall comparisons are complete at each new site. Every included cohort admission is checked against exact study boundaries. One assignment per eligible hospitalization, sequence lengths, unchanged event totals, and landmark/terminal boundaries are checked.',
'Aggregates in this folder are local exploratory results and have not passed disclosure review for external sharing. Patient-level records and review queues remain in ignored local output directories.','',
'FILES',
'aggregate/site_comparison.csv -- principal site metrics',
'aggregate/support_waterfall_raw_arf_clusters.csv -- local profiles, outcomes, stability',
'aggregate/all_analyses_comparison.csv -- all sensitivity comparisons',
'aggregate/unknown_by_year.csv -- annual denominators and unknown hours',
'aggregate/source_conflict_checks.csv -- source versus waterfall review flags',
'figures/unknown_support_by_site_year.png -- annual coverage comparison',
'figures/*_support_waterfall_raw_arf_clusters.png -- support profiles of joint-state clusters',
'local/comparison_manifest.json -- input runs and shared waterfall source hash']
(root/'comparison_report.txt').write_text('\n'.join(lines)+'\n')
print(root/'comparison_report.txt')
