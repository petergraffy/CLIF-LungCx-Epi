from pathlib import Path
import csv,json,os
root=Path(os.environ.get('CLIF_UCMC_REFRESH_OUTPUT','output/run_UCMC_etl_refresh_audit_20261026'))
def read(f):
 with (root/'aggregate'/f).open() as h:return list(csv.DictReader(h))
rows=read('final_comparison.csv');co=json.loads((root/'aggregate/cohort_comparison.json').read_text());manifest=json.loads((root/'aggregate/audit_manifest.json').read_text());rec=read('previously_entirely_unknown_recovery.csv')[0];res=read('residual_unknown_source_evidence.csv')[0]
lines=['UCMC CLIF 2.1 respiratory ETL refresh audit — 2026-10-26','',
'The newer Box source shows a major restoration of respiratory device coverage. The configured Desktop source is the old copy: its matched source rows equal the saved prior snapshot after sorting. The Box respiratory file is dated October 23, 2026 and is about 225 MB versus about 54 MB on Desktop.',
'Updated source directory: '+manifest['source_path'],
'Currently configured directory: '+manifest['current_analysis_config_path'],'',
'PAIRED COMPARISON — SAME 650 LANDMARK HOSPITALIZATIONS',
'Method | unknown hours / 46,800 | unknown % | entirely unknown stays | ARF-positive stays']
for x in rows:lines.append(f"{x['method']} | {x['unknown_hours']} | {float(x['unknown_pct']):.3f}% | {x['entirely_unknown_stays']} | {x['arf_positive_stays']}")
lines+=['',
'All 64 previously entirely unknown stays now have identifiable source device categories during the window. None remains entirely unknown even before waterfall. The raw-source-plus-six-hour-carry comparison improves from 31.04% to 6.91% unknown hours; this is not only an effect of applying waterfall. With the same bounded CLIFpy 2.1.0 waterfall, unknown hours decrease from 6,166 to 30 (13.18% to 0.064%). The 30 residual hours span 11 hospitalizations: 29 have no same-hour respiratory row; one has a source FiO2 value; none has positive source oxygen flow.',
'The post-2019 discontinuity is substantially reduced. In 2024, raw-source/six-hour-carry unknown hours fall from 57.85% to 8.86%; in 2023 they fall from 49.48% to 10.32%. Remaining intermittent gaps are mostly resolved by the existing bounded waterfall approach.', '',
'COHORT / TIMING / QUALITY CHECKS',
'Both old and updated anchor tables produce the same 849 eligible hospitalizations and the same 650 +72h landmark hospitalizations. In the matched cohort, ICU origin, hospital admission/discharge, age, event type/time/validity, sex, race and ethnicity have zero changes. Patient counts and outcome timings therefore do not explain this respiratory coverage improvement.',
'The respiratory timestamp schema changed from timezone-annotated UTC to unannotated timestamps. After the existing UTC normalization, old/new event-state keys align most strongly at zero shift (19,759 matches versus roughly 4,800–5,000 with +/-5/6-hour shifts). No heuristic shift was applied. This supports unchanged relative timing but does not independently establish the producing system’s timezone semantics.',
'The new table contains additional rows and substantially more populated device labels: within the fixed 72h windows, known-category source rows increase from 14,912 to 25,469 and room-air rows from 1,425 to 3,676. Every populated device category maps successfully. Source-table coverage is improved; whether any newly populated values were inferred during ETL requires producing-code provenance and must not be described automatically as new direct bedside observations.',
'Applying the updated respiratory source changes 13,700 hourly support states and 4,892 physiological states in the paired raw comparison. Detected ARF rises from 215 to 237 hospitalizations using the same prior vitals/labs and unchanged detection rules. The rest of the clinical data is held fixed here to isolate the respiratory-table change; a full fresh-source run may additionally reflect other updated tables.',
'Waterfall physiology is NOT used for the conservative updated comparison: the refreshed raw-evidence ARF states are preserved when adding waterfall support. Numeric bfill remains disabled, and all source/output rows are bounded before ICU+72h. Review flags remain for 13 source-positive-flow rows assigned room air by waterfall, and 5,866 nasal-cannula rows with filled FiO2=0.21. These are preprocessing review flags, not proof of erroneous documentation, and support retaining the raw-evidence physiology policy.', '',
'ANALYSIS STATUS',
'All old UCMC, NU, Rush and federated primary artifacts remain unchanged; frozen hashes were checked. This is a paired ETL audit, not a new clustering or pooled-outcome result. No patients are excluded for missing support.',
'The next analysis must use a fresh UCMC run ID and the Box directory above, rebuild UCMC clinical derivation and sensitivities, regenerate its fingerprints, then relearn/version shared templates and rerun mapping/outcomes/pooling/transportability. Old assignments and pooled estimates should not simply be relabeled as refreshed results. Cached preparation is not reused across source versions.', '',
'FILES',
'aggregate/final_comparison.csv — matched raw and waterfall comparisons.',
'aggregate/previously_entirely_unknown_recovery.csv — recovery of the previous 64 stays.',
'aggregate/cohort_comparison.json and anchor_field_changes.csv — cohort/timestamp/outcome consistency.',
'aggregate/unknown_by_year.csv and figures/unknown_support_by_year.png — annual coverage improvement.',
'aggregate/source_summary.csv, source_state_counts.csv, timestamp_alignment.csv — direct source checks.',
'aggregate/residual_unknown_source_evidence.csv and waterfall_conflict_flags.csv — remaining review points.',
'aggregate/final_audit.json — retention, boundary and frozen-artifact checks.',
'All identifiers, source rows and hourly panels remain under local/. These aggregate diagnostics are local review outputs.']
(root/'etl_refresh_report.txt').write_text('\n'.join(lines)+'\n')
print(root/'etl_refresh_report.txt')
