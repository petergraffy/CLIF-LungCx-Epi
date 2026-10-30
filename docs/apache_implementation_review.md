# CLIF hospice-project APACHE review — October 30, 2026

The referenced project contains an APACHE II-named implementation, not APACHE IV.
The reviewed commit is `e71e6a36130e4e011c652fd93ce9433d22a86cb6`.
The [helper](https://github.com/parkegan/clif_hospice_trend/blob/e71e6a36130e4e011c652fd93ce9433d22a86cb6/code/utils/function_apache_II.R)
and its [caller](https://github.com/parkegan/clif_hospice_trend/blob/e71e6a36130e4e011c652fd93ce9433d22a86cb6/code/04_aim2a3_risk_adjusted_variation.R#L796-L816)
were read along with setup and documentation. The caller explicitly uses24h;
the helper's48h default must not be confused with the actual caller. Both anchor
to hospital admission. Our severity window must instead start at first ICU arrival.

## Decision

A corrected, explicitly CLIF-adapted APACHE-II-like comparator is a feasible
supplementary extension. This public code does not establish a fully observed,
validated APACHE II score or an APACHE IV comparison for our current extracts.
The current full-score availability statements therefore remain applicable;
they should not be interpreted as ruling out a separately labeled proxy.
No source scores, frozen outcome models, figures or held-out predictions were
changed during this review. The upstream cohort function was not executed on
clinical records. Its source was not imported into the distributed application.

## Findings and required adaptation

| Finding | Why it matters for this study |
|---|---|
| Hospital-admission anchor | Change to ICU arrival and a bounded first24h window. Retain before-landmark-only information. |
| Hemoglobin replaces hematocrit | This is a modified component. UCMC has no audited hematocrit; NU has limited arterial hematocrit under `other`; Rush has canonical hematocrit. A hemoglobin substitution must be labeled, not described as measured hematocrit. |
| Hemoglobin thresholds90–180 scale with no unit conversion | All three analyzed sites' audited hemoglobin labels use g/dL. The helper assigns4 points to10 and0 points to100. Those are the same concentration under a factor10 unit change, so direct ingestion is incompatible. |
| Missing GCS becomes zero points/GCS15 | Unknown CNS function becomes presumed normal, especially consequential at Rush. Preserve missingness and report observed component coverage. |
| Other missing components summed with `na.rm=TRUE` | A returned numerical score is not evidence of a complete score. Primary measured-component and explicitly imputed sensitivities must be distinct. |
| Chronic-health points set to0 | The optional diagnosis block is commented out. No validated chronic-health or elective/emergency surgical context is recovered by this function. Do not equate missing context with absence. |
| Creatinine ARF doubling not used by caller | The function supports a scalar modifier, but the main calculation never supplies it. A verified acute-renal-failure modifier remains unavailable; any proxy must state this omission. |
| Completeness means at least10 point columns | The intended total has14 components; defaulted GCS/chronic points also count as present. This is not full physiological/clinical completeness. Return true observed components, missing patterns and score variants. |
| Valid decimals fall between closed ranges | For example temperature38.45, MAP109.5 and sodium129.5 produceNA; subsequent summation can turn those gaps into0 contribution. Use complete continuous intervals and boundary tests. |
| Oxygenation join can select future FiO2 | A synthetic ABG1h before the first FiO2 matches that later value with the current `rollends` settings. It is not strictly prior-only despite the comment. Explicitly validate pair times, bounds and maximum distance. |
| Oxygenation requires both PaO2 and PaCO2 even in low-FiO2 branch | Audit branch-specific requirements and missingness rather than unnecessarily turning usable measurements into absent components. Review the A–a equation against a primary specification before reuse. |
| Source numeric values/units and invalid GCS not fully validated | Prefer verified numeric fields and reference units; reject invalid observations rather than clamping them to normal or extreme valid scores. |

The original [APACHE II publication](https://pubmed.ncbi.nlm.nih.gov/3928249/)
describes12 physiology measurements plus age and previous health status. A helper
with substituted/omitted/defaulted components is not that fully observed score.
In lung cancer, immune status and treatment context are particularly material
omissions. Billing metastasis or an unspecified cancer proxy does not supply them.

## Proposed incorporation

1. Implement a separately named modified score in the first24h after ICU arrival,
   with explicit units, complete point intervals, validated pairing and component
   provenance. Prefer observed hematocrit when available; keep any hemoglobin
   proxy in a distinct, consistently defined variant across sites.
2. Keep absent GCS/physiology unknown. Export12-component coverage and missing
   patterns, and distinguish complete observed physiology, proxy substitutions,
   missing-domain assumptions, age points, and unavailable chronic/renal context.
3. Preserve existing phenotypes. Use a new versioned sensitivity for severity-only,
   baseline+severity, +respiratory summaries, and +phenotype comparisons on the
   same evaluation cases. Repeat death/composite held-out AUC, Brier, calibration,
   paired increments and decision curves using only other-site training inputs.
4. Identify the result as a CLIF-adapted APACHE-II-like score, not standard APACHE II
   or APACHE IV. A useful proxy comparison can partially address the severity
   concern, but cannot substantiate superiority over validated APACHE IV.

## Reproducibility

`diagnostics/apache_review/check_upstream.R` loads only reviewed function
definitions and exercises point/join behavior on synthetic inputs. It never
calls the cohort function. The decoded source snapshot, commit/blob identifiers,
probe tables and status are in ignored `output/run_APACHE_review_20261030/`.

```sh
Rscript diagnostics/apache_review/check_upstream.R \
  output/run_APACHE_review_20261030/local/upstream_function_apache_II.R \
  output/run_APACHE_review_20261030/reports
```

The three-site lab units and hematocrit availability cited here come from the
frozen October26 source audits, not a newly executed APACHE validation study.
