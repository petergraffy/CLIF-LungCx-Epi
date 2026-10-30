# Implemented CLIF-adapted APACHE-II sensitivity

This extends the frozen October26 three-site72h landmark analysis. It uses
APACHE-II-inspired acute-physiology point thresholds plus age. It does not
implement validated standard APACHE II or APACHE IV and is not a replacement
for the primary phenotype/outcome analysis. The external hospice-project helper
was reviewed but not imported. See `apache_implementation_review.md`.

## Score definition

The window is [first ICU arrival, first ICU arrival+24h), using laboratory result
times and recorded vital/assessment times. Take the greatest point contribution
within that window for each of12 physiology components: temperature, MAP, heart
rate, respiratory rate, oxygenation, acid-base, sodium, potassium, creatinine,
hematology, WBC, and GCS. Add conventional age points. Closed upper-bound gaps
are replaced by complete continuous intervals. Invalid inputs remain unknown.

- Temperature is CLIF temp_c. MAP is recorded valid MAP; only if unavailable,
  calculate (SBP+2*DBP)/3 from unambiguous simultaneous valid measurements.
- Oxygenation uses arterial PaO2 and normalized FiO2. FiO2 is same-specimen lab
  evidence or a prior respiratory value within1h; documented room air permits
  0.21. Future respiratory FiO2 is not carried backward. PaCO2 is required only
  for the FiO2>=0.5 A–a branch and must identify the same gas specimen. Valid
  within-window collection time defines specimen identity when available;
  otherwise use result time. The approximation is FiO2*713-PaCO2/0.8-PaO2,
  explicitly assuming sea level and respiratory quotient0.8.
- Acid-base uses observed arterial pH, with serum bicarbonate substituted when
  no valid arterial pH is available. Substitution is flagged. This remains an
  adapted rule and does not establish that a gas was never drawn.
- Hematology consistently uses estimated hematocrit percent=3*measured
  hemoglobin(g/dL) at all sites. It is an approximation, not measured hematocrit.
  Observed-Hct scores are calculated separately for availability diagnostics;
  sites do not silently mix different hematology definitions in the comparator.
- Creatinine uses the usual concentration thresholds without the unverified
  acute-renal-failure doubling modifier. No modifier is inferred from generic
  billing AKI or cancer diagnoses.
- GCS must be an observed integer3–15 or a complete simultaneous E/V/M tuple.
  Missing GCS is unknown, not15. Sedation-free assessment cannot be verified.
- Numeric values/reference units are normalized explicitly (including Hb g/L
  versus g/dL, creatinine umol/L, hematocrit fractions, oxygen partial pressures
  in kPa and micro-symbol WBC units). Urine-WBC /hpf values are rejected as blood
  counts. Missing unit fields assume canonical CLIF category units and are
  reported; unrecognized explicit units are rejected.

The adapted total is available only when all12 physiology components and age
are available. The maximum for this APS+age definition is62. Unverified chronic-
health/surgical points and acute-renal-failure modifier are recorded as unknown,
not claimed to be absent. The score is therefore not the standard APACHE total.
A separate partial-sum plus12 missing-indicator analysis retains every primary
case; it never changes primary clustering physiology or silently fills its score.

## Frozen measurements and source metadata

Numeric values, result timestamps, cohort members and outcomes come from the
original prepared caches. Source names, units, LOINCs and specimen timestamps
can enrich only exact hospitalization/result-time/category/numeric-value keys.
Conflicting metadata are flagged; newer source observations are not added.
NU tables are now under CLIF-NU/CLIF_2_1_1. Its638,304 frozen lab rows have metadata
matches for611,381 rows;26,923 do not match, and70,996 newer source rows are ignored.
Unmatched metadata retain documented canonical-unit assumptions. UCMC and Rush
match all filtered cached lab rows. These are provenance qualifications, not an
independent validation of ETL or of clinical measurement accuracy.

## Predictive comparison

Use three distinct evaluation cohorts: complete adapted-component cases, all
primary cases with partial-sum/missing indicators, and cases complete for both
the adapted score and the original modified SOFA. Every compared model uses the
same patients within a cohort/endpoint. No score-performance difference between
different samples is described as a head-to-head comparison.

Fit score alone; baseline; baseline+phenotype; baseline+adapted score; baseline+
respiratory summaries; and nested score/summary/phenotype combinations. In the
matched cohort also fit modified SOFA alone and the corresponding baseline/
summary/phenotype variants. Summaries include hour72 support, modality durations,
ARF and support/physiology missingness. The baseline matrix retains the prior
retrospective POA billing covariates; these are not all prospectively available
bedside inputs. Scores contain age points while baseline also adjusts age.

Original held-out phenotype assignments/templates remain frozen. Site ridge
models train locally, and only other-site coefficients enter each held-out
predictor. Synthesis is coefficientwise inverse variance, with flagged zero
priors for unsupported terms. Required score and intercept terms must have
training support before evaluation. This is exploratory internal-external
validation after data inspection, not a newly validated APACHE risk equation.

Endpoints are hospital death and hospital death-or-hospice within28 days after
ICU+72h. Report AUC/Brier, calibration intercept/slope, paired increments,
calibration bins and exploratory decision curves. Two hundred patient bootstraps
condition on frozen learners; calibration intercept/slope are point estimates.
Rare endpoint partitions are explicitly withheld. No held-out recalibration is
used. Complete/matched-score cases are selected by measurement availability.

## Completed run and interpretation

`output/run_APACHE_adapted_20261030_v2/` contains636 complete adapted-component
cases (171 UCMC,391 NU,74 Rush), all3,310 cases in the missing-indicator analysis,
and341 matched-score cases (96/198/47). All score/case/prediction/export audits
pass; original frozen phenotype and outcome artifacts are unchanged.

Among complete adapted-component cases, phenotype additions beyond baseline+
adapted score+respiratory summaries do not improve held-out AUC at any site;
all paired intervals include zero. The all-case missing-indicator analysis shows
small gains at NU/Rush but not UCMC. Matched adapted-score-versus-SOFA AUC
differences are uncertain at every site. These findings do not establish
consistent prognostic superiority, inadequacy of validated scores, or clinical
benefit from deploying the phenotypes.

## Site return contract

The normal site supplement now captures component completeness/provenance/status
under `share/supplement/apache_adapted_*.csv`; candidate score rows stay local.
The full reviewer runner automatically executes the adapted-score model phase
and returns all allowlisted tables under `share/federated/apache_adapted/`.
The coordinator checks common specifications/consensus and rejects mixed phases.
Current versioned site packets are under the completed run's
`site_packets/SITE/share/`, with combined supplement tables in
`reports/site_supplement/apache_adapted/`. Parent packets are preserved.

New executions and scripts are documented in
`diagnostics/apache_adapted/README.txt`. Standard APACHE IV remains unavailable;
this extension addresses only a clearly qualified adapted-score sensitivity.
