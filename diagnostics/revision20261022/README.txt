Reviewer-gap extension: 2026-10-22

Completed results: output/run_revision_gaps_20261022/reports/revision_gap_report.txt
This extends the frozen three-site primary rehearsal. It preserves all original
primary phenotype/model artifacts. It does not rerun the other seven sites or
replace the submitted ten-site manuscript automatically.

Reproduce from this repository with a FRESH output directory:
  CLIF_REVISION_OUTPUT=output/run_revision_gaps_NEW \
    bash diagnostics/revision20261022/run.sh
The runner rejects an existing output root. It uses the validated site-local
waterfall and baseline caches and original frozen primary rehearsal under
output/run_federated_rehearsal_20261022. These are LOCAL inputs, not exchange
packages. The sensitivity input roots currently follow the three local cache
names; this is a rehearsal runner, not a portable institutional ETL installer.

Completed work:
- Code-status admission/window-end changes, newly documented DNI/comfort,
  outcomes by code, advanced-cancer components, ARF type/timing/late evidence.
- Primary frozen-classifier subgroup and adjustment sensitivities, including
  strict documented-full/no-category-change, COVID and first hospitalization.
- ARF-channel fractions jointly modeled, with interval/rare-partition safeguards.
- Direct modified-SOFA-only and nested severity/summary/phenotype predictors,
  with identical evaluation cases within cohort/endpoint, held-out sites,
  AUC/Brier/calibration/paired increments and exploratory decision curves.
- Complete measured-component SOFA and separate all-case missing-domain-zero
  sensitivity. Neither is claimed to be validated standard SOFA. Renal scoring
  uses creatinine only; adjusted proxy models include missing-domain indicators.
- APACHE II/IV feasibility and explicit unavailable-comparison rows. No actual
  APACHE score can be calculated from the current extracts. Hematocrit source
  labels are audited, including NU's arterial values mapped to other.
- Eight fully refederated inferential sensitivities: 24h,48h,first-detected-ARF,
  support-only, raw-source low-flow FiO2, 0/12h carry and raw-support/no-waterfall.
  Each uses 100 local patient bootstraps and full plus three held-out consensuses.
  Earlier windows use twelve 2/4h bins; known-support assignment threshold scales
  to one sixth of the window. No support-coverage cohort exclusion.
- Separately federated absorbing-state description (absorbing_v2). D groups
  include terminal states and are NEVER used for prognostic regression. The
  initial diagnostic is preserved but superseded; v2 prevents terminal events
  from being ranked as ordinary support changes.
- Reconciled local denominators, aggregate tables, figures, reporting crosswalk,
  analysis-backed reviewer response and manuscript prose drafts.

Role boundaries:
  sites/SITE/local/     records, sequences, models, individual predictions
  sites/SITE/export/    aggregate local returns
  exchange/            staged CSV/JSON only
  coordinator/         frozen severity predictors
  sensitivities/       independent local/exchange/coordinator workflows
  reports/             aggregate results and drafts
No real institutional data exchange occurred. The linked aggregate package
still needs institutional disclosure review; safeguards are not differential
privacy. Do not distribute local/ files or the entire output tree.

Validation:
  Rscript tests/revision_gaps.R
  Rscript diagnostics/revision20261022/site_audit.R UCMC  (also NU/RUSH)
  Rscript diagnostics/revision20261022/central_audit.R
The central audit stages new audits once; run in a fresh root via run.sh.
Local audits check all cases retained, twelve bins, window/event boundaries,
independent frozen score predictions and unchanged original primary artifacts.
Central audits verify 32 consensus content hashes and all 96 outcome versions.
Separate absorbing checks validate exact terminal-event encoding and no outcome
regression. Bootstrap uncertainty conditions on fixed trained models.

Important limitations:
No valid APACHE II/IV scores; incomplete modified SOFA; no stage/ECOG/therapy
intent; substantial code-status/physiology missingness; unclassifiable outcomes;
sparse death/full-code contrasts; only three previously examined sites. Results
favor descriptive phenotyping, not universal prognostic superiority or causal
inference. Manuscript drafts have explicit author/governance placeholders and
have not overwritten the submitted Word files. See the full report.
