Three-site refresh after UCMC ETL correction — 2026-10-26

Active sites: UCMC, NU and Rush. MIMIC is retired and never included here.

Inputs and immutable output locations are set in environment.sh:
- Fresh UCMC source preparation from Box 03-CLIF-2.1/2.1.0, run ID
  ccx_box_20261026_joint. All required source tables are reloaded; no old UCMC
  clinical preparation or respiratory panel is substituted.
- NU/Rush use their validated Oct22 CLIF/waterfall caches. Their local OM fits,
  bootstrap stability, shared assignments, outcome fits and sensitivities are
  rerun. Source panels are independently checked unchanged.
- New primary federation: output/run_federated_rehearsal_20261026
- New reviewer sensitivities: output/run_revision_gaps_20261026
- Cross-run comparison/report: output/run_three_site_refresh_20261026/reports

Original runs and clinical input files are preserved. The default private site
configuration was switched after successful completion to the updated Box source
and fresh UCMC run ID; its previous contents are backed up in the ignored
orchestrator local directory. Do not reuse the old run ID with new source tables.

Execution in this session:
  Rscript diagnostics/multisite/run_site.R \
    output/run_three_site_refresh_20261026/local/ucmc_config.json
  bash diagnostics/refresh20261026/run_after_baseline.sh
run_after_baseline.sh requires that fresh baseline and refuses to overwrite
existing primary/revision outputs via the underlying runners. It executes the
fresh bounded waterfall, full primary federation, all reviewer-gap analyses,
eight refederated inferential sensitivities, descriptive absorbing-state
federation, paired cross-run comparisons and final audits.

Methods are preserved: common-cost OM/Ward local clustering, 100 patient
bootstraps, interval time-resolved fingerprints, outcome-free global k selection,
frozen per-patient assignments, post-ICU+72h cause-specific models, REML modified
Hartung–Knapp pooling, three held-out-site tests, 200 patient-bootstrap intervals,
modified-SOFA and exact respiratory-summary comparators, and retained unknown
states/uncertain assignment statuses. Model/report code now supports variable k;
it never forces the old four-group solution. Separate fits' P-numbers must not
be interpreted as equivalent phenotypes without comparing profiles.

Parameterization:
CLIF_REVISION_PRIMARY_INPUT selects the completed primary federation consumed by
reviewer-gap workers. CLIF_FEDERATED_INPUT_SITE selects each local source cache
for discovery and all sensitivity workers (including absorbing states and source
catalogues). Score-model phenotype terms are taken from the corresponding frozen
consensus, including versions with fewer/more than four groups. Historical
numerical manuscript/response generators are skipped; new reports read actual
aggregate results only.

Validation includes full source fingerprints, exact landmark-member retention,
unchanged NU/Rush input panels, frozen training/site separation and prediction
reproduction, primary/sensitivity hashes, matched complete-score comparisons,
terminal-state boundaries, original primary artifact preservation and export
identifier checks. Cases with unknown support are never excluded.

Remaining scientific limits: missing stage/ECOG/therapy intent; incomplete
modified SOFA (creatinine-only renal); no valid APACHE II/IV comparison;
code-status coverage/semantics; sparse endpoint partitions; only three previously
examined sites. No institutional data transfer or disclosure approval occurs.
Patient data stay in ignored local/ directories. Actual exchange still needs
local review of the linked aggregate package.
