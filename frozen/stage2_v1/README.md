# Frozen Stage 2 bundle — v1

This folder ships with the repository. It contains aggregate phenotype templates,
the fixed analysis plan, checksums and original Stage 1 acceptance records.
See [site instructions](../../code/README.md) before running.

- Primary: four joint respiratory-support/ARF phenotypes, learned from all ten
  sites; 72-hour landmark, followed by 28 days of hospital outcome assessment.
- Discovery sensitivities: 24h, 48h, ARF onset, support-only, low-flow FiO2,
  ARF carry 0h/12h, raw support, and descriptive absorbing states.
- Full and leave-one-site-out mappings are frozen for every variant. All ten
  sites participate in every Stage 2 analysis, including Penn and Rush.
- Penn's empty uploaded h24/carry12 discovery files contribute no fabricated
  fingerprints. Those two mappings use nine available discovery contributions;
  Penn applies the same frozen mappings and returns all local models.
- Existing legacy discovery snapshots are accepted by explicit coordinator
  instruction. Acceptance is not a claim of source-code validation where source
  changes were unavailable. Original cached cohorts/instants are preserved.
- Each site runs one command and returns one aggregate packet. Prediction uses
  local patient-grouped cross-validation; LOO mappings assess phenotype
  transportability, not external validation of an outcome predictor.
- SOFA, adapted APACHE and summary comparators, four cause-specific Cox endpoints,
  clinical sensitivities, model diagnostics, Table 1/2 descriptors, exact hourly
  trajectories, Sankey inputs, ROC/tdAUC and pooled KM/AJ inputs are included.

The manifest binds the exact analysis source, R lock and file inventory. Do not
edit this folder or the Stage 1 manifests. The R environment is restored using
`scripts/setup.R`; R 4.4.2 and Python 3.12 are required. No patient data, local
configurations, clinical caches or institutional upload folders are distributed.
A late site can be admitted using `scripts/admit_stage2_site.R` without changing
the full frozen prototypes; see the coordinator instructions.
