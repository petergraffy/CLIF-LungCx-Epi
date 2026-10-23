# CLIF lung cancer respiratory trajectories

This is a new federated CLIF run for the CCX revision. It rederives respiratory trajectories from source CLIF tables. Existing returns under `sites/` and prior manuscript outputs are historical and are not inputs to this workflow.

The submitted implementation, including local edits present on October 21, 2026, is preserved under [archive/2026-10-21](archive/2026-10-21/README_ARCHIVE.md). New code does not source it.

## Site workflow

Run from the repository root. Create `config/config.json` using [the template](config/config_template.json). Keep scientific analysis settings consistent across sites and choose a new `run_id` whenever settings or the source version change.

```sh
Rscript code/00_setup.R
Rscript code/02_derive.R
```

Preparation is included in derivation; `Rscript code/01_prepare.R` can run it separately. Review aggregate audits and return **only**:

```text
output/run_<SITE>_<RUN_ID>/share/
```

**Never return `local/` or the whole run directory.** Local caches contain patient identifiers, dates, clinical data, and sequences. They are required for the second pass.

After the coordinator returns the frozen consensus mapping:

```sh
Rscript code/03_model.R path/to/consensus_mapping.csv
```

Return the resulting `share/models_<CONSENSUS_ID>/` after local disclosure review. `run.sh` and `run.bat` support `setup`, `prepare`, `derive`, `model <mapping>`, and `test`.

## Coordinator workflow

Supply explicit directories for new returns. Historical `sites/` directories are never discovered automatically.

```sh
Rscript code/postprocessing/01_consensus.R output/consensus_v1 returns/siteA/share returns/siteB/share
# Return output/consensus_v1/consensus_mapping.csv to sites.
Rscript code/postprocessing/02_pool.R output/pooled_revision_v1 returns/siteA/models_ID returns/siteB/models_ID
```

Consensus matching uses hourly respiratory state profiles and preserves time ordering. Matching, labels, and references use no mortality, hospice, SOFA, vasopressor, or length-of-stay information. Terminal-state sensitivity phenotypes are descriptive and excluded from outcome modeling. Pooled effects use REML and modified Hartung–Knapp inference, with site estimates, heterogeneity, and prediction intervals.

## Analysis specification

The primary analysis derives joint respiratory-support × ARF phenotypes among adults with primary lung cancer present on admission who are alive and hospitalized at first ICU admission plus 72 hours. Follow-up starts at that ICU landmark. Early deaths, hospice discharges, and live discharges are counted separately and represented in a full-cohort absorbing-state sensitivity.

The workflow includes 24- and 48-hour landmarks, support-only clustering, first-detected-ARF anchoring with follow-up after its own 72-hour landmark, estimated low-flow FiO₂ sensitivity, ARF persistence sensitivities, patient bootstrap stability, and leave-one-site-out harmonization. Missing support and physiology are explicit states. No hospitalization is excluded for insufficient support documentation, including entirely unknown support sequences. HFNC, CPAP, and NIPPV remain separate.

Read the [analysis protocol](docs/analysis_protocol.md), [reviewer implementation map](docs/reviewer_implementation.md), and [export contract](docs/export_contract.md) before distributing a run. These document estimands, assumptions, validation boundaries, and analyses requiring additional data.

## Verification

```sh
Rscript tests/run_tests.R
Rscript tests/test_descriptives.R
Rscript tests/run_integration.R
```

Tests use synthetic data, including ward-to-ICU transfers, early events, missing observations, variable cluster counts, suppressed cells, and a complete two-site round trip. The integration test also runs the regression checks and takes several minutes. Real site execution and clinical validation remain necessary before interpreting new results.

## Respiratory support missingness audit

Run `Rscript diagnostics/01_audit_support.R` after a completed derivation. This diagnostic traces unknown hours to source records, care location, and calendar year, then tests an explicitly exploratory oxygen-unspecified state using measured flow and bounded forward carry. It preserves the primary derivation and stores all outputs locally under an ignored `output/run_*/` directory. `CLIF_SUPPORT_AUDIT_OUTPUT` can select a different local audit directory. The diagnostic is separate from the prespecified analysis registry.

## CLIF waterfall comparison

`diagnostics/02_compare_waterfall.sh` compares the configured completed baseline with official CLIFpy waterfall preprocessing, preserving patient-specific landmark and terminal-event boundaries. It reruns all nine derivations and adds a comparison that holds ARF states fixed while changing support. Outputs remain in a separate ignored run directory; baseline files and the production primary specification remain unchanged. See `diagnostics/waterfall/README.txt` for the pinned environment, exact boundaries, provenance, and source-context QA.

Reviewer-gap and direct severity comparisons are implemented in
[`diagnostics/revision20261022/README.txt`](diagnostics/revision20261022/README.txt).
This extends the three-site development rehearsal with refederated sensitivities,
ARF/code-status characterization, direct **modified-SOFA** comparisons and an
explicit APACHE feasibility audit. APACHE II/IV cannot be validly calculated from
the available extracts; absent performance is never replaced by a proxy result.
