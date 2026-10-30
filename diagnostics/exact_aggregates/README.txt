Exact aggregate counts (current project policy)

Exact counts, including rare cells and numbers at risk, are now the default.
No disclosure masking or rounding is applied to these newly generated returns.
Patient IDs, hospitalization IDs, raw event timestamps, individual sequences and
predictions remain local. Non-estimable or unavailable quantities remain missing.
Historical banded exports are retained solely for provenance; do not infer exact
counts from their endpoints or midpoints. A legacy R option exists to reproduce
those historical exports/tests; it is not enabled in the current workflow.

Current results: output/run_exact_aggregates_20261030/
Only site_packets/SITE/share/ is intended for exchange. It contains:
- federated/descriptives: Table1/Table2, denominator/missingness/moment tables.
- federated/survival_tdauc: KM, risk sets, AUC(t), fixed-grid tdROC and methods.
- federated/competing_events: daily Aalen-Johansen estimates and exact risk sets.
- exact_counts: screening/annual/landmark flows, all primary/sensitivity cohort
  counts, saved Cox and score model denominators, original severity completeness,
  exact hourly support/ARF counts, cluster sizes/stability and exact supplement.

These are extensions of frozen results: no phenotype learning, classifier or
clinical model is refit. Earlier fitted coefficient/prediction results are kept
separately; removing count suppression does not make invalid fits estimable.
Interval-valued phenotype templates remain the frozen scientific representation;
new discovery exports additionally retain exact underlying aggregate state counts.

Reproduce each site, opening only that site's local caches:
Rscript diagnostics/descriptives/site_capture.R SITE NEW_DESCRIPTOR_ROOT
Rscript diagnostics/manuscript_figures/site_curves.R SITE PRIMARY_ROOT NEW_CURVE_ROOT
Rscript diagnostics/manuscript_figures/site_survival_tdauc.R SITE PRIMARY_ROOT NEW_SURVIVAL_ROOT
Rscript diagnostics/exact_aggregates/site_counts.R SITE NEW_COUNT_ROOT

The three-site rehearsal defaults point to the verified October26/30 caches.
The package.py script stages these fresh phases into the versioned exact packet.
Coordinator collection (aggregate files only):
python3 diagnostics/exact_aggregates/collect.py NEW_OUTPUT SITE_A_SHARE SITE_B_SHARE SITE_C_SHARE

The collector verifies site/version/file hashes, excludes identifier columns and
rejects count bands or disclosure-masked cells. Exact continuous moments pool
means/SD; categorical partitions retain total/observed/missing denominators.
Do not average site medians or sum site-distinct patient counts into a global
unique-patient count. ICU/hospital clocks and primary/sensitivity cohorts stay
explicitly labeled. Whole-stay descriptors do not become prognostic covariates.

Plot the exact versions with:
CLIF_SURVIVAL_OUTPUT=output/run_exact_aggregates_20261030 python3 diagnostics/manuscript_figures/plot_survival_tdauc.py
CLIF_FIGURE_OUTPUT=output/run_exact_aggregates_20261030 CLIF_DESCRIPTOR_REPORTS=output/run_exact_aggregates_20261030/reports/publication_tables python3 diagnostics/manuscript_figures/refinements.py

Tests: Rscript tests/exact_aggregates.R
Historical masking tests explicitly enable legacy mode in their own R sessions.
