Initial respiratory support to cluster/phenotype Sankey

No new patient-level export is needed. The existing descriptor tables contain
this joint distribution: categorical.csv rows with variable=first_hour_support,
cohort=primary_phenotypes (shared) or provisional_local_clusters (phase1), and
group=the cluster label. Exclude Overall rows to avoid double-counting.
The corresponding denominators.csv and availability.csv establish complete
partitions and distinguish structural zero from unavailable information.

Phase1: derive/supplement automatically captures these categories by local
cluster through capture_descriptive_tables. Each site can have different local
clusters: do not sum local labels across sites. The plot script generates only
site-specific figures for these returns.

Phase2: after frozen shared mapping and local per-hospitalization assignment,
recapture descriptors by shared phenotype using diagnostics/descriptives/
site_capture.R (or diagnostics/site_return/capture.R). Returning phase1 local
cross-tabs alone does not reconstruct final shared assignments. Return the new
aggregate descriptors under federated/descriptives/ in the final site packet.
The plotting script verifies a common frozen consensus before summing sites.

Run from repository root:
python3 diagnostics/sankey/plot.py NEW_OUTPUT SITE_SHARE [SITE_SHARE ...]

Accepted inputs are a complete site share (federated/descriptives/), a descriptor
run root (share/descriptives/), or the descriptor folder itself. All declared
source hashes, site labels, availability and cluster denominators are checked.
Only exact observed counts are plotted; unavailable or masked cells are not
imputed. Unknown initial support remains an explicit source node. P1-P4 and I
remain destinations; no support-insufficient stays are dropped.

Initial support is the first ICU-hour (h=0) state in the analysis panel after
configured support processing, not the first later known respiratory label.
The plot is restricted to the primary72h landmark cohort. Early deaths/discharges
are handled in separate descriptive sensitivities, not assigned future primary
phenotypes. This diagram shows initial support to trajectory classification;
it does not display hourly support transitions or outcome flows.

Current pooled and site examples plus exact count tables:
output/run_nearest_assignment_20261030/sankey/
Pooled N3310; unknown initial support215; I130. Source and destination margins
reconcile to the original descriptor denominators. Models are not refitted.

Tests:
python3 tests/sankey.py
Rscript tests/publication_descriptives.R
These verify count conservation, unknown retention, absence of Overall double
counting, suppression/unavailability rejection and phase1 versus shared labels.
