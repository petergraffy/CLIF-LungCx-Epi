#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
export CLIF_FEDERATED_OUTPUT=output/run_nearest_assignment_20261030/primary
export CLIF_REVISION_PRIMARY_INPUT="$CLIF_FEDERATED_OUTPUT"
export CLIF_APACHE_PRIMARY="$CLIF_FEDERATED_OUTPUT"
export CLIF_REVISION_OUTPUT=output/run_nearest_assignment_20261030/revision
export CLIF_APACHE_OUTPUT=output/run_nearest_assignment_20261030/apache_adapted
export CLIF_EXACT_OUTPUT=output/run_nearest_assignment_20261030/exact
mkdir -p "$CLIF_EXACT_OUTPUT/logs"
NEAREST_PIDS=()
for NEAREST_SITE in UCMC NU RUSH; do
 (
  Rscript diagnostics/descriptives/site_capture.R "$NEAREST_SITE" "$CLIF_EXACT_OUTPUT/sites/$NEAREST_SITE/descriptives"
  Rscript diagnostics/manuscript_figures/site_curves.R "$NEAREST_SITE" "$CLIF_FEDERATED_OUTPUT" "$CLIF_EXACT_OUTPUT/sites/$NEAREST_SITE/competing_events"
  Rscript diagnostics/manuscript_figures/site_survival_tdauc.R "$NEAREST_SITE" "$CLIF_FEDERATED_OUTPUT" "$CLIF_EXACT_OUTPUT/sites/$NEAREST_SITE/survival"
  Rscript diagnostics/exact_aggregates/site_counts.R "$NEAREST_SITE" "$CLIF_EXACT_OUTPUT/sites/$NEAREST_SITE/counts"
  Rscript diagnostics/site_return/model_diagnostics.R "$NEAREST_SITE" "$CLIF_FEDERATED_OUTPUT" "$CLIF_EXACT_OUTPUT/sites/$NEAREST_SITE/model_diagnostics"
 ) > "$CLIF_EXACT_OUTPUT/logs/$NEAREST_SITE.log" 2>&1 &
 NEAREST_PIDS+=("$!")
done
for NEAREST_PID in "${NEAREST_PIDS[@]}"; do wait "$NEAREST_PID"; done
python3 diagnostics/exact_aggregates/package.py
Rscript diagnostics/descriptives/collect.R "$CLIF_EXACT_OUTPUT/reports/publication_tables" "$CLIF_EXACT_OUTPUT/site_packets/UCMC/share" "$CLIF_EXACT_OUTPUT/site_packets/NU/share" "$CLIF_EXACT_OUTPUT/site_packets/RUSH/share"
python3 diagnostics/exact_aggregates/collect.py "$CLIF_EXACT_OUTPUT/reports/aggregate_collection" "$CLIF_EXACT_OUTPUT/site_packets/UCMC/share" "$CLIF_EXACT_OUTPUT/site_packets/NU/share" "$CLIF_EXACT_OUTPUT/site_packets/RUSH/share"
export CLIF_DESCRIPTOR_REPORTS="$CLIF_EXACT_OUTPUT/reports/publication_tables"
export CLIF_SURVIVAL_OUTPUT="$CLIF_EXACT_OUTPUT"
export CLIF_FIGURE_OUTPUT="$CLIF_EXACT_OUTPUT"
python3 diagnostics/manuscript_figures/plot_survival_tdauc.py
python3 diagnostics/manuscript_figures/refinements.py
