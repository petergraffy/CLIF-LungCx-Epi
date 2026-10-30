#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
case "${1:-derive}" in
  setup) Rscript code/00_setup.R ;;
  prepare) Rscript code/01_prepare.R ;;
  derive) Rscript code/02_derive.R ;;
  supplement) Rscript code/04_supplement.R ;;
  model)
    if [ "$#" -ne 2 ]; then echo "Usage: ./run.sh model path/to/consensus_mapping.csv" >&2; exit 1; fi
    Rscript code/03_model.R "$2"
    ;;
  test)
    Rscript tests/run_tests.R
    Rscript tests/apache_adapted.R
    Rscript tests/sofa_recovery.R
    Rscript tests/site_supplement.R
    Rscript tests/publication_descriptives.R
    Rscript tests/competing_event_curves.R
    Rscript tests/survival_tdauc.R
    Rscript tests/exact_aggregates.R
    Rscript tests/pooled_survival.R
    Rscript tests/pooled_roc.R
    Rscript tests/nearest_assignment.R
    Rscript tests/site_return_collection.R
    Rscript tests/model_supplement_diagnostics.R
    ;;
  *) echo "Usage: ./run.sh [setup|prepare|derive|supplement|model MAPPING|test]" >&2; exit 1 ;;
esac
