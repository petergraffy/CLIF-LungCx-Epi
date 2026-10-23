#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
case "${1:-derive}" in
  setup) Rscript code/00_setup.R ;;
  prepare) Rscript code/01_prepare.R ;;
  derive) Rscript code/02_derive.R ;;
  model)
    if [ "$#" -ne 2 ]; then echo "Usage: ./run.sh model path/to/consensus_mapping.csv" >&2; exit 1; fi
    Rscript code/03_model.R "$2"
    ;;
  test) Rscript tests/run_tests.R ;;
  *) echo "Usage: ./run.sh [setup|prepare|derive|model MAPPING|test]" >&2; exit 1 ;;
esac
