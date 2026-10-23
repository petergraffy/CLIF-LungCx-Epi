#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
wf_out="${CLIF_WATERFALL_OUTPUT:-output/run_UCMC_waterfall_20261022}"
wf_python="${CLIF_WATERFALL_PYTHON:-${wf_out}/venv/bin/python}"
export CLIF_WATERFALL_OUTPUT="$wf_out"
if [[ -e "$wf_out/local/waterfall_derivation.rds" ]]; then
  echo "This comparison already exists. Select a fresh CLIF_WATERFALL_OUTPUT directory." >&2
  exit 1
fi
if [[ ! -x "$wf_python" ]]; then
  echo "Set CLIF_WATERFALL_PYTHON to an isolated Python environment with clifpy==2.1.0." >&2
  exit 1
fi
"$wf_python" diagnostics/waterfall/test_boundaries.py
Rscript diagnostics/waterfall/prepare.R
"$wf_python" diagnostics/waterfall/run_waterfall.py
Rscript diagnostics/waterfall/compare.R
Rscript diagnostics/waterfall/qa.R
Rscript diagnostics/waterfall/plot.R
