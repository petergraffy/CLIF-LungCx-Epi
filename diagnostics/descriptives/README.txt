Publication descriptor capture and aggregate collection

Automatic portable route: code/04_supplement.R captures overall and provisional
local tables; diagnostics/site_return/capture.R captures frozen shared Table2.
See docs/site_return_workflow.md. All exact records remain local.

Rehearsal-only extension of frozen three-site caches:
Rscript diagnostics/descriptives/site_capture.R SITE NEW_DESCRIPTOR_ROOT
Rscript diagnostics/descriptives/add_to_return.R SITE DESCRIPTOR_ROOT BASELINE_ROOT OLD_SHARE NEW_SHARE
Rscript diagnostics/descriptives/collect.R NEW_REPORT_ROOT SITE_SHARE [SITE_SHARE ...]

The extension retains parent packet hashes and does not refit models. Collect
reads only aggregates and verifies file hashes, consensus and denominators.
Tests: Rscript tests/publication_descriptives.R
Current run: output/run_descriptives_20261030_v2.
