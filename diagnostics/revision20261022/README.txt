Reviewer sensitivity and severity-comparison workers

This directory name records the implementation origin; its reusable workers
are used by the portable two-stage workflow. Follow code/README.md. Do not
invoke a historical three-site driver or load another site's clinical cache.

Stage 1 discovers h24, h48, onset, support_only, low_flow, carry0, carry12 and
raw_support sequences, plus descriptive absorbing-event sequences. The
coordinator freezes a separate full/held-out consensus for each variant.

Stage 2 fits post-landmark associations, code-status and ARF characterization,
COVID/full-code/first-hospitalization subgroups, severity and simple-summary
comparators, and held-out prediction. Unknown evidence remains explicit.
Absorbing groups are descriptive and are not regressed on their defining exits.

The supplied site list controls reusable workers. Source-label and code-status
vocabulary audits report unavailable optional fields without requiring a manual
site metadata form. Scientific definitions are in docs/analysis_protocol.md.
Local clinical caches and generated reports are excluded from Git. Only the
runner's designated aggregate return folders are transferred.
