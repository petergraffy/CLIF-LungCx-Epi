# Analysis configuration

Use one site configuration: copy `config_template.json` to `config.json`.
The complete setup and run instructions are in [the code README](../code/README.md).
Run commands from the repository root, not this directory.

Set `site_name`, `repo`, `tables_path`, `python`, and the coordinator-issued
`analysis.run_id`. Use absolute paths and retain the template's scientific
settings. The adopted two-stage waterfall workflow requires CLIF parquet tables.
The earlier comparison loader also supports CSV/FST; it is not the deployment
entry point.

`config.json` is private and ignored by Git. There is no separate site template
or metadata form. Administrative/reporting information is maintained outside
the analysis repository. Automatic environment and source-contract hashes still
accompany returns.

`coordinator_template.json` is a different configuration for the coordinator:
it lists Stage 1 return folders and shared-cluster choices. Sites do not use it.
Copy it to the ignored `coordinator.json` only when coordinating a run.

For Stage 2, retain the original `site_name`, `analysis.run_id` and analysis
settings. If using a new checkout, set `repo` to that checkout and `stage1_root`
to the original `output/run_<SITE>_<RUN_ID>_federated/` directory. `stage1_root`
is optional for Stage 1 and when the default original run location is available.
The runner creates a separate Stage 2 workspace from the bundle's release ID;
do not change the Stage 1 run ID to create that workspace.

All participating sites are expected to contribute to all analyses. Explicit `accepted_stage1`
records bind approved legacy discovery snapshots to their original manifest
hashes. They do not waive checksums for files used to learn templates. Sites
should use only their coordinator-issued released bundle.

The end date is exclusive: `2025-01-01` includes all of 2024. Change a run ID
when inputs or scientific parameters change. Sites must not modify the
coordinator-issued analysis settings or reuse caches under a different contract.

Stage 2 is one site pass: `scripts/site.R stage2 config/config.json` uses the
tracked `frozen/stage2_v1/` bundle by default and produces one `stage2_return/`.
`stage2_analysis.json` is the coordinator's fixed specification, copied into the
bundle at freeze time; sites must not edit it. A coordinator may explicitly name
missing discovery contributions in `discovery_exceptions`; every site still runs
all Stage 2 analyses under the shared frozen mappings.
