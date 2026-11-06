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

Stage 2 needs only a site config, the source CLIF tables, pinned environment and
tracked frozen bundle. It rebuilds all local input data without clustering.
`stage1_root` is ignored; no original Stage 1 run ID or output folders are needed.
Use a fresh `analysis.run_id` if starting another Stage 2 run or changing inputs.
The runner creates an isolated Stage 2 workspace and records the input-file
sizes/mtimes and frozen timestamp policy for safe resume.

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
