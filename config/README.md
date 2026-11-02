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

The end date is exclusive: `2025-01-01` includes all of 2024. Change a run ID
when inputs or scientific parameters change. Sites must not modify the
coordinator-issued analysis settings or reuse caches under a different contract.
