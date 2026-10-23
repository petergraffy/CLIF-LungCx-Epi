# Site configuration

Copy `config_template.json` to `config.json` and set `site_name`, `repo`, and `tables_path`. The loader identifies `.csv`, `.parquet`, and `.fst` files directly; `file_type` is retained only for compatibility with site configuration tools.

Do not commit `config.json`. It is excluded by this directory's `.gitignore`.

Analysis defaults are defined in `code/lib/common.R`. Keep the scientific settings consistent across sites. The end date is exclusive, so `2025-01-01` includes all of 2024. Select a new `analysis.run_id` for any change to the analysis parameters or source code. Existing local caches and model exports cannot be silently reused under a different specification.

Optional `analysis.consensus_k` is a positive integer overriding the median local k at the coordinator; normally omit it. Set `min_cell` to the approved site disclosure threshold. The default is ten. Missing observations and suppressed cells remain missing in central analysis.

Additional optional inputs are `clif_code_status` and `clif_oncology` (also accepted without the `clif_` prefix). Oncology has one row per hospitalization and may include `histology`, `stage`, `ecog`, `therapy_line`, and `therapy_intent`; these are descriptive until their timing and coding are validated.
