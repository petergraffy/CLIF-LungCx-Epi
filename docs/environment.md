# Reproducible site environment

Use R 4.4.2, Python 3.12 and the coordinator-issued source commit at every stage.
R versions are not installed by renv. Run from the repository root:

```sh
Rscript --vanilla scripts/setup.R
python3.12 -m venv .venv
.venv/bin/python -m pip install -r config/waterfall_requirements.txt
Rscript --vanilla scripts/site.R preflight config/config.json
```

`setup.R` bootstraps renv 1.1.5 if needed and restores `renv.lock` into
`renv/library/`. Setup changes project libraries, not the global library.
The site/coordinator entry points explicitly add the project library and verify
all locked package versions. Child R workers inherit that library through
R_LIBS. Thus --vanilla does not bypass the project's package selection.
Sites do not call renv::snapshot(), update.packages(), or install latest package
versions to bypass an error. Reopen/restart R if packages were already loaded
before setup. Native RStudio users should run the documented Rscript commands
in the Terminal, which avoids loaded namespaces and personal startup settings.

The lockfile includes direct and transitive R dependencies, including Arrow,
fst, TraMineR, survival, comorbidity, metafor, ggplot2 and timeROC. CRAN access,
archived package access and sometimes native build tools/system libraries are
needed on a clean machine. Locked package versions may still require an Arrow
source build on a platform without matching binaries; retain the version and
resolve system-library/compiler requirements. renv library/cache contents are
ignored; only the lockfile and restoration scripts belong in Git.

Python requirements capture the dependency closure used by official CLIFpy
2.1.0 waterfall, excluding unrelated manuscript/agent tooling. The runtime
check verifies Python major/minor, pinned package versions and records the
waterfall function source hash. The existing upstream Polars runtime check can
fall back to DuckDB; the actual waterfall code remains the pinned upstream
function. Do not replace it with custom imputation.

A Stage 1 return includes a copy of the R lockfile, R session information,
Python requirements and runtime-check output. Its manifest also records source
and lockfile hashes. Stage 2 rejects returns/bundles from a different source
or environment contract. Source CLIF tables and numeric local diagnostics are
never part of the environment return.
