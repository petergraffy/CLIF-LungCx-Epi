# Python Parquet fallback for Stage 2

Sites should follow the [site instructions](../code/README.md). This separate
pipeline replaces R Arrow Parquet I/O with the pinned Python environment's
PyArrow. R still performs the existing analysis. It requires R 4.4.2 and the
same locked R dependencies except `arrow`, plus Python 3.12 and the existing
`config/waterfall_requirements.txt` packages. R Arrow can be absent or any
version, including 25. Do not edit `renv.lock`.

```sh
Rscript --vanilla python_fallback/setup.R
Rscript --vanilla python_fallback/site.R stage2-preflight config/config.json
Rscript --vanilla python_fallback/site.R stage2 config/config.json
```

Return only `output/run_<SITE>_stage2_python_v1/analysis/stage2_return/` after
`PYTHON FALLBACK COMPLETE`. No Stage 1 folders are needed. The run includes
all original Stage 2 models and sensitivities, full/LOO phenotype assignment,
local validation, survival histories, ROC/tdAUC inputs, SOFA/APACHE audits,
Table 1/2 descriptors, and supplemental exports.

## How isolation and provenance work

The runner verifies the original frozen bundle and source contract before
building a local analysis copy in `output/run_<SITE>_stage2_python_v1/local/pipeline/`.
Explicit versioned overrides change only I/O and environment checks. They
cover the cohort reader, timestamp schema audit, waterfall input/output,
code-status catalogues, and all SOFA supplement readers. Production source
files, existing native run folders, and `frozen/stage2_v1/` remain unchanged.
There is no runtime namespace patching or evaluation of edited source text.

The generated bundle contains identical numerical templates and analysis
plans. Its distinct source hash and manifest identify the PyArrow pipeline;
it also records the original frozen bundle's hash. All generated source is
checked against `release.json`. Native and fallback return manifests retain
their actual source and environment contracts.

The reader filters identifiers and projects columns in Python before loading
clinical data into R. Text and integers use escaped CSV; floating-point and
timestamp values use binary doubles and separate null masks. This preserves
large signed identifiers, empty strings, nulls, NaN and numerical precision.
Temporary clinical transport files stay in a private local directory and
are removed after each request. Nested/nonstandard projected field types
fail explicitly. No patient records are sent to the coordinator.

Preflight tests actual Python Dataset/Parquet reading and writing as well as
the pinned environment. A partially completed native run is not reused.
Use the same fallback command/config to resume its completed steps. Preserve
the generated workspace until your return is accepted. If the fallback
release or configuration changes, use a fresh `CLIF_FALLBACK_OUTPUT` directory
as directed by the coordinator; do not edit completion markers.

## Coordinator: collect native and fallback packets together

When any sites used the Python fallback, use the compatible collector for
the complete roster, passing one return per site:

```sh
Rscript --vanilla python_fallback/collect.R frozen/stage2_v1 output/run_pooled_stage2_v1 .venv/bin/python /path/to/site1/stage2_return /path/to/site2/stage2_return
```

Continue the command with every remaining site's return path. The collector
requires the frozen roster; it rejects duplicate or missing sites. It
validates every packet inventory and accepts only the exact native contract
or the issued fallback contract tied to the same parent bundle, lockfile,
and consensus. It does not rewrite the received manifests. Pooling methods
and the aggregate output layout are the existing Stage 2 implementation.
The command pools associations and prediction histograms and reconstructs
pooled AJ/KM curves, ROC and tdAUC outputs. The collection manifest records
each site's Parquet backend. The standard
native collector remains unchanged and deliberately rejects fallback packets.

For association tables only, `python_fallback/central_collect.R` takes
`BUNDLE NEW_OUTPUT SITE_RETURN [...]`; the wrapper above also runs the
existing central survival and prediction programs.

## Validation and existing input semantics

The isolated three-site trial compared 30 cohort-filtered source tables,
nine raw and ten waterfall analysis panels per site, frozen primary assignments,
and 12 primary Cox fits per site. The final reader uses the same transport.

The fallback preserves the current source-time representation and code-status
tie handling rather than changing the scientific analysis. In the trial,
one Rush stay with conflicting simultaneous orders had a different landmark
code-status label in three sensitivities because of scan order; baseline
code status, all hourly panels and primary results matched. Timezone-naive
conversion under a non-UTC source policy can depend on the host timezone in
the existing analysis. The fallback return records the host timezone. These
are existing analysis portability issues, not reasons to alter frozen
phenotypes inside a site run.

## Maintaining the release

Coordinator only: after reviewing changes to an override or launcher, run
`Rscript --vanilla python_fallback/make_release.R` and validate a fresh
workspace. Commit `release.json` alongside the fallback. The generator refuses
to bind to a modified original source contract. Do not regenerate a release
locally at a participating site to bypass an error.
