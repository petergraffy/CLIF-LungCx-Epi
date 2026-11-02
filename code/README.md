# Site and coordinator run instructions

All commands in this guide run from the **repository root**, one directory above `code/`.

## Sites start here

**Run Stage 1 only. Return its aggregate results, then STOP and wait for the
coordinator’s Stage 2 parameter bundle.** This is the same local discovery →
central fingerprints/consensus → frozen local assignment/outcomes workflow
rehearsed with UCMC, NU and Rush. Sites do not choose a shared phenotype count,
use another site's private data, or run coordinator commands.

| Who | What to run | What happens next |
|---|---|---|
| Each site, now | Environment setup, preflight, `stage1` | Return `stage1_return/`; stop |
| Coordinator | Review Stage 1 returns, choose parameters, `freeze` | Send the frozen `stage2_parameters/` bundle to sites |
| Each site, after that bundle | `stage2-fit` | Return `stage2_fit_return/` |
| Coordinator | Freeze held-out predictors and common survival grid | Send `stage2_evaluation/` |
| Each site | `stage2-evaluate` | Return `stage2_final_return/` |
| Coordinator | `collect` | Pool associations, competing risks, ROC and supplemental returns |

There are **two scientific stages**. Stage 2 includes a required coordinator
round trip: other-site prediction coefficients cannot be learned before sites
fit the frozen phenotypes, and exact pooled AJ curves require a common event-time
grid. Neither exchange shares patient identifiers or individual predictions.

## 1 Set up the same environment at every site

Use the coordinator-issued commit of this repository and **R 4.4.2** and
**Python 3.12**. Do not update packages or edit `renv.lock` at a site.
Run all commands from the repository root. The R setup restores exact package
versions into the ignored project library `renv/library/`; it does not upgrade
your global library or install R itself.

```sh
Rscript --vanilla scripts/setup.R
python3.12 -m venv .venv
.venv/bin/python -m pip install -r config/waterfall_requirements.txt
```

On Windows, use `py -3.12 -m venv .venv` and
`.venv\Scripts\python.exe -m pip install -r config\waterfall_requirements.txt`.
The `Rscript --vanilla ...` commands below work on Windows too. Set the Python
executable in the site config to the corresponding absolute path.

`renv.lock` pins the R packages and dependency closure. Python requirements pin
CLIFpy 2.1.0 and its waterfall dependencies. A compiler/system dependencies may
be needed if packages such as Arrow must build from source; inspect setup errors
rather than substituting package versions. [Environment details](../docs/environment.md).

## 2 Configure your site

Copy `config/site_template.json` to **`config/site.json`**. Change only:

- `site_name`: the coordinator-approved site ID, using letters, numbers, `_` or `-`;
- `repo`: absolute path to this checkout;
- `tables_path`: absolute path containing your CLIF parquet tables;
- `python`: absolute path to the Python virtual-environment executable;
- `analysis.run_id`: the coordinator-issued run label.

Keep all other scientific settings unchanged. The study period is 2018–2024.
The portable waterfall workflow currently requires **parquet** inputs. Required
tables: patient, hospitalization, ADT, hospital_diagnosis (including POA),
respiratory_support, vitals, labs, medication_admin_continuous and
patient_assessments. Code status and optional oncology data have explicit
availability reporting. Stage/histology/treatment/ECOG are not inferred.

Copy `config/site_metadata_template.csv` to `config/site_metadata.csv` and fill
in institutional/source metadata where known. Missing confirmations remain
explicit. Both site config files are git-ignored.

## 3 Run Stage 1 and return the results

```sh
Rscript --vanilla scripts/site.R preflight config/site.json
Rscript --vanilla scripts/site.R stage1 config/site.json
```

Preflight checks R/package versions, the Python waterfall runtime and required
table files. Stage 1 prepares the cohort, applies bounded CLIF waterfall support
while retaining source-evidence ARF, discovers local clusters, and assesses
stability. It runs the primary, eight inferential sensitivity discoveries and
the descriptive absorbing-event discovery. It also captures provisional local
characteristics, cohort flow, missing-support/ARF diagnostics and SOFA/APACHE
component/source completeness for the supplement. No shared phenotype assignment
or post-landmark association/prediction fitting occurs in Stage 1.

**Send this directory only**, using your institution’s approved transfer method:

```text
output/run_<SITE>_<RUN_ID>_federated/stage1_return/
```

The command prints the exact path. A complete return contains
`return_manifest.json`, `discovery/`, `descriptives/`, `cohort_audit/`,
`supplement/` and `environment/`. The coordinator uses `discovery/` to learn
phenotypes; clinical/outcome diagnostics are separate and do not enter matching.
Local cluster labels are provisional and cannot be pooled as shared phenotypes.

**STOP after Stage 1.** Preserve the checkout, config and **all** corresponding
`output/run_<SITE>_<RUN_ID>*` folders on the same machine: the prepared baseline,
waterfall panels, sequences and discovery caches are required for Stage 2.
Never send the whole run directory, `local/`, CLIF tables, RDS caches, logs or
individual predictions. Counts in the designated return are exact aggregates,
including small cells. The runner validates inventory, hashes and CSV headers;
institutional review of the aggregate return still applies.

Completed steps can be resumed with the same command and unchanged code/config.
Failed workers may leave partial files: inspect `local/steps/*.log`; if the worker
cannot resume, use a fresh coordinator-approved run ID. Do not delete caches or
mix files from different runs. A completion manifest is written only after all
required work succeeds.

## 4 Stage 2 only after the coordinator releases parameters

Do **not** use `config/active_analysis.json` or the old three-site results as your
Stage 2 inputs. The coordinator will send a new bundle derived from the actual
Stage 1 returns. Keep the original Stage 1 checkout and site config.

```sh
Rscript --vanilla scripts/site.R stage2-fit config/site.json /absolute/path/to/stage2_parameters
```

The runner checks that your exact Stage 1 return was used, verifies code and
lockfile compatibility, then applies the frozen full and leave-one-site-out
phenotype templates. Every support-sufficient stay gets the closest phenotype;
insufficient-support stays remain in the analysis. Sites fit four post-landmark
cause-specific Cox endpoints, reviewer sensitivities and severity/simple-summary
comparators, and export training aggregates plus relative event-time counts.

Return only the printed `stage2_fit_return/` directory. The coordinator then
sends a frozen predictor/shared-grid bundle. Run:

```sh
Rscript --vanilla scripts/site.R stage2-evaluate config/site.json /absolute/path/to/stage2_evaluation
```

Return only the printed **`stage2_final_return/`** directory. It includes the
complete model/supplement packet under `share/`, exact hourly/count exports,
pooled ROC score-bin counts, and survival/tdAUC summaries. Evaluation uses frozen
other-site predictors without held-out-site refitting. Risk-set exports include
200 patient-cluster bootstrap count replicates for pooled AJ/KM uncertainty.
Read [the site return guide](../docs/site_return_workflow.md) for output definitions.

## Coordinator only

Sites must not run these commands. Copy `config/coordinator_template.json` to
`config/coordinator.json` and list explicit Stage 1 return directories. Use the
same commit/environment as sites. **No historical `sites/` folder is scanned.**

```sh
Rscript --vanilla scripts/coordinator.R inspect config/coordinator.json
# Review fingerprints, selected local k/stability, missingness and descriptors.
# Set selected_k per analysis if desired; null selects the supported silhouette maximum.
Rscript --vanilla scripts/coordinator.R freeze config/coordinator.json
```

`freeze` learns outcome-free full and leave-one-site-out consensus templates
and writes `stage2_parameters/` under your configured coordinator output.
Send that entire bundle to all listed sites. You can select a supported shared k
in the config; unsupported choices fail. Never choose it using outcome performance.
At least three sites are needed so each held-out discovery still has two sites.
Changing selected parameters after freezing requires a new coordinator output.

After all sites return their Stage 2 fit aggregates:

```sh
Rscript --vanilla scripts/coordinator.R predictors config/coordinator.json /returns/A/stage2_fit_return /returns/B/stage2_fit_return /returns/C/stage2_fit_return
# Send the generated stage2_evaluation/ bundle to sites.
```

After all final returns:

```sh
Rscript --vanilla scripts/coordinator.R collect config/coordinator.json /returns/A/stage2_final_return /returns/B/stage2_final_return /returns/C/stage2_final_return
```

Extend the return lists for all participating sites; names are not hardcoded.
The coordinator only reads aggregate packets. The full primary consensus,
held-out predictors, common survival grid and site inventories are versioned and
hash-checked. Pooling uses REML/modified Hartung–Knapp associations and exact
risk/event count pooling for AJ curves. Set the coordinator config `python` if
your plotting environment is not `.venv/bin/python` (Windows: `.venv/Scripts/python.exe`).

## Analysis and verification

Primary: joint respiratory support–ARF sequences over the first 72 ICU hours,
among hospitalizations alive and in hospital at ICU +72 hours; outcomes start
then. Separate hospital death, hospice, composite and live non-hospice discharge
models follow through 28 post-landmark days. Absorbing death/discharge groups
are descriptive, never regressed on their defining outcomes. Unknown physiology
is not negative. HFNC/CPAP/NIPPV are distinct. Modified SOFA and adapted APACHE
are qualified adaptations; these are not validated standard APACHE IV comparisons.

[Protocol](../docs/analysis_protocol.md) ·
[Reviewer requirements](../docs/site_return_requirements.csv) ·
[Site return details](../docs/site_return_workflow.md) ·
[Two-stage deployment details](../docs/two_stage_run.md).

Developer synthetic checks (no clinical data):

```sh
bash run.sh test
Rscript --vanilla tests/deployment.R
python3 tests/arf_display.py
python3 tests/sankey.py
```

Earlier `code/02_derive.R` → raw-profile consensus/model commands and dated
`diagnostics/*/run.sh` files are comparison or rehearsal workflows, not site
instructions for this deployment. The submitted implementation is preserved in
[archive/2026-10-21](../archive/2026-10-21/README_ARCHIVE.md). Historical local figures
and `config/active_analysis.json` describe the three-site rehearsal, not parameters
that new sites should apply before the coordinator completes Stage 1 review.
