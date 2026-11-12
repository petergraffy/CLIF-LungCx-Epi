# Site run instructions

**Stage 2 is ready.** It rebuilds the analysis data from your CLIF tables and
applies the frozen phenotypes. **You do not need Stage 1 folders or results.**

Run all commands from the **repository root** (the folder containing `scripts/`
and `renv.lock`). Use the coordinator-issued repository release.

## Stage 2

### 1. Set up the environment

Use **R 4.4.2**, **Python 3.12**, and **Git**. Restore the pinned R environment:

```sh
Rscript --vanilla scripts/setup.R
```

This restores the project R packages using `renv.lock`; it does not install R.

Create the Python environment, or reuse an existing environment with the same
pinned packages.

**macOS / Linux:**

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r config/waterfall_requirements.txt
```

**Windows:**

```sh
py -3.12 -m venv .venv
.venv\Scripts\python.exe -m pip install -r config\waterfall_requirements.txt
```

Keep the pinned versions. See [environment help](../docs/environment.md) if setup fails.

### 2. Configure your site

Copy [config/config_template.json](../config/config_template.json) to
**`config/config.json`**. Edit only these fields:

| Field | What to enter |
|---|---|
| `site_name` | Your coordinator-assigned site ID, with the exact capitalization |
| `repo` | Absolute path to this repository checkout |
| `tables_path` | Absolute path to your CLIF parquet tables |
| `python` | Absolute path to the Python executable in your pinned environment |

The run name is already set to **`stage2_v1`** in the template; leave it unchanged.
Your site ID keeps each site’s output folder separate. Existing configs from a
completed run can keep their original run name.

Keep all other settings unchanged. In Windows JSON paths, use `/` or escaped
backslashes (`\\`). Your config stays local and is excluded from Git.

The required CLIF tables are: patient, hospitalization, ADT, hospital_diagnosis
(with POA flags), respiratory_support, vitals, labs, medication_admin_continuous,
and patient_assessments.

### 3. Run Stage 2

```sh
Rscript --vanilla scripts/site.R stage2-preflight config/config.json
Rscript --vanilla scripts/site.R stage2 config/config.json
```

Preflight checks your setup. The second command rebuilds the local data, applies
the frozen mappings, and runs all models, sensitivities, diagnostics and exports.
It uses the **`frozen/stage2_v1/` bundle already included in the repository**.
No local clustering is performed.

### OHSU Stage 2

Set `site_name` to `OHSU` in `config/config.json`. Use the same pinned
environment and the OHSU bundle included in this checkout:

```sh
Rscript --vanilla scripts/site.R stage2-preflight config/config.json frozen/stage2_v1_ohsu
Rscript --vanilla scripts/site.R stage2 config/config.json frozen/stage2_v1_ohsu
```

This bundle runs directly from OHSU's CLIF tables using the original frozen
phenotype prototypes. Return only the `stage2_return` folder printed by the
runner. If R Arrow is unavailable at OHSU, contact the coordinator for an
OHSU-compatible fallback bundle.

### Alternative: Python fallback when R Arrow fails

Use this separate route if R Arrow will not install, lacks Dataset/Parquet
support, or your machine has a different R Arrow version. It uses the same
frozen phenotypes and runs the complete Stage 2 analysis and supplement.
**R Arrow is not installed, loaded, or version-checked by this route.** An
installed Arrow 25 can remain installed. Keep the repository `renv.lock`
unchanged.

Use the same `config/config.json` and pinned Python environment described
above. Run from the repository root:

```sh
Rscript --vanilla python_fallback/setup.R
Rscript --vanilla python_fallback/site.R stage2-preflight config/config.json
Rscript --vanilla python_fallback/site.R stage2 config/config.json
```

Use `python_fallback/setup.R` for this route. The standard `scripts/setup.R`
restores R Arrow from the lockfile.

The fallback rebuilds from CLIF with no Stage 1 folders and no clustering.
Its fixed run name is **`stage2_python_v1`**. It creates its own workspace;
rerun the same fallback command to resume. It does not resume a partially
completed standard R Arrow run.

When it prints **`PYTHON FALLBACK COMPLETE`**, upload only:

```text
output/run_<SITE>_stage2_python_v1/analysis/stage2_return/
```

The packet contains the same model, sensitivity, descriptor, supplement,
survival, ROC and tdAUC export families, with explicit PyArrow provenance.
Do not upload `local/pipeline/` or the whole run folder.

[Fallback details and coordinator collection](../python_fallback/README.md)

### 4. Return the results

When the command prints **`STAGE 2 COMPLETE`**, send the coordinator **only the
printed `stage2_return/` folder**, using the agreed transfer method:

```text
output/run_<SITE>_<RUN_ID>_stage2_<RELEASE_ID>/stage2_return/
```

**One run, one return. No second site pass is needed.** The packet contains the
aggregate results needed for the manuscript and supplement. Patient records,
identifiers, individual predictions and RDS caches remain local; do not send the
whole run directory. Keep the new Stage 2 workspace until your return is accepted.

## Stage 1 — only if you have not submitted it yet

Use the coordinator-issued Stage 1 checkout and run label. Follow **steps 1–2
above** to set up the environment and config, setting `analysis.run_id` to the
Stage 1 label supplied by the coordinator. Then run:

```sh
Rscript --vanilla scripts/site.R preflight config/config.json
Rscript --vanilla scripts/site.R stage1 config/config.json
```

When the command prints **`STAGE 1 COMPLETE`**, send **only**:

```text
output/run_<SITE>_<RUN_ID>_federated/stage1_return/
```

**Stop and wait for the coordinator to admit your site to Stage 2.** Keep the
submitted aggregate packet for your records. Stage 1 output folders are not
required to run Stage 2.

## If a command fails

Open the log file named in the error message. Contact the coordinator with the
error if you need help; do not upload a partial return.

To resume, rerun the same command with the same config, code and bundle. If the
runner reports changed inputs or code, contact the coordinator for restart
instructions. Do not rename a completed run, edit hashes or delete completion
markers to bypass an error. A site-admission error requires a bundle
update from the coordinator.

[Coordinator instructions and workflow details](../docs/two_stage_run.md)
