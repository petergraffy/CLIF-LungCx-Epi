# Site and coordinator run instructions

Run every command from the **repository root**, one directory above `code/`.

## Stage 2 — start here if you already completed Stage 1

**Sites with completed Stage 1 should now prepare for Stage 2 using the
coordinator-issued frozen bundle and its specified code commit.** Stage 2 applies
the shared phenotypes locally and analyzes post-landmark outcomes. Patient-level
records and identifiers remain at your site.

**Still need to run Stage 1?** Go to [Stage 1 instructions below](#stage-1--for-a-site-that-has-not-yet-returned-results).
Complete Stage 1, return its packet, and wait for coordinator admission to Stage 2.
The coordinator will include the remaining site in the study workflow; sites do
not learn or change the shared mapping themselves.

### Before running Stage 2

1. Preserve the original Stage 1 checkout, private config and **all** original
   `output/run_<SITE>_<RUN_ID>*` folders. Do not delete caches, edit manifests,
   change the original run ID, or rerun Stage 1 just to start Stage 2.
2. Use the **released** `frozen/stage2_v1/` folder included in your clone.
   It must contain `bundle_manifest.json`. Use the released checkout containing
   that bundle; the manifest enforces its analytical source and R-lock hashes. A candidate/`NOT_RELEASED` folder is not an executable bundle.
3. Use a separate Stage 2 checkout if the original checkout contains local fixes.
   Copy the original config into its ignored `config/config.json`, set `repo` to
   the Stage 2 checkout, and add `stage1_root` pointing to the original federated
   run directory. Keep `site_name`, `analysis.run_id` and scientific settings
   exactly as in Stage 1. The site ID is case-sensitive.

Example config fields (retain all the other original fields):

```json
{
  "repo": "/absolute/path/to/STAGE2_CHECKOUT",
  "stage1_root": "/absolute/path/to/ORIGINAL_CHECKOUT/output/run_SITE_ORIGINAL_RUN_ID_federated"
}
```

`stage1_root` is a directory, not the `stage1_return/` folder. It is optional when
using the same checkout and the original run remains at its default location.
Keep the original source tables and waterfall caches accessible at their original
paths; the prepared caches retain those local references.

Use **R 4.4.2** and **Python 3.12**. Restore the pinned R packages in the Stage 2
checkout; `renv` does not install R itself:

```sh
Rscript --vanilla scripts/setup.R
```

You may reuse the original pinned Python environment by retaining its absolute
`python` config path. For a new Python environment:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r config/waterfall_requirements.txt
```

On Windows use `py -3.12 -m venv .venv` and
`.venv\Scripts\python.exe -m pip install -r config\waterfall_requirements.txt`.
Set `python` in your private config to that executable. Do not substitute package
versions or update `renv.lock`. [Environment details](../docs/environment.md).

### Run Stage 2 once and return one packet

The frozen bundle is included in your clone at **`frozen/stage2_v1/`**. Use the
released code checkout and restore the R/Python environments above. Then run:

```sh
Rscript --vanilla scripts/site.R stage2-preflight config/config.json
Rscript --vanilla scripts/site.R stage2 config/config.json
```

The default bundle is `frozen/stage2_v1`. An explicit third argument may select a
later coordinator-issued release. Preflight checks bundle hashes, site admission,
original Stage 1 lineage, settings and caches without fitting models.

Stage 2 creates an isolated workspace, preserves the original cached instants and
cohort membership, applies the frozen full and leave-one-site-out phenotypes,
and runs **all** local modeling, sensitivities, diagnostics and exports.
Every support-sufficient stay receives its closest phenotype; insufficient
support remains explicitly represented. Death, hospice, death or hospice, and
live non-hospice discharge each have separate cause-specific Cox models.

**Return only the printed `stage2_return/` directory:**

```text
output/run_<SITE>_<ORIGINAL_RUN_ID>_stage2_<RELEASE_ID>/stage2_return/
```

That single packet includes Table 1 and phenotype Table 2 descriptors, severity
components and missingness, all model coefficients and diagnostics, exact hourly
support/ARF counts, Sankey counts, phenotype transport agreement, survival event
and censor counts, and prediction evaluation aggregates. Retain the complete
workspace locally. No patient records, identifiers, RDS caches or individual
predictions should be uploaded. **No second site pass or second bundle is needed.**

Prediction uses five-fold **local patient-grouped cross-validation**: all
hospitalizations from one patient stay in one fold. The 28-day out-of-fold risk
markers are evaluated over days 1–28; ROC curves are pooled at days 7, 14 and 28.
Baseline, baseline + measured-component SOFA, baseline + adapted APACHE and
baseline + phenotype are compared on the same complete-score cohort. Summary
comparators and all-case missing-indicator models are retained. Empty or
nonestimable comparisons return explicit status/denominator records.

These are internally validated local models. Leave-one-site-out templates test
**phenotype transportability**, not external validation of an outcome predictor
trained at other sites. Bootstrap uncertainty resamples patients and their fixed
out-of-fold predictions, conditional on the frozen phenotypes and fitted folds.
Adapted APACHE is neither standard APACHE II nor APACHE IV.

The coordinator reconstructs exact pooled KM/Aalen–Johansen curves from aggregate
relative-duration event and censor histories, including patient-cluster bootstrap
counts; no shared event-grid exchange with sites is required.

Completed steps resume with the same command, code, config and bundle. Read
failed worker logs in the Stage 2 workspace. Do not mix bundles or edit completion
markers. A site missing from the bundle should contact the coordinator for an
admission update; the shared prototypes remain frozen.

## Stage 1 — for a site that has not yet returned results

**The Stage 1 workflow remains available for the remaining site.** Use the
coordinator-issued Stage 1 commit and run label, follow the steps below, and return
`stage1_return/`. Stop there. The coordinator will tell you which released bundle
to use for Stage 2 after reviewing your return; do not run coordinator commands.

### 1. Set up the Stage 1 environment

Use a **Git checkout** of the coordinator-issued commit (Git must be installed),
**R 4.4.2**, and
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

### 2. Configure your site

Copy `config/config_template.json` to **`config/config.json`**. Change only:

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

### 3. Run Stage 1 and return the results

```sh
Rscript --vanilla scripts/site.R preflight config/config.json
Rscript --vanilla scripts/site.R stage1 config/config.json
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

If `07_descriptors_supplement` reports `sh: ~/.../audit.log: No such file or
directory`, an older runner passed an unexpanded home-directory path to shell
redirection. The corrected runner expands configuration paths and canonicalizes
the audit output before starting child processes, including paths with spaces.
As a workaround on the original checkout, use full absolute paths in the site
config (for example `/Users/YOUR_ACCOUNT/Desktop/CLIF_Projects/CLIF-LungCx-Epi`
for `repo`) and a fresh coordinator-issued run ID, then repeat preflight and
Stage 1. Preserve the previous run folders: a failed descriptor step may already
have written an immutable partial package, so simply repeating the old run ID
may fail with `Descriptor package exists`. Updating source also requires a fresh
run ID under the code-contract checks; do not edit deployment hashes or markers.

### Stage 1 refresh only when requested by the coordinator

**This refresh procedure is not the Stage 2 setup. Sites with completed Stage 1
should use the Stage 2 instructions above and preserve their original run ID.**

The corrected code accepts the UCUM blood-gas unit `mm[Hg]` as `mmHg` for the
adapted APS-plus-age comparator. It does not substitute normal values for absent
measurements. Each site now returns `supplement/apache_adapted_lab_units.csv`,
which counts first-24-hour lab rows by source unit and acceptance/rejection.
Review this alongside component completeness before interpreting a missing score.

Source and R lockfile contracts normalize Windows CRLF to LF. All other text
differences remain detectable. Returns include per-file source hashes so the
coordinator can identify changed files. Uploaded artifacts retain exact byte-level
checksums; empty files are rejected. The coordinator does not waive older hash
mismatches or rewrite their manifests.

For sites that already ran Stage 1, preserve the original checkout/config and all
run folders. Obtain the same corrected checkout from the coordinator, keep the
issued `renv.lock`, and set a **new common coordinator-issued run ID** in each
site's private config (for example, `federated_20261103_v2`). Then run:

```sh
Rscript --vanilla scripts/setup.R
Rscript --vanilla scripts/site.R preflight config/config.json
Rscript --vanilla scripts/site.R stage1 config/config.json
```

Return the new complete `stage1_return/` folder. This refresh regenerates local
severity components and discovery diagnostics under one code contract; the
coordinator must use the refreshed returns together. Simply rerunning the old
completed run ID cannot refresh scores. Do not edit hashes, remove step markers,
replace individual result files, or combine old and refreshed packets.

Sites must investigate truly missing PaO2, bilirubin or discharge disposition in
their source/ETL locally. Unknown measurements and dispositions remain explicitly
unknown; a central aggregate upload cannot recover them. Preserve their counts
in the supplement and cohort reconciliation.


## Coordinator workflow

Inspect the original Stage 1 returns, choose supported phenotype counts and record
explicit legacy-source acceptance in a private coordinator config. All sites
participate in Stage 2. Any missing discovery contribution must be named and
justified in `discovery_exceptions`; it does not remove that site from modeling.

```sh
Rscript --vanilla scripts/coordinator.R inspect config/coordinator.json
Rscript --vanilla scripts/coordinator.R freeze config/coordinator.json
```

`freeze` learns outcome-free full and leave-one-site-out templates for every
variant, copies the fixed analysis plan, and writes a checksummed
`stage2_parameters/` bundle. Publish those aggregate CSV/JSON templates under
`frozen/stage2_v1/` in Git; source tables, institutional uploads, private configs
and patient caches remain outside Git. The shipped release records exactly which
sites contributed to each discovery, including explicitly accepted omissions.

Sites run `stage2` once and upload `stage2_return/`. After every participating site
returns, collect all packets with the same private coordinator config:

```sh
Rscript --vanilla scripts/coordinator.R collect config/coordinator.json /path/SITE_A/stage2_return /path/SITE_B/stage2_return /path/SITE_C/stage2_return
```

The coordinator checks the release identity and source hashes, pools association
estimates and prediction histograms, reconstructs pooled survival/competing-risk
curves, and retains all site tables for the supplement. No outcome-predictor
bundle is sent back to sites.

A late site completes the retained Stage 1 workflow first. Create an admission revision without relearning the full prototypes:

```sh
Rscript --vanilla scripts/admit_stage2_site.R frozen/stage2_v1 /path/LATE_SITE/stage1_return frozen/stage2_v1_admission2 "Reviewed Stage1 discovery; accept original cached cohort and instants."
```

Issue that revision to all sites before their runs. The late site uses the full frozen templates for its
external phenotype transport assessment because it contributed no discovery data.
