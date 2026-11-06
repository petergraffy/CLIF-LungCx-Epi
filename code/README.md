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
2. Obtain the **released** `stage2_parameters/` folder from the coordinator.
   It must contain `bundle_manifest.json`. Use the Git commit specified in that
   manifest. A candidate/`NOT_RELEASED` folder is not an executable bundle.
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

### A. Validate your Stage 2 setup

```sh
Rscript --vanilla scripts/site.R stage2-preflight config/config.json /absolute/path/to/stage2_parameters
```

This checks the released bundle, current code/environment, site admission,
approved original Stage 1 manifest, original settings, and the required local
caches. It prints the original Stage 1 directory and the new Stage 2 workspace.
Preflight does not fit models or alter the original run.

If your site is not in the bundle, stop and contact the coordinator. The coordinator-issued bundle must name your site before execution. Every
participating site is expected to complete all Stage 2 analyses. If a cache or
local discovery export is missing, recover the original local run; do not replace
labels or rebuild a manifest to make the check pass.

### B. Fit the frozen phenotypes and local outcomes

```sh
Rscript --vanilla scripts/site.R stage2-fit config/config.json /absolute/path/to/stage2_parameters
```

The runner creates an **isolated Stage 2 workspace**, copies and verifies the
original local discovery caches, and records their lineage. It preserves cached
instants and landmark membership. It does not rerun clustering, reinterpret
source timezones, or edit Stage 1 deployment hashes. The frozen mapping is applied
locally: every support-sufficient stay gets its closest phenotype, and
insufficient-support stays remain explicitly represented.

This pass fits separate post-landmark cause-specific Cox models for **death,
hospice, death or hospice, and live non-hospice discharge**, plus the specified
sensitivities and severity/simple-summary comparators. It returns aggregates for
held-out prediction and the shared survival grid.

**Return only the printed `stage2_fit_return/` directory, then STOP.** Its path is:

```text
output/run_<SITE>_<ORIGINAL_RUN_ID>_stage2_<RELEASE_ID>/stage2_fit_return/
```

Keep the entire Stage 2 workspace locally. Do not send `local/`, RDS files, source
tables, identifiers or individual predictions. Do not run the evaluation pass
until the coordinator sends the second bundle.

### C. Evaluate after the coordinator sends the second bundle

The coordinator pools the fit aggregates and returns `stage2_evaluation/` with
frozen predictors and the common event-time grid. Then run:

```sh
Rscript --vanilla scripts/site.R stage2-evaluate config/config.json /absolute/path/to/stage2_evaluation
```

Use the same Stage 2 checkout/config and retain both bundles. Return only the
printed **`stage2_final_return/`** directory:

```text
output/run_<SITE>_<ORIGINAL_RUN_ID>_stage2_<RELEASE_ID>/stage2_final_return/
```

The final packet includes model/supplement outputs, exact trajectory counts,
survival risk/event counts, ROC score-bin counts and time-dependent AUC returns.
Other-site predictors are evaluated without refitting on the held-out site.
The coordinator pools associations and curves after receiving these returns.
[Site return definitions](../docs/site_return_workflow.md).

Stage 2 has two local passes with a coordinator exchange between them:

| Who | Action | Return or next step |
|---|---|---|
| Site | `stage2-preflight`, then `stage2-fit` | Send `stage2_fit_return/`; wait |
| Coordinator | Freeze predictors and event-time grid | Send `stage2_evaluation/` |
| Site | `stage2-evaluate` | Send `stage2_final_return/` |
| Coordinator | Collect final packets | Pool associations, AJ curves and ROC outputs |

Completed Stage 2 steps resume with the same command, code, config and bundle.
Read failed worker logs in the **Stage 2** workspace. Do not delete completion
markers or mix bundles from different releases.

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


## Coordinator only

Sites must not run these commands. Copy `config/coordinator_template.json` to
`config/coordinator.json` and list explicit Stage 1 return directories. Use the
issued Stage 2 commit/environment. **No historical `sites/` folder is scanned.**

```sh
Rscript --vanilla scripts/coordinator.R inspect config/coordinator.json
# Review fingerprints, selected local k/stability, missingness and descriptors.
# Set selected_k per analysis if desired; null selects the supported silhouette maximum.
Rscript --vanilla scripts/coordinator.R freeze config/coordinator.json
```

`freeze` learns outcome-free full and leave-one-site-out consensus templates
and writes `stage2_parameters/` under your configured coordinator output.
Send that entire bundle to all participating sites. **All sites are expected to
participate in the primary analysis, every sensitivity, and the outcome analyses.**
Do not predesignate sites for exclusion because of outstanding file repairs.
Missing or damaged discovery files are repair items; the coordinator verifies
all required inputs before issuing the final bundle.

Leave-one-site-out models are validation folds: each participating site is left
out of training in turn and then evaluated. This does not remove any site from
the study or from the main full-cohort analysis.

`accepted_stage1` approvals bind the exact original return-manifest hash and
source/environment identity, with a documented reason for accepting a legacy
Stage 1 discovery snapshot. These scoped approvals do not rewrite source hashes
or accept corrupted discovery inputs. The current Stage 2 code/environment must
match the released bundle exactly. Source timestamp reinterpretation or cohort
rederivation requires a new reviewed Stage 1 run, not a Stage 2 migration.

The primary mapping uses exact aggregate fingerprints, equal site weighting and
four shared phenotypes. All participating sites contribute to all discovery
analyses; sensitivity mappings have their own selected parameters. Finalize the
participating-site roster and required inputs before issuing the bundle. The
remaining site still runs Stage 1 and then follows the same Stage 2 workflow.

 You can select a supported shared k
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
python3 tests/repository_hygiene.py
Rscript --vanilla tests/deployment.R
python3 tests/arf_display.py
python3 tests/sankey.py
```

Earlier `code/02_derive.R` → raw-profile consensus/model commands are retained
for comparison and synthetic tests, not site instructions for this deployment.
Retired rehearsal runners, local result pointers, site returns, exposome data
and old implementation snapshots are excluded from the shared tree. Use the
coordinator-issued bundle, never historical local outputs, for Stage 2.

### Source/timestamp audit only when requested

If requested by the coordinator, run the standalone `site_release_audit.R`
provided in the reconciliation tools package **from the original Stage 1
checkout**, using its original config and R library:

```sh
Rscript --vanilla /absolute/path/to/tools/site_release_audit.R config/config.json /absolute/path/to/tools/reference_contract.json
```

Keep the tools outside that original checkout. The audit writes a separate
`output/run_<SITE>_<RUN_ID>_release_audit/` directory and never changes the Stage 1
packet or caches. Return that audit directory when requested. It contains code
differences, package versions and timestamp **schema metadata**, with no patient
rows, source-table rows, private configs or RDS caches. It diagnoses the installed
environment without installing or updating packages.

A timezone encoded in an Arrow timestamp schema describes actual instants.
Naive timestamps require the timezone of the **source extract** and a reviewed
conversion/DST policy. Do not infer that timezone from the hospital's location,
the current R session, or the computer's timezone. Changes to cohort/timestamp
code require source review; an unexplained hash cannot be reconciled by changing
its value. Submit the actual changed source with the audit.

The coordinator can accept a source-hash ordering difference only when every
named source-file hash matches and the submitted aggregate hash reproduces from
the submitted order. New source inventories use locale-independent ordering.
The missing-package setup guard affects installation checks; it does not relax
the lockfile or change the analysis. Arrow/R version differences still require
review independently of that guard.

Candidate numerical templates are not a release. Sites must wait for a complete
validated `stage2_parameters/` bundle with its production `bundle_manifest.json`;
folders labeled `candidate` or `NOT_RELEASED` cannot be used for `stage2-fit`.
All sites are expected to contribute to every analysis. Missing descriptors do
not invalidate a verified primary fingerprint, but required discovery files must
be repaired before final release. Pending files remain documented and must never
be made to pass by regenerating a manifest for altered results.
