# Coordinator postprocessing

Run this from the repository root after collecting completed Stage 2 aggregate
returns. It accepts the available sites and uses the issued frozen bundle.
Provide the original return directories so that their manifests can be verified.
The command uses the existing R analysis environment plus Python 3.10 or later
with NumPy, Matplotlib, and python-docx. The cohort flow also uses DiagrammeR,
DiagrammeRsvg, and rsvg in R.

```sh
python3 postprocessing/run.py output/run_stage2_pool \
  "/path/to/site-a/stage2_return" "/path/to/site-b/stage2_return"
```

Supply one argument per site. Use a fresh output directory for each pool.
The script verifies inventory checksums, checks for patient-level columns, and
checks native and PyArrow fallback returns against their issued contracts. The
documented UCSF operational difference is flagged in the audit, and the script
also pools the primary associations using only returns that match their issued
contract. The folder is outside the site's frozen source inventory.

- `packet_audit.csv` and `site_counts.csv` record the return checks and Stage 2
  cohort sizes. `combined_cohorts.csv` compares Stage 1 and Stage 2 counts when
  both are available.
- `pooled_associations/` contains the primary, sensitivity, and exploratory
  interval Cox estimates. `publication_tables/` contains pooled descriptors,
  and `tables/` contains editable main manuscript tables.
- `survival/reports/` contains pooled Kaplan–Meier and Aalen–Johansen estimates
  and their risk sets. `prediction/` contains pooled ROC curves, time-dependent
  AUC, and paired AUC differences from local cross-validation.
- `figures/` contains the three main figures, cohort flow, calibration and net
  benefit, and the frozen 24- and 48-hour reference profiles.

Pooling uses only existing aggregate returns and leaves the frozen assignments
in place. The prediction estimates use locally cross-validated scores.
Keep the institutional returns, run folders, and assembled manuscript files
in the designated project storage.
