## Code Directory

This directory contains the scripts used to run the CLIF lung cancer ICU trajectory analysis. Each site should run the scripts locally in the order listed below.

### Workflow

1. **Restore the R environment**

Run `00_renv_restore.R`

This script restores the project environment using **renv** and installs all required R packages.

---

2. **Build the lung cancer ICU cohort**

Run `01_lungcx_cohort.R`

This script:

- Identifies ICU admissions
- Applies lung cancer diagnosis and resection criteria
- Defines the ICU index time (`t0`)
- Constructs ARF physiologic flags
- Assigns air pollution exposure
- Produces the analysis-ready cohort dataset

Outputs include:

- `cohort_lung`
- `analysis_ready`

---

3. **Run federated clustering and modeling**

Run `02_federated_clusters.R`

This script:

- Builds 72-hour respiratory trajectories
- Performs sequence clustering to identify trajectory phenotypes
- Characterizes clusters using severity and comorbidity metrics
- Fits pollution exposure models
- Generates PHI-safe federated outputs

---

### Output

Final results are written to:

`output/final`

These outputs include cluster summaries, model results, and figures suitable for federated pooling across CLIF sites.