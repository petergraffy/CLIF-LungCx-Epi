# CLIF Lung Cancer ICU Trajectories and Air Pollution

## Overview

This repository contains a federated analysis pipeline for studying respiratory failure trajectories, ICU outcomes, and long-term ambient air pollution exposure among patients with lung cancer admitted to the ICU across CLIF sites.

The project is designed so that each site runs the same analytic code locally and returns only PHI-safe, aggregated outputs for central pooling. No patient-level data are transferred.

## Main Aims

The current workflow supports three linked aims:

1. Define a multicenter ICU cohort of hospitalizations with primary lung cancer present on admission.
2. Use unsupervised state sequence clustering to identify distinct early respiratory support and acute respiratory failure trajectory phenotypes during the first 72 hours of ICU care.
3. Evaluate whether long-term PM2.5 and NO2 exposure are associated with trajectory severity, trajectory phenotype, and key ICU outcomes.

The project also produces descriptive epidemiology outputs to support a manuscript focused on lung cancer ICU epidemiology, including baseline cohort summaries, cluster-specific clinical profiles, early support patterns, and year-level trend summaries.

## High-Level Workflow

Sites should run the scripts in [code](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/code) in this order:

1. [00_renv_restore.R](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/code/00_renv_restore.R)
2. [01_lungcx_cohort.R](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/code/01_lungcx_cohort.R)
3. [02_federated_clusters.R](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/code/02_federated_clusters.R)

### Step 1: Environment restore

`00_renv_restore.R` restores the project R environment and installs required packages.

### Step 2: Cohort construction

`01_lungcx_cohort.R` builds the analysis cohort and creates the core objects used downstream:

- `cohort_lung`
- `analysis_ready`
- `flow_lung`
- `exclusion_lung`

This step:

- identifies adult hospitalizations with an ICU stay
- restricts to hospitalizations with primary lung cancer present on admission
- assigns the ICU index time
- links long-term county-level PM2.5 and NO2 exposure
- constructs cohort-level outcomes and baseline variables

### Step 3: Federated clustering and modeling

`02_federated_clusters.R` performs local trajectory construction, clustering, site-level modeling, and federated export generation.

This step:

- reconstructs hourly respiratory support and ARF states over the first 72 hours
- computes optimal matching distances with `TraMineR`
- selects the number of clusters using peak average silhouette width
- reorders local clusters so local `Cluster 1` is the healthiest reference cluster
- summarizes clusters using SOFA, Charlson, vasopressors, and cancer severity proxies
- fits site-level pollution and outcome models
- exports aggregated tables, coefficients, response-curve summaries, and site figures

## Cohort Definition

The intended analytic cohort is:

- adult hospitalizations
- lung cancer present on admission
- ICU stay present
- valid geography for exposure assignment

The current code uses a primary lung cancer definition based on diagnosis prefixes corresponding to malignant neoplasm of lung and bronchus. The cohort code was tightened to exclude non-primary codes such as respiratory failure diagnoses, personal history codes, and secondary pulmonary metastasis codes.

## Trajectory Phenotyping

Trajectories are built over the first 72 hours after ICU admission.

Each hospitalization is represented as an hourly sequence of:

- respiratory support tier
- acute respiratory failure subtype

Respiratory support is collapsed to:

- `ROOM AIR`
- `LOW_O2`
- `NIV`
- `IMV`
- `OTHER`

ARF subtype is collapsed to:

- `NO_ARF`
- `ARF_HYPOX`
- `ARF_HYPER`
- `ARF_MIXED`

Sequence clustering uses:

- optimal matching distance
- transition-rate substitution costs
- Ward hierarchical clustering

Sites choose the local number of clusters based on peak average silhouette width. Local clusters are then relabeled by severity so that the healthiest phenotype serves as the reference cluster in downstream models.

## Clinical and Exposure Variables

The workflow currently uses or derives:

- age
- sex
- race
- ethnicity
- admission year
- county-based PM2.5 exposure
- county-based NO2 exposure
- Charlson comorbidity score
- advanced cancer proxies
- SOFA total and domain scores
- respiratory support at ICU admission
- invasive mechanical ventilation exposure
- vasopressor exposure
- in-hospital death
- hospice discharge
- death or hospice discharge
- 72-hour landmark survival outcome

Advanced cancer proxies include:

- metastatic disease present on admission
- malignant pleural effusion present on admission
- cachexia present on admission

## Site-Level Models

Each site fits models locally using standardized exposures:

- multinomial logistic regression for trajectory phenotype
- logistic regression for death or hospice discharge
- generalized linear modeling for ICU length of stay
- ordinal regression for trajectory severity rank
- Cox proportional hazards modeling for post-72-hour death or hospice discharge

The pipeline also fits pollutant-by-cluster interaction models and generates exposure-response curves for:

- predicted mortality by cluster
- predicted ICU length of stay
- predicted trajectory severity

## Federated Outputs

The site return package is designed to be federated-friendly. Sites should return only aggregated outputs, coefficient tables, or fixed-grid predicted summaries.

Core exports include:

- `consort_flow_counts`
- `consort_exclusion_reasons`
- `table1_overall_continuous`
- `table1_overall_categorical`
- `table1_by_cluster_continuous`
- `table1_by_cluster_categorical`
- `cluster_silhouette_ra`
- `cluster_centroids_federated`
- `cluster_static_summary`
- `cluster_signature_rs_hourly`
- `cluster_signature_arf_hourly`
- `cluster_severity_rank`
- `cluster_sofa_summary`
- `cluster_charlson_summary`
- `model_outputs_standardized`
- `model_metadata`
- `mediation_style_decomposition_no2`

New epidemiology-oriented federated outputs include:

- `epi_overall_summary`
- `epi_device_t0_summary`
- `epi_admission_year_summary`
- `epi_cancer_severity_summary`
- `epi_start_support_to_cluster_counts`

These are sufficient for central pooling of:

- overall cohort epidemiology summaries
- initial respiratory support distributions
- year-level trends
- cancer burden summaries
- alluvial figures linking initial support to final trajectory phenotype

For the alluvial specifically, central pooling should combine each site’s `epi_start_support_to_cluster_counts` file, remap local clusters to pooled phenotypes using the central mapping table, and then sum counts across sites.

## Figures Generated Locally at Sites

The site script also generates local figures for quality control and presentation, including:

- sequence plots
- cluster SOFA profiles
- respiratory support and ARF signature figures
- predicted mortality by pollutant and cluster
- predicted ICU LOS response curves
- predicted severity response curves
- initial respiratory support distribution
- admission-year trend figure
- start-support to final-cluster alluvial

These figures are useful for local review, but central pooling should rely on the returned CSV tables rather than the site PNG files whenever possible.

## Privacy and Governance

This project is intended for federated execution. Please do not return patient-level files.

The current workflow avoids county-level spatial exports in the standard site return package because several sites expressed concern about sharing low-count county results. County maps and county-level cluster summaries are therefore not part of the recommended federated deliverables.

## Central Pooling

Central postprocessing scripts in [code/postprocessing](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/code/postprocessing) combine site returns into manuscript-facing pooled outputs, including:

- pooled phenotype harmonization
- manuscript Table 1 and Table 2
- pooled hourly trajectory figures
- pooled SOFA summaries
- pooled air pollution meta-analysis
- pooled interaction figures
- pooled epidemiology figures and tables

These scripts assume that each site has returned the standard PHI-safe output bundle generated by `02_federated_clusters.R`.

## Required CLIF Inputs

The analysis relies on the following CLIF tables or equivalent local sources:

- `patient`
- `hospitalization`
- `adt`
- `respiratory_support`
- `vitals`
- `labs`
- `medication_admin_continuous`
- `patient_assessments`
- `hospital_diagnosis`

Please refer to:

- [CLIF Data Dictionary](https://clif-icu.com/data-dictionary)
- [CLIF Tools](https://clif-icu.com/tools)
- [CLIF ETL Guide](https://clif-icu.com/etl-guide)

## Repo Structure

- [code](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/code): site-side analysis scripts
- [code/postprocessing](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/code/postprocessing): central pooling and manuscript figure scripts
- [config](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/config): configuration files
- [exposome](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/exposome): county-year exposure inputs
- [output](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/output): local and pooled outputs
- [sites](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/sites): returned site output folders used for central pooling

## Notes

- Michigan’s overall categorical Table 1 export has been missing in some returned site folders; pooled categorical summaries should be checked against site availability.
- Hopkins contributed an updated rerun with the corrected cluster reference logic.
- If sites rerun after code updates, central pooled outputs should be regenerated from the refreshed `sites/` folder.
