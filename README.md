# CLIF Lung Cancer ICU Trajectory Phenotyping

CLIF version: 2.1

## Overview

This repository contains a federated analysis pipeline for trajectory-based phenotyping of respiratory failure among patients with lung cancer admitted to the ICU across CLIF sites.

The main scientific focus of this repo is the clinical and trajectory component of the project: identifying reproducible early respiratory phenotypes, describing their organ dysfunction and cancer severity profiles, and quantifying their associations with ICU length of stay and mortality-related outcomes. Sites run the same code locally and return only PHI-safe aggregated outputs for central pooling.

## Study Focus

The workflow in this repository supports a manuscript framed around:

1. defining a multicenter ICU cohort of adult hospitalizations with lung cancer present on admission
2. reconstructing hourly respiratory support and acute respiratory failure states during the first 72 hours after ICU admission
3. using state sequence clustering to identify consensus respiratory trajectory phenotypes
4. characterizing those phenotypes using SOFA, Charlson, vasopressor exposure, cancer severity proxies, initial support, and clinical outcomes
5. evaluating associations between phenotype membership and ICU length of stay, in-hospital death or hospice discharge, and post-72-hour survival

Some legacy code and outputs related to long-term air pollution exposures remain in the repository because they were part of the broader project workflow, but they are not the primary emphasis of this README.

## High-Level Workflow

Sites should run the scripts in [code](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/code) in this order:

1. [00_renv_restore.R](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/code/00_renv_restore.R)
2. [01_lungcx_cohort.R](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/code/01_lungcx_cohort.R)
3. [02_federated_clusters.R](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/code/02_federated_clusters.R)

### Step 1: Environment restore

`00_renv_restore.R` restores the project R environment and installs required packages.

### Step 2: Cohort construction

`01_lungcx_cohort.R` builds the analytic cohort and creates the core objects used downstream:

- `cohort_lung`
- `analysis_ready`
- `flow_lung`
- `exclusion_lung`

This step:

- identifies adult hospitalizations in the study window
- restricts to hospitalizations with lung cancer present on admission
- requires an ICU stay
- assigns the ICU index time
- derives respiratory support summary features
- constructs cohort-level outcomes and baseline variables

### Step 3: Trajectory construction and federated modeling

`02_federated_clusters.R` performs local trajectory construction, sequence clustering, cluster characterization, and federated export generation.

This step:

- reconstructs hourly respiratory support and ARF states over the first 72 hours
- computes optimal matching distances with `TraMineR`
- selects the local number of clusters using peak average silhouette width
- reorders local clusters so the healthiest phenotype is the reference cluster
- summarizes clusters using SOFA, Charlson, vasopressors, and cancer severity proxies
- fits site-level outcome models
- exports aggregated tables, model coefficients, trajectory summaries, and site figures

## Cohort Definition

The intended analytic cohort is:

- adult hospitalizations
- lung cancer present on admission
- ICU stay present

The current code defines lung cancer present on admission as at least one diagnosis flagged POA that matches a primary lung cancer code. In the current codebook this corresponds to ICD-10-CM `C34.xx` malignant neoplasm of bronchus/lung codes: `C34.00`, `C34.01`, `C34.02`, `C34.10`, `C34.11`, `C34.12`, `C34.20`, `C34.21`, `C34.30`, `C34.31`, `C34.32`, `C34.80`, `C34.90`, `C34.91`, and `C34.92`. The cohort function also supports ICD-9-CM primary lung cancer codes with prefix `162` if present in a site's source data.

The cohort code was tightened to exclude non-primary codes such as personal history codes and secondary pulmonary metastasis codes.

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

Sites choose the local number of clusters based on peak average silhouette width. Local clusters are then relabeled by severity so that the healthiest local phenotype serves as the reference cluster in downstream models. Central pooling harmonizes those local clusters into pooled consensus phenotypes.

## Clinical Characterization

The trajectory phenotypes are characterized using variables such as:

- age
- sex
- race
- ethnicity
- admission year
- Charlson comorbidity score
- advanced cancer proxies
- SOFA total and domain scores
- respiratory support at ICU admission
- invasive mechanical ventilation exposure
- vasopressor exposure
- ICU length of stay
- in-hospital death
- hospice discharge
- death or hospice discharge
- post-72-hour survival

Advanced cancer proxies include:

- metastatic disease present on admission
- malignant pleural effusion present on admission
- cachexia present on admission

## Site-Level Models

Each site fits models locally for the manuscript-facing clinical analyses, including:

- multinomial logistic regression for trajectory phenotype membership
- logistic regression for death or hospice discharge
- generalized linear modeling for ICU length of stay
- ordinal regression for trajectory severity rank
- Cox proportional hazards modeling for post-72-hour death or hospice discharge

These models are designed to support pooled interpretation of how respiratory trajectory phenotype relates to severity, resource use, and survival.

## Federated Outputs

The site return package is designed to be federated-friendly. Sites should return only aggregated outputs, coefficient tables, or fixed-grid summaries.

Core trajectory and clinical exports include:

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
- `model_outputs_standardized`
- `model_metadata`

Epidemiology-oriented exports used for pooled figures and manuscript tables include:

- `epi_overall_summary`
- `epi_device_t0_summary`
- `epi_admission_year_summary`
- `epi_cancer_severity_summary`
- `epi_start_support_to_cluster_counts`

These are sufficient for central pooling of:

- overall cohort epidemiology summaries
- phenotype-specific clinical profiles
- initial respiratory support distributions
- start-support to phenotype alluvial figures
- year-level trend summaries
- cancer burden summaries

## Figures Generated Locally at Sites

The site script also generates local figures for quality control and presentation, including:

- sequence plots
- respiratory support signature figures
- ARF signature figures
- cluster SOFA profiles
- cluster severity summaries
- initial respiratory support distribution
- start-support to final-cluster alluvial
- admission-year trend figures

These figures are useful for local review, but central pooling should rely on the returned CSV tables rather than the site PNG files whenever possible.

## Privacy and Governance

This project is intended for federated execution. Please do not return patient-level files.

The workflow is structured so that each site runs the same code locally and shares only PHI-safe aggregated outputs for central pooling.

## Central Pooling

Central postprocessing scripts in [code/postprocessing](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/code/postprocessing) combine site returns into manuscript-facing pooled outputs, including:

- pooled phenotype harmonization
- manuscript Table 1 and Table 2
- pooled hourly trajectory figures
- pooled SOFA summaries
- pooled epidemiology figures and tables
- pooled start-support to phenotype alluvial figures

Some central scripts also retain broader project outputs from prior analyses, including air-pollution-related files, but the main purpose of the pooled workflow described here is the trajectory and clinical paper.

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
- [output](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/output): local and pooled outputs
- [sites](/Users/saborpete/Desktop/Peter/Postdoc/CLIF-LungCx-Epi/sites): returned site output folders used for central pooling

## Notes

- This README is intentionally centered on the trajectory and clinical manuscript framing.
- Air-pollution-related code remains in the repo as part of the broader project history, but it is not the primary focus of this documentation.
