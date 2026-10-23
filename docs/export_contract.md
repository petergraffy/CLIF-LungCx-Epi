# Federated export contract

All new CSVs carry `schema_version`, `site`, and `run_id`. Source manifests record the scientific parameters, code MD5 values, R version, and package versions. The central derivation pass rejects incompatible source versions/parameters and duplicate site packages. Model packages refer to a single hashed `consensus_id` and a shared reference.

## Derivation pass

| Family | Content |
|---|---|
| `analysis_plan` | Primary/sensitivity role, trajectory origin, window, clustering inputs, outcome-modeling status |
| `arf_detection_audit` | Aggregate observed ARF detection and initial-evidence categories |
| `cohort_exclusions` | Mutually exclusive cohort exclusion counts |
| `outcome_timestamp_audit` | Disposition mapping, event-time source, validity, counts |
| `source_availability` | Availability of optional code-status and oncology sources |
| `h*_landmark_flow` | Early death, hospice, live discharge, unknown outcome, eligible counts |
| `h*_sequence_flow` | Included counts; missing support never excludes a hospitalization |
| `h*_signature` | Local cluster, hour, joint state, state count, cluster N, probability, suppression status |
| `h*_silhouette` | Candidate k, mean silhouette, selected k |
| `h*_stability` | Completed patient bootstrap replicates, cluster Jaccard and adjusted Rand summaries |
| `h*_continuous` | Observed/missing N, sums, sum of squares, mean by variable and cluster |
| `h*_categorical` | Counts, explicit denominator, percentage by variable and cluster |
| `h*_missingness` | Support/physiology coverage and missing severity counts |
| `derivation_status` | Explicit success/failure by analysis |

`local/arf_detection.rds`, `local/prepared.rds`, and `local/derivation.rds` are patient-level and **must remain local**. The latter stores sequences and assignments so the second pass uses the same derivation. Deleting local caches requires a new derivation; do not attempt to reconstruct patient assignments centrally.

## Consensus pass

The coordinator returns `consensus_mapping.csv`, keyed by site, run, analysis, and local cluster. It includes phenotype, shared reference, matching distance, matching status, and consensus version. Preserve unmatched clusters. The package also retains source signatures, prototypes, LOSO results, and consensus status. A consensus with outcomes as states is explicitly labeled descriptive-only.

The consensus package also contains pooled phenotype continuous and categorical tables. Continuous means and standard deviations use summed counts, sums, and sums of squares. Categorical denominators are counted once per site and variable. If a contributing site partition is suppressed or a required site variable is missing, pooled totals remain missing and carry an explicit incomplete status. Site-level cohort, annual, silhouette, stability, and code-status outcome tables are retained for audit.

## Model pass

Each `share/models_<ID>/` folder contains the frozen mapping, phenotype main effects, covariance matrices, model formulas and actual N/events, covariate missingness, mapping coverage, fixed-grid competing-risk curves, and local cross-validation audit/performance/calibration/decision tables. Sensitivities use separate identifiers. Model failures produce metadata rather than fabricated zero effects. The coordinator retains site validation metrics separately; it does not average calibration slopes into a global validation claim.

## Suppression and reconciliation

The default minimum nonzero cell count is ten; use the locally approved threshold. If a component count in a partition is small, all numeric result values in that partition are masked, including complementary cells. For hourly state profiles, the entire distribution for that cluster-hour is masked. Small clusters have their cluster size and full profile masked. Zero counts remain genuine zero only in unsuppressed partitions. CSV blanks never mean zero.

Exact denominators, model sizes, percentages, and linked output families can create disclosure risks across tables. Automated rules are conservative safeguards, not a proof against reconstruction across every output. Local disclosure review of the complete package remains required. If suppression leaves insufficient common profile coverage for consensus, the workflow stops or leaves profiles unmatched; the coordinator must agree on a permissible aggregate representation with sites, not infer hidden counts.

Central code never rescales site counts, imputes suppressed cells, interpolates missing site-years, or averages site medians. Counts and percentages should be reconciled against unsuppressed denominators locally. Historical return formats are rejected rather than repaired.
