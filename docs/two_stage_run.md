# Two stage deployment contract

This deployment packages the adopted three-site rehearsal method. It does not
use the legacy raw-profile one-to-one consensus commands. The run instructions are in `code/README.md`. The initial site job
is scripts/site.R stage1. The coordinator alone releases Stage 2 parameters.

Stage 1 produces local OM/Ward discoveries for primary, h24, h48, onset,
support_only, low_flow, carry0, carry12, raw_support and descriptive absorbing
sequences. Fingerprints retain time order through twelve relative-time bins.
Local k and stability are measured locally; shared k and templates do not exist
until the coordinator reviews and freezes the aggregate returns.

The coordinator inspect command checks sites and contracts without learning or
assigning shared phenotypes. The freeze command uses only discovery fingerprints,
cluster eligibility/stability and respiratory features. Provisional descriptors,
outcomes, SOFA, APACHE, code status and vasopressors are not matching features.
Coordinator-selected k must be an eligible candidate with cross-site support.
Null k uses the supported silhouette maximum. Full and held-out variants are
independently frozen. The current release selects four primary phenotypes; sensitivity k values are
recorded separately. New releases still require a supported cross-site solution.

Two scientific stages require three site return points:

1. Stage 1 discovery and diagnostics -> coordinator parameters/templates.
2. Stage 2 fit associations and predictive training aggregates -> coordinator
   frozen other-site predictors and common relative-event-time grid.
3. Stage 2 evaluation and risk counts -> final coordinator pooling.

The latter exchange is necessary to preserve held-out validation and exact
pooled risk sets. It is not a second discovery pass or an instruction to change
local clusters. Sites do not fit the final held-out predictor on their own data.
The final risk/event counts support both death/hospice and live non-hospice
Aalen-Johansen incidence, with the other endpoint competing. No identifiers,
individual sequences, individual predictions or calendar timestamps are exchanged.

All stages have explicit CSV/JSON allowlists, hashed inventories, source/lockfile
contracts and immutable step markers. Local stdout/stderr logs stay in local/steps.
A successful site return is aggregate capture, not a claim that every requested
clinical variable or model is available. Nonestimable effects and absent source
fields remain explicit. The source commit and local run caches must be preserved.

Stage 1 may take substantial time: exact optimal matching is quadratic in local
cohort size, and the primary/sensitivity discovery passes include 100 patient
bootstraps. Run on an institutional machine with enough memory for distance
matrices and keep it available for Stage 2. There is no automatic subsampling or
silent relaxation of bootstrap/cluster settings. Contact the coordinator if a
site cannot complete a prespecified analysis; do not submit a partial directory
as a complete return.

Stage 2 is the first workflow in the code README. A remaining late site follows
the retained Stage 1 instructions and waits for coordinator admission; it does
not update the frozen mapping. Stage 2 uses a separate, release-specific local
workspace and verifies the original manifest and local discovery exports before
copying caches. It preserves cached instants and exposure membership without
reinterpreting source timestamps or modifying Stage 1 source hashes.

The coordinator may explicitly accept legacy Stage 1 snapshots. All participating
sites are expected to contribute to all primary, sensitivity and outcome analyses.
Required discovery inputs must pass their original checksums; pending files are
repair items, not reasons to predesignate site exclusions. A missing ancillary
descriptor does not invalidate a valid primary fingerprint. Leave-one-site-out
validation excludes each evaluation site from training in turn; it does not
exclude sites from the study. Current Stage 2 code and environment remain
strictly matched to the released bundle.
