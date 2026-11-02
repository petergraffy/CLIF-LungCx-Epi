Federated discovery and outcome workers

Sites and coordinators use code/README.md and scripts/site.R / coordinator.R.
These are reusable worker modules, not independent site entry points.

contract.R defines the respiratory fingerprint, common-cost sequence method,
nearest-template assignment and supplied participating-site list. Site workers
read only their local clinical caches. Central discovery uses aggregate
fingerprints and cluster stability; outcome and vasopressor data do not enter
matching, labels or reference selection. Templates are frozen before modeling.

Stage 1 returns local discovery for coordinator assessment. Stage 2 applies the
coordinator-issued full/held-out templates, estimates site associations and
trains aggregate predictive coefficients. Subsequent evaluation uses frozen
other-site coefficients; no held-out-site refitting occurs.

Local records, sequences, fits and individual predictions stay in the ignored
output directory. Only designated aggregate return folders are transferred.
Retired three-site runners/results are not shipped in the shared tree.

Synthetic worker checks: Rscript tests/federated_rehearsal.R
Portable handoff checks: Rscript --vanilla tests/deployment.R
