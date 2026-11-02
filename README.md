# CLIF lung cancer respiratory trajectories

Federated temporal phenotyping of respiratory support and acute respiratory
failure evidence in ICU patients with lung cancer, using CLIF data.

The primary analysis characterizes the first 72 hours after ICU arrival among
hospitalizations alive and in hospital at the landmark, then evaluates subsequent
hospital death, hospice discharge and live non-hospice discharge. Local sequence
discovery and aggregate fingerprints support centrally learned, frozen shared
phenotypes. Patient-level data remain at each institution.

**Sites and coordinators: see [the code README](code/README.md) for environment
setup, configuration, Stage 1 and Stage 2 commands, and required returns.**

Supporting documentation:

- [Analysis protocol](docs/analysis_protocol.md)
- [Reviewer requirements](docs/site_return_requirements.csv)
- [Site return reference](docs/site_return_workflow.md)
- [Reproducible environment](docs/environment.md)
- [Archived submitted implementation](archive/2026-10-21/README_ARCHIVE.md)

The current workflow is a new analysis for the manuscript revision. Historical
returns and figures are preserved separately and are not automatically used as
inputs to a new federation. Modified SOFA and adapted APACHE comparisons are
explicitly qualified adaptations rather than validated standard scores.
