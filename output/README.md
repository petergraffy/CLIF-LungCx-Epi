# New run outputs

New runs use `output/run_<SITE>_<RUN_ID>/`.

- `local/` contains patient-level caches. Keep it at the site; never upload it.
- `share/` contains aggregate derivation returns and run manifests for disclosure review.
- `share/models_<CONSENSUS_ID>/` contains second-pass model returns against the frozen mapping.

Historical output directories are not new-run inputs. Do not rescale counts, fill suppressed cells with zero, or combine old and new returns. The coordinator supplies explicit return directories.
