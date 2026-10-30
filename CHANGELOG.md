# Changelog

Version numbers have the form MAJOR.MINOR.PATCH. MAJOR changes when existing calls may stop working.

## 1.0.0 (not released yet)

First numbered version.

- Version function `model_balancing_version` and dependency check `model_balancing_dependencies`
- Improved documentation; bug fixes in demo scripts
- SBtab 1.1 options tables
- `state_data_to_data` moved here from the Metabolic Network Toolbox
- `model_balancing` checks the dependencies automatically, once per Matlab session
- `cmb_estimation`: fixed misspelled variable names that crashed model balancing when the initial point had to be computed
- `model_balancing_dependencies` also checks Enzyme Cost Minimization
- Minimal test suite (folder `tests`)

The state before the cleanup of October 2026 is available under the git tag `v0.9-legacy`.
