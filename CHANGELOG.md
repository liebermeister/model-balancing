# Changelog

Version numbers have the form MAJOR.MINOR.PATCH. MAJOR changes when existing calls may stop working.

## 1.0.0 (2026-10-31)

First numbered version.

- HTML documentation removed from the repository (it can be regenerated with m2html)
- Version function `model_balancing_version` and dependency check `model_balancing_dependencies`
- Improved documentation; bug fixes in demo scripts
- SBtab 1.1 options tables
- `state_data_to_data` moved here from the Metabolic Network Toolbox
- `model_balancing` checks the dependencies automatically, once per Matlab session
- `cmb_estimation`: fixed misspelled variable names that crashed model balancing when the initial point had to be computed
- `model_balancing_dependencies` also checks Enzyme Cost Minimization
- Minimal test suite (folder `tests`)
- `cmb_load_and_convert_model_and_data` passes the prior table of this package (`cmb_prior_file`) to `sbtab_load_network_model_and_data`

The state before the cleanup of October 2026 is available under the git tag `v0.9-legacy`.
