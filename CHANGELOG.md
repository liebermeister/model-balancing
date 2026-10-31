# Changelog

Version numbers have the form MAJOR.MINOR.PATCH. MAJOR changes when existing calls may stop working.

## Unreleased

Bug fixes that change results for real data:

- `cmb_prepare_posterior`: metabolite and enzyme levels without data now keep their prior mean and standard deviation. Before, metabolites without data kept a tightened standard deviation (when no error bars were given) and enzymes without data got a standard deviation of about 1e15
- `state_data_to_data`: the logarithm of the kinetic data errors was taken twice (`log(std_ln)` instead of `std_ln`), so kinetic data had a wrong weight and JSON files contained negative standard deviations. Example JSON files and results made before this fix are affected

Other bug fixes:

- Networks with closed cycles (rank-deficient stoichiometric matrix) can now be balanced: `cmb_define_parameterisation`, `cmb_estimation`, `cmb_qall_to_kinetics`, `cmb_make_bounds`, `cmb_prepare_posterior`
- `cmb_default_options`: upper bound of KA was taken from the row of the prior table for A
- New function `cmb_save_json_file`, used in `demo_cmb_experimental_data` (the function `save_json_file` is not part of this package)
- Tests for these fixes (`tests/test_data_handling.m`)

Known issue: with precise kinetic data, the initial value variant `preposterior_mode` can make `fmincon` fail (flag -2); the default `polytope center` works.

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
