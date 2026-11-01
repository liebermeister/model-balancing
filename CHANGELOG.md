# Changelog

Version numbers have the form MAJOR.MINOR.PATCH. MAJOR changes when existing calls may stop working.

## Unreleased

All functions and files starting with `cmb_` now start with `mb_` (for example `mb_estimation`, `mb_default_options` and `mb_prior.tsv`), and the struct `cmb_options` is now `mb_options`. The folders `cmb-variables`, `cmb-file-input-output` and `cmb-cvxpy` and the demo scripts `demo_cmb_...` now use `mb` as well.

- New folder `python-sbtab-to-json`: stand-alone conversion of SBtab input files to JSON files for the python implementation
- New folder `matlab-python-comparison`: test comparing the matlab and python implementations

Bug fixes that change results for real data:

- `mb_prepare_posterior`: metabolite and enzyme levels without data now keep their prior mean and standard deviation. Before, metabolites without data kept a tightened standard deviation (when no error bars were given) and enzymes without data got a standard deviation of about 1e15
- `state_data_to_data`: the logarithm of the kinetic data errors was taken twice (`log(std_ln)` instead of `std_ln`), so kinetic data had a wrong weight and JSON files contained negative standard deviations. Example JSON files and results made before this fix are affected

Other bug fixes:

- Networks with closed cycles (rank-deficient stoichiometric matrix) can now be balanced: `mb_define_parameterisation`, `mb_estimation`, `mb_qall_to_kinetics`, `mb_make_bounds`, `mb_prepare_posterior`
- `mb_default_options`: upper bound of KA was taken from the row of the prior table for A
- New function `mb_save_json_file`, used in `demo_mb_experimental_data` (the function `save_json_file` is not part of this package)
- Tests for these fixes (`tests/test_data_handling.m`)

Known issue: with precise kinetic data, the initial value variant `preposterior_mode` can make `fmincon` stop at the initial point, with flag -2 or without an error. This point lies on the boundary of the allowed region, where the interior-point algorithm cannot move. The default `polytope center` works.

## 1.0.0 (2026-10-31)

First numbered version.

- HTML documentation removed from the repository (it can be regenerated with m2html)
- Version function `model_balancing_version` and dependency check `model_balancing_dependencies`
- Improved documentation; bug fixes in demo scripts
- SBtab 1.1 options tables
- `state_data_to_data` moved here from the Metabolic Network Toolbox
- `model_balancing` checks the dependencies automatically, once per Matlab session
- `mb_estimation`: fixed misspelled variable names that crashed model balancing when the initial point had to be computed
- `model_balancing_dependencies` also checks Enzyme Cost Minimization
- Minimal test suite (folder `tests`)
- `mb_load_and_convert_model_and_data` passes the prior table of this package (`mb_prior_file`) to `sbtab_load_network_model_and_data`

The state before the cleanup of October 2026 is available under the git tag `v0.9-legacy`.
