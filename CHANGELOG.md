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

- `mb_estimation`: with the initial value variant `preposterior_mode`, `fmincon` could stop at the initial point (flag -2, or a wrong result without an error), because this point lay on the boundary of the allowed region. It is now moved slightly towards the polytope centre (test: `tests/test_preposterior_mode.m`)
- `mb_estimation`: the variant `flat_objective` crashed, and the variant `average_sample` ignored the result of its first run. The first run of `average_sample` also did not receive the enzyme data (`model_balancing`)
- `mb_save_results`: `true_model.tsv` contained the balanced kinetic constants instead of the true ones. The Gibbs free energies in `balanced_states.tsv` had the wrong sign, and those in `true_states.tsv` had the wrong sign for negative fluxes
- `mb_generate_artificial_data`: true reaction affinities are now given in kJ/mol, like the balanced ones (before: in units of RT)
- `mb_graphics`: option `plot_true_vs_data` crashed
- `mb_log_posterior`: infeasible points gave the score NaN instead of -inf; an unknown enzyme score type now gives an error
- `mb_objective`: wrong value for the score `neg_log_posterior` when a gradient was requested

Changes that can affect results:

- The lower bound for Keq is now `quantities.Keq.min` from the prior table (before: 1 / `quantities.Keq.max`), both here and in `python-sbtab-to-json`
- Without a constraint file, the upper bound for metabolite concentrations is now `quantities.c.max` (before: 100 times this value), both here and in `python-sbtab-to-json`

Other changes:

- `sbtab_to_json.py`: kinetic data are limited to the bounds in the prior table, as in Matlab; new option `--km-max` (default 1000 mM) for data and bounds of KM
- Example files with `e_coli` or `ecoli` in their names renamed to `escherichia-coli` (`escherichia_coli` in Matlab function and script names)
- The 36 E. coli results in `examples/results/escherichia-coli_noor_2016` were recomputed with the corrected code
- Installation notes and demo descriptions moved from `INSTALLATION` and the inner `README.md` files to the main `README.md`, which also had a wrong folder in the Quick start
- `Contents.m` moved to the folder `model-balancing`, so that `help model-balancing` works again
- New function `mb_log_posterior_gradient`: analytic gradient, used by default (option `use_gradient`, now 1; test: `tests/test_gradient.m`). The old gradient code in `mb_log_posterior` was wrong and has been removed
- The variable `true` is now called `true_model` (also in the file `results.mat`)
- Parameterisation `KV_KM_KA_KI` removed (it crashed), together with the option `parameterisation`
- Unused code commented out: `mb_read_prepared_model_and_data`, `mb_variables_to_log_posterior`, `mb_qall_to_q`, the tighter bounds in `mb_estimation`, and the computation of `gradient_V`
- Folder `resource-data` renamed to `resources`, folder `model-balancing/model-balancing` to `model-balancing/model-balancing-functions`
- Corrections in function headers and comments

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
