Model balancing results for test models
=======================================

This directory contains model balancing results (computed with the matlab implementation)
for the example models in `resource-data/models`. There is one subdirectory per model and data type:

* `branch_point_model_artificial`, `double_branch_model_artificial`, `three_chain_model_artificial`, `e_coli_artificial`:
  toy models and E. coli model with artificial data
* `e_coli_noor_2016`: E. coli model with experimental data

## Names of the result folders

Each result folder corresponds to one model balancing run (estimation scenario and settings).

Models with artificial data: `[STATE]_state_[KINETIC]_alpha_[ALPHA]_beta_[BETA]`

* `[STATE]`: `noisefree` or `noisy` artificial state data (metabolite concentrations, enzyme concentrations, fluxes)
* `[KINETIC]`: kinetic data used
  * `no_kinetic`: no kinetic data
  * `noisefree_Keq` or `noisy_Keq`: equilibrium constants only (without or with noise)
  * `noisefree_kinetic` or `noisy_kinetic`: all kinetic constants (without or with noise)
* `[ALPHA]`: stringency parameter alpha for the enzyme fit (option `cmb_options.enzyme_score_alpha`)
* `[BETA]`: stringency parameter beta for c/KM ratios (option `cmb_options.beta_ln_c_over_km`)
* Suffix `_with_crosscovariances` (where present): run with cross-covariances between types of kinetic constants (option `cmb_options.use_crosscovariances`)

E. coli model with experimental data: `[KINETIC]_data_[KINETIC_DATA_SET]_alpha[ALPHA]`

* `[KINETIC]`: `no_kinetic`, `with_Keq`, or `with_kinetic`
* `[KINETIC_DATA_SET]`: name of the kinetic data set (`original` or `balanced`; option `cmb_options.kinetic_data_set`);
  see `resource-data/models/e_coli_noor_2016/e_coli_ccm_kinetic_data`
* `[ALPHA]`: stringency parameter alpha (for example, `alpha0_01` means alpha = 0.01)

## Files in each result folder

* `[MODEL_NAME]_[SCENARIO].json`: model balancing problem in JSON format (input for the python implementation;
  see `examples/json/README.md`)
* `data/report.txt`: summary (numbers of data points and fitted variables, options used, estimation errors, calculation time)
* `data/options.tsv`: options used (SBtab)
* `data/balanced_model.tsv`, `data/balanced_states.tsv`: balanced kinetic model and metabolic states (SBtab)
* `data/data_kinetic.tsv`, `data/data_states.tsv`: kinetic data and state data used (SBtab)
* `data/true_model.tsv`, `data/true_states.tsv`: true kinetic model and states (SBtab; only for artificial data)
* `data/results.mat`: all results as matlab variables (see `help model_balancing`)
* `ps-files/`: scatter plots (PDF format; written by `cmb_graphics`)
  * `data_[QUANTITY].pdf`: balanced values compared to data values
  * `true_[QUANTITY].pdf`: balanced values compared to true values (only for artificial data)
  * `diagnostic_errors.pdf`, `diagnostic_flux_signs.pdf`: diagnostic plots

  Here `[QUANTITY]` is one of `x` (metabolite concentrations), `e` (enzyme concentrations), `KV`, `KM`, `Keq`,
  `Kcatf`, `Kcatr` (kinetic constants), or `Kappmaxf`, `Kappmaxr` (maximal apparent catalytic rates, compared to kcat values).

Note: `cmb_graphics` writes the graphics as EPS files; in this directory, they were converted to PDF to save space.
