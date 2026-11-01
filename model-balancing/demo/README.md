Matlab demo scripts for model balancing
=======================================

* `demo_mb_artificial_data.m`: reads a model from an SBML or SBtab file, generates artificial data
  (kinetic constants and metabolic state data), and runs model balancing to try and recover the true model and states.

* `demo_mb_experimental_data.m`: reads a model and data (kinetic data and metabolic state data) from SBtab files,
  runs model balancing, and saves the results. This script is the best starting point for balancing your own model
  (see the section "Your model" in the script). It also shows how to export a model balancing problem to a JSON file
  for the python implementation.

* `model_balancing_figure_2.m`: generates Figure 2 of the article (Liebermeister and Noor, 2021), showing the enzyme
  posterior score for different values of the stringency parameter alpha.

By default, the demo scripts use the model `branch_point_model`. The other example models (from `resource-data/models`)
can be used by uncommenting the respective lines in the scripts. Output files are written to matlab's `tempdir`
(variable `result_dir`).
