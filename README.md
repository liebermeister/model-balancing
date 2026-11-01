Model balancing
===============

[Model balancing](https://www.metabolic-economics.de/model-balancing/index.html) is a computational method to
determine plausible kinetic constants and metabolic states in kinetic metabolic models. It integrates flux,
metabolite, protein, and kinetic constant data, using prior distributions for all these variables, and
computes the joint posterior mode.
Model balancing can be run in matlab or python. Data tables can be provided in [SBtab](https://www.sbtab.net)
format, models can be provided in [SBML](https://sbml.org) or [SBtab](https://www.sbtab.net) format.

This repository contains the matlab implementation, example models and data, and example results.

## Bug fixes of October 2026: effect on earlier results

Two bugs in versions up to 1.0.0 affected results for real data. Both are fixed in the current code (see `CHANGELOG.md`).

* **Missing data** (`mb_prepare_posterior`): metabolites without data kept a too narrow distribution if the data file had no error bars, which pulled them towards the prior mean. Enzymes without data got a standard deviation of about 1e15 instead of the prior value.
* **Kinetic data errors** (`state_data_to_data`): the logarithm of the error was taken twice. For example, a geometric standard deviation of 1.5 acted like a standard deviation of 0.90 on the log scale, instead of 0.405, so kinetic data had about one fifth of their intended weight. JSON files written by the old code contain negative standard deviations.

The results for real data in the article (Liebermeister & Noor, 2021) were computed with the old code; the article has not been changed. Results for artificial data are not affected, since they have no missing values and set the kinetic data errors directly. The results in the folder `examples/results` and the files in `examples/json` were also made with the old code, so please regenerate them before using them.

The table shows the effect on the *E. coli* model with real data (balanced kinetic data, alpha = 0.5), as the median of |ln(fit/data)|, old code -> fixed code.

| Scenario          | Keq          | kcat+        | kcat-        | KM           | Metabolites  |
|-------------------|--------------|--------------|--------------|--------------|--------------|
| With kinetic data | 0.46 -> 0.22 | 0.54 -> 0.32 | 0.24 -> 0.09 | 0.43 -> 0.15 | 0.10 -> 0.22 |
| Keq data only     | 0.18 -> 0.03 | 1.78 -> 2.14 | 2.48 -> 2.31 | 1.67 -> 1.68 | 0.05 -> 0.06 |
| No kinetic data   | 3.28 -> 3.20 | 1.57 -> 1.57 | 1.30 -> 1.33 | 1.75 -> 1.75 | 0.01 -> 0.01 |

With kinetic data, the model now fits the kinetic data better, and estimates moved by a median factor of about 1.3 (Keq: 1.65). Metabolites without data moved most, by a median factor of about 2, and up to about 40 (up to about 100 with Keq data only). Without kinetic data, the estimates hardly changed.

## Quick start (matlab)

1. Clone this repository and the four required toolboxes (see [Installation](#installation)):
   ```
   git clone https://github.com/liebermeister/model-balancing.git
   git clone https://github.com/liebermeister/matlab-utils.git
   git clone https://github.com/liebermeister/metabolic-network-toolbox.git
   git clone https://github.com/liebermeister/sbtab-matlab.git
   git clone https://github.com/liebermeister/enzyme-cost-minimization.git
   ```
2. In matlab, add all of them (including subdirectories) to your matlab path, for example
   ```matlab
   addpath(genpath('model-balancing/model-balancing'));
   addpath(genpath('matlab-utils'));
   addpath(genpath('metabolic-network-toolbox'));
   addpath(genpath('sbtab-matlab'));
   addpath(genpath('enzyme-cost-minimization'));
   ```
3. Run one of the demo scripts:
   ```matlab
   demo_mb_artificial_data      % generate artificial data for a small model and balance it
   demo_mb_experimental_data    % balance a model with data read from SBtab files
   ```
   By default, output files are written to matlab's `tempdir`.

## Installation
The following installation instructions are only for the matlab version.
For python, please refer to our [Read The Docs](https://model-balancing.readthedocs.io/en/latest/index.html) page.

- Clone the following [GitHub](https://github.com/liebermeister) repositories
    - [`matlab-utils`](https://github.com/liebermeister/matlab-utils) - utility functions
    - [`metabolic-network-toolbox`](https://github.com/liebermeister/metabolic-network-toolbox) - metabolic network toolbox
    - [`sbtab-matlab`](https://github.com/liebermeister/sbtab-matlab) - SBtab toolbox
    - [`enzyme-cost-minimization`](https://github.com/liebermeister/enzyme-cost-minimization) - enzyme cost minimization toolbox
- [SBML toolbox](https://sbml.org/software/sbmltoolbox/) (optional - needed only if SBML files are used;
  the version needs to match your matlab version and requires a matching version of libSBML)
- Make sure that all these directories and their subdirectories are included in your matlab path.

The code was tested with Matlab R2025b on Linux. The optimisation uses `fmincon`
(Optimization Toolbox).

## Contents of this repository

| Folder          | Contents |
|-----------------|----------|
| `model-balancing` | Matlab code: main function `model_balancing.m`, functions (`model-balancing/model-balancing-functions`), demo scripts (`model-balancing/demo`) |
| `resources` | Example models (`resources/models`) and the prior table (`resources/data/data-prior`) |
| `examples`      | Input files for the python implementation (`examples/json`), model balancing results (`examples/results`), and result graphics (`examples/graphics`) |
| `python-sbtab-to-json` | Stand-alone python script that converts SBtab input files into the JSON files for the python implementation |
| `matlab-python-comparison` | Test that compares the matlab implementation, the converter, and the python implementation |

### Matlab
To run model balancing in matlab, please adjust the demo scripts in `model-balancing/demo` to your model and data.
For help, please type `help model-balancing` in matlab, or see the documentation at
[www.metabolic-economics.de/model-balancing/matlab-doc/](https://www.metabolic-economics.de/model-balancing/matlab-doc/).

Demo scripts (folder `model-balancing/demo`):

* `demo_mb_artificial_data.m`: reads a model from an SBML or SBtab file, generates artificial data
  (kinetic constants and metabolic state data), and runs model balancing to try and recover the true model and states.
* `demo_mb_experimental_data.m`: reads a model and data (kinetic data and metabolic state data) from SBtab files,
  runs model balancing, and saves the results. This script is the best starting point for balancing your own model
  (see the section "Your model" in the script). It also shows how to export a model balancing problem to a JSON file
  for the python implementation.
* `demo_model_balancing_figure_2.m`: generates Figure 2 of the article (Liebermeister and Noor, 2021), showing the enzyme
  posterior score for different values of the stringency parameter alpha.

By default, the demo scripts use the model `branch_point_model`. The other example models (from `resources/models`)
can be used by uncommenting the respective lines in the scripts. Output files are written to matlab's `tempdir`
(variable `result_dir`). Matlab scripts for generating the example models can be found in the subfolders `matlab`
of the model folders.

Implementation notes:

* Optimisation is performed by the fmincon function with the interior-point method.
* For numerical reasons, the matlab implementation does not tolerate equal lower and upper bounds for metabolites. If c_min == c_max, 10^-5 is added automatically to ln c_max. This also speeds up the calculation.
* With the option `initial_values_variant = 'average_sample'`, the algorithm starts by running model balancing on an average model state (with metabolic state data given by the geometric mean over the metabolic states). The resulting kinetic constants are used as initial values for the full calculation with several metabolic states.

### Python
For using model balancing in python, please refer to our [Read The Docs](https://model-balancing.readthedocs.io/en/latest/index.html) page;
the code itself can be found on [our GitLab repository](https://gitlab.com/elad.noor/model-balancing).
The script `python-sbtab-to-json/sbtab_to_json.py` makes the JSON input files from the same SBtab files as used in matlab, so matlab is not needed for this step.

### Example models
Example model and data files can be found in the folder `resources`.
The folder `resources/models` contains files for four different example models.
Each model comes as an SBML (`.xml`) and an SBtab (`.tsv`) file.
The folders `examples/results/[MODEL_NAME]` contain balancing results for a large number of model balancing problems
with these models: artificial data for all 4 models, experimental data for the E. coli model (folder `escherichia-coli_noor_2016`).

### How to balance your own model
To balance your own model, please have a look at the matlab demo script `model-balancing/demo/demo_mb_experimental_data.m`.
You need to prepare three input files in SBtab format, describing

1. the network structure and metabolite constraints (tables Reaction, Compound, Position (optional), Parameter (optional))
2. kinetic data (table ParameterData)
3. state data (tables FluxData, MetaboliteConcentrationData, EnzymeConcentrationData)

For an example, please have a look at the files `artificial_network_true.tsv`, `artificial_kinetic_data.tsv`,
and `artificial_state_data.tsv` in the folder `resources/models/branch_point_model/data/`.

To balance your model in matlab, you can then use the code from `demo_mb_experimental_data.m`.
Please note that there are many more settings that you can change. For an overview, type `help mb_default_options` in matlab.

To balance your model in python, you first need to convert your model and data into a JSON file.
This can be done in matlab and is shown in `demo_mb_experimental_data.m`. For an overview of the JSON format,
type `help cvxpy_problem_data_structure` in matlab, or see `examples/json/README.md`.

## Tests

From the repository folder, run

```matlab
results = runtests('tests')
```

## Versions

To see the installed version, type `model_balancing_version`. Older versions can be downloaded from the [releases page](https://github.com/liebermeister/model-balancing/releases) of this repository, and changes between versions are listed in `CHANGELOG.md`.

## License
This package is released under the [MIT License](LICENSE).

## Contact
Please contact [Wolfram Liebermeister](mailto:wolfram.liebermeister@gmail.com)
and [Elad Noor](mailto:elad.noor@weizmann.ac.il) with any questions or comments.

## How to cite
If you use model balancing, please cite:

Liebermeister W. and Noor E. (2021), *Model balancing: a search for in-vivo kinetic constants and consistent
metabolic states*, Metabolites 11(11):749,
[doi:10.3390/metabo11110749](https://doi.org/10.3390/metabo11110749)

Citation information is also provided in the file `CITATION.cff`.

## References
- Liebermeister W. and Noor E. (2021), *Model balancing: a search for in-vivo kinetic constants and consistent
  metabolic states*, Metabolites 11(11):749, [doi:10.3390/metabo11110749](https://doi.org/10.3390/metabo11110749)
  (preprint: [bioRxiv doi:10.1101/2019.12.23.887166](https://www.biorxiv.org/content/10.1101/2019.12.23.887166v2))
- Website: [www.metabolic-economics.de/model-balancing/](https://www.metabolic-economics.de/model-balancing/index.html)
