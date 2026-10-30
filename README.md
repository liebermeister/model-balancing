Model balancing
===============

[Model balancing](https://www.metabolic-economics.de/model-balancing/index.html) is a computational method to
determine plausible kinetic constants and metabolic states in kinetic metabolic models. It integrates flux,
metabolite, protein, and kinetic constant data, using prior distributions for all these variables, and
computes the joint posterior mode.
Model balancing can be run in matlab or python. Data tables can be provided in [SBtab](https://www.sbtab.net)
format, models can be provided in [SBML](https://sbml.org) or [SBtab](https://www.sbtab.net) format.

This repository contains the matlab implementation, example models and data, and example results.

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
   addpath(genpath('model-balancing/matlab'));
   addpath(genpath('matlab-utils'));
   addpath(genpath('metabolic-network-toolbox'));
   addpath(genpath('sbtab-matlab'));
   addpath(genpath('enzyme-cost-minimization'));
   ```
3. Run one of the demo scripts:
   ```matlab
   demo_cmb_artificial_data      % generate artificial data for a small model and balance it
   demo_cmb_experimental_data    % balance a model with data read from SBtab files
   ```
   By default, output files are written to matlab's `tempdir`.

## Contents of this repository

| Folder          | Contents |
|-----------------|----------|
| `model-balancing` | Matlab code: main function `model_balancing.m`, functions (`model-balancing/model-balancing`), demo scripts (`model-balancing/demo`) |
| `doc`           | HTML documentation of the Matlab code |
| `resource-data` | Example models (`resource-data/models`) and the prior table (`resource-data/data/data-prior`) |
| `examples`      | Input files for the python implementation (`examples/json`), model balancing results (`examples/results`), and result graphics (`examples/graphics`) |

### Matlab
To run model balancing in matlab, please adjust the demo scripts in `model-balancing/demo` to your model and data.
For help, please type `help model-balancing` in matlab, or see the documentation at
[www.metabolic-economics.de/model-balancing/matlab-doc/](https://www.metabolic-economics.de/model-balancing/matlab-doc/)
(also contained in this repository in `doc`).

### Python
For using model balancing in python, please refer to our [Read The Docs](https://model-balancing.readthedocs.io/en/latest/index.html) page;
the code itself can be found on [our GitLab repository](https://gitlab.com/elad.noor/model-balancing).

### Example models
Example model and data files can be found in the folder `resource-data`.
The folder `resource-data/models` contains files for four different example models.
Each model comes as an SBML (`.xml`) and an SBtab (`.tsv`) file.
The folders `examples/results/[MODEL_NAME]` contain balancing results for a large number of model balancing problems
with these models: artificial data for all 4 models, experimental data for the E. coli model (folder `e_coli_noor_2016`).

### How to balance your own model
To balance your own model, please have a look at the matlab demo script `model-balancing/demo/demo_cmb_experimental_data.m`.
You need to prepare three input files in SBtab format, describing

1. the network structure and metabolite constraints (tables Reaction, Compound, Position (optional), Parameter (optional))
2. kinetic data (table ParameterData)
3. state data (tables FluxData, MetaboliteConcentrationData, EnzymeConcentrationData)

For an example, please have a look at the files `artificial_network_true.tsv`, `artificial_kinetic_data.tsv`,
and `artificial_state_data.tsv` in the folder `resource-data/models/branch_point_model/data/`.

To balance your model in matlab, you can then use the code from `demo_cmb_experimental_data.m`.
Please note that there are many more settings that you can change. For an overview, type `help cmb_default_options` in matlab.

To balance your model in python, you first need to convert your model and data into a JSON file.
This can be done in matlab and is shown in `demo_cmb_experimental_data.m`. For an overview of the JSON format,
type `help cvxpy_problem_data_structure` in matlab, or see `examples/json/README.md`.

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

The code was tested with Matlab R2019b on Linux. The optimisation uses `fmincon`
(Optimization Toolbox).

## Tests

From the repository folder, run

```matlab
results = runtests('tests')
```

## Versions

To see the installed version, type `model_balancing_version`. Older versions can be downloaded from the [releases page](https://github.com/liebermeister/model-balancing/releases) of this repository, and changes between versions are listed in `CHANGELOG.md`.

## License
This package is released under the [GNU General Public License, version 2](LICENSE).

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
