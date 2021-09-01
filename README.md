Model balancing
===============

[Model balancing](https://www.metabolic-economics.de/model-balancing/index.html) is a computational method to
determine plausible kinetic constants and metabolic states in kinetic metabolic models. It integrates flux,
metabolite, protein, and kinetic constant data, using prior distributions for all these variables, and
computes the joint posterior mode.
Model balancing can be run in matlab or python. Data tables can be provided in [SBtab](https://www.sbtab.net)
format, models can be provided in  [SBML](http://sbml.org) or  [SBtab](https://www.sbtab.net) format.

### Matlab
To run model balancing in matlab, please have a look at the demo scripts in `matlab/demo` and adjust them to your model and data.
For HTML documentation of the matlab code, see `matlab/doc`. For help in matlab, please type `help model-balancing`.

### Python
For using model balancing in python, please refer to our [Read The Docs](https://model-balancing.readthedocs.io/en/latest/index.html) page,
and the code itself is found on [our GitLab repository](https://gitlab.com/elad.noor/model-balancing)

### Example models
Example model and data files can be found in the folder "resources". 
The folder `resources/models` contains files for four different example models. 
Each model comes as an SBML (`.xml`) and an SBtab (`.tsv`) file. 
The folders `examples/results/[MODEL_NAME]` contains balancing results for a large number of model balancing problems
with these models: artificial data for all 4 models, experimental data for the E coli model (folder `e_coli_noor_2016`).

### How to balance your own model
To balance your own model, please have a look at the matlab demo script `demo_cmb_experimental_data.m`.
You need to prepare three input files in SBtab format, describing

1. the network structure and metabolite constraints (tables Reaction, Compound, Position (optional), Parameter (optional))
2. kinetic data (table ParameterData)
3. state data (tables MetabolicFluxData, MetaboliteConcentrationData, EnzymeConcentrationData)

For an example, please have a look at the files `artificial_network_true.tsv`, `artificial_kinetic_data.tsv`,
and `artificial_state_data.tsv` in the folder `./resourcedir/models/branch_point_model/data/`.

To balance your model in matlab, you can then use the code from `demo_cmb_experimental_data.m`.
Please note that there are many more settings that you can change. For an overview, type `help cmb_default_options` in matlab.

To balance your model in python, you first need to convert your model and data into a JSON file.
This can be done in matlab and is shown in `demo_cmb_experimental_data.m`. For an overview of the JSON format,
type `help cvxpy_problem_data_structure` in matlab.

## Installation
The following installation instructions are only for the Matlab version.
For Python, please refer to our [Read The Docs](https://model-balancing.readthedocs.io/en/latest/index.html) page.

- [SBML toolbox](http://sbml.org/Software/SBMLToolbox) (optional - needed only if SBML files are used;
  version needs to match your matlab version and  requires a matching version of libSBML)
- Clone the following [GitHub](https://github.com/liebermeister) repositories
    - [`matlab-utils`](https://github.com/liebermeister/matlab-utils) - utility functions
    - [`metabolic-network-toolbox`](https://github.com/liebermeister/metabolic-network-toolbox) - metabolic network toolbox
    - [`sbtab-matlab`](https://github.com/liebermeister/sbtab-matlab) - SBtab toolbox
    - [`enzyme-cost-minimization`](https://github.com/liebermeister/enzyme-cost-minimization) - enzyme cost minimization toolbox
- Make sure all the directories and subdirectories are included in your Matlab path.

The code was tested with Matlab R2019b on Linux. 

## License
This package is released under the [GNU General Public License](LICENSE).

## Contact
Please contact [Wolfram Liebermeister](mailto:wolfram.liebermeister@gmail.com)
and [Elad Noor](mailto:elad.noor@weizmann.ac.il) with any questions or comments.

## References
Liebermeister W. and Noor E. (2021), *Model balancing: in search of consistent
metabolic states and in-vivo kinetic constants*
[bioRxiv doi:10.1101/2019.12.23.887166v2](https://www.biorxiv.org/content/10.1101/2019.12.23.887166v2)

[www.metabolic-economics.de/model-balancing/](https://www.metabolic-economics.de/model-balancing/index.html)
