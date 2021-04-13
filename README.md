Model balancing
===============

[Model balancing](https://www.metabolic-economics.de/model-balancing/index.html) is a computational method to determine plausible kinetic constants and metabolic states in kinetic metabolic models. It integrates flux, metabolite, protein, and kinetic constant data, using prior distributions for all these variables, and computes the joint posterior mode.
Model balancing can be run in matlab or python. Data tables can be provided in [SBtab](https://www.sbtab.net) format, models can be provided in  [SBML](http://sbml.org) or  [SBtab](https://www.sbtab.net) format.

### Matlab
To run model balancing in matlab, please have a look at the demo scripts in matlab/demo and adjust them to your model and data. For HTML documentation of the matlab code, see matlab/doc. For help in matlab, please type 'help model-balancing'.

### Python and JSON file format
To run model balancing in python, your model and input data need to be encoded in a JSON file. To generate this file format from model and data files (in SBML / SBtab formats), you can use the matlab function cvxpy_problem_data_structure.m (as shown in the demo script matlab/demo/demo_cmb_experimental_data.m). For a description of the JSON format, type 'help cvxpy_problem_data_structure' in matlab.

### Example models
Example model and data files can be found in the folder "resources". The folder resources/models contains files for four different example models. Each model comes as an SBML (.xml) and an SBtab (.tsv) file. The folders cvxpy/examples/[MODEL_NAME] contains balancing results for a large number of model balancing problems with these models: artificial data for all 4 models, experimental data for the E coli model (folder e_coli_noor_2016).

### How to balance your own model
To balance your own model, we recommend to get acquainted with the matlab demo script demo_cmb_experimental_data.m. You need to prepare three input files in SBtab format, describing

1. the network structure and metabolite constraints (tables Reaction, Compound, Position (optional), Parameter (optional))
2. kinetic data (table ParameterData)
3. state data (tables MetabolicFluxData, MetaboliteConcentrationData, EnzymeConcentrationData)

For an example, please have a look at the files 'artificial_network_true.tsv', 'artificial_kinetic_data.tsv', and 'artificial_state_data.tsv' in the folder ./resourcedir/models/branch_point_model/data/.

To balance your model in matlab, you can then use the code from demo_cmb_experimental_data.m. Please note that there are many more settings that you can change. For an overview, type 'help cmb_default_options' in matlab.

To balance your model in python, you first need to convert your model and data into a JSON file. This can be done in matlab and is shown in demo_cmb_experimental_data.m. For an overview of the JSON format, type 'help cvxpy_problem_data_structure' in matlab.

## Dependencies
### Matlab
- [SBML toolbox](http://sbml.org/Software/SBMLToolbox) (optional - needed only if SBML files are used; version needs to match your matlab version and  requires a matching version of libSBML)
- Clone the following [GitHub](https://github.com/liebermeister) repositories
    - [`matlab-utils`](https://github.com/liebermeister/matlab-utils) - utility functions
    - [`metabolic-network-toolbox`](https://github.com/liebermeister/metabolic-network-toolbox) - metabolic network toolbox
    - [`sbtab-matlab`](https://github.com/liebermeister/sbtab-matlab) - SBtab toolbox
    - [`enzyme-cost-minimization`](https://github.com/liebermeister/enzyme-cost-minimization) - enzyme cost minimization toolbox
-  Make sure all the directories and subdirectories are included in your Matlab path
The code was tested with Matlab R2019b on Linux. 
### Python
- Install using `pip`:
    - `cvxpy` ~= 1.1
    - `pint` ~= 0.16
    - `numpy` ~= 1.20
    - `scipy` ~= 1.6
    - `pandas` ~= 1.2
    - `Mosek` ~= 9.2

## License
This package is released under the [GNU General Public License](LICENSE).

## Contact
Please contact [Wolfram Liebermeister](mailto:wolfram.liebermeister@gmail.com) and [Elad Noor](mailto:elad.noor@weizmann.ac.il) with any questions or comments.

## References
Liebermeister W. and Noor E. (2021), *Model balancing: in search of consistent metabolic states and in-vivo kinetic constants*
[bioRxiv doi:10.1101/2019.12.23.887166v2](https://www.biorxiv.org/content/10.1101/2019.12.23.887166v2)

[www.metabolic-economics.de/model-balancing/](https://www.metabolic-economics.de/model-balancing/index.html)