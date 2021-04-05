Model balancing
===============

[Model balancing](https://www.metabolic-economics.de/model-balancing/index.html) is a computational method to determine plausible kinetic constants and metabolic states in kinetic metabolic models. It integrates flux, metabolite, protein, and kinetic constant data, using prior distributions for all these variables, and computes the joint posterior mode.
Model balancing can be run in matlab or python. Data tables can be provided in [SBtab](https://www.sbtab.net) format, models can be provided in  [SBML](http://sbml.org) or  [SBtab](https://www.sbtab.net) format.

### Matlab
To run model balancing in matlab, please have a look at the demo scripts in matlab/demo and adjust them to your model and data. For help in matlab, see 'help model-balancing'. For HTML documentation of the matlab code, see matlab/doc.

### Python
To run model balancing in python, your model and input data need to be encoded in a JSON file. To generate this file format from model and data files (in SBML / SBtab formats), you can use the matlab function cvxpy_problem_data_structure.m (as shown in the demo script matlab/demo/demo_cmb_experimental_data.m). For a description of the JSON format structure, see 'help matlab/model-balancing/cmb-cvxpy/cvxpy_problem_data_structure' in matlab.

### Example models
Example model and data files can be found in the folder "resources".

## Dependencies
### Matlab
- [SBML toolbox](http://sbml.org/Software/SBMLToolbox) (optional - needed only if SBML files are used)
- Clone the following [GitHub](https://github.com/liebermeister) repositories
    - [`matlab-utils`](https://github.com/liebermeister/matlab-utils) - utility functions
    - [`metabolic-network-toolbox`](https://github.com/liebermeister/metabolic-network-toolbox) - metabolic network toolbox
    - [`sbtab-matlab`](https://github.com/liebermeister/sbtab-matlab) - SBtab toolbox
    - [`enzyme-cost-minimization`](https://github.com/liebermeister/enzyme-cost-minimization) - enzyme cost minimization toolbox
-  Make sure all the directories and subdirectories are included in your Matlab path
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