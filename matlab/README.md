Matlab functions for model balancing
====================================

For documentation, see www.metabolic-economics.de/model-balancing/matlab-doc/

## Contents of this folder

Main function

* `model_balancing.m` (calls functions in subfolder `model-balancing`)

Subfolders

* `model-balancing` Matlab functions for model balancing. 

* `demo` Demo scripts for running model balancing with artificial or experimental data

* `doc` HTML source code for documentation at www.metabolic-economics.de/model-balancing/matlab-doc/

Matlab scripts for generating the example models can be found in the folder `../resource-data/models`

## Implementation notes

* The code has been tested with Matlab R2019b on Linux.

* Optimisation is performed by the fmincon function with the interior-point method.

* For numerical reasons, the matlab implementation does not tolerate equal lower and upper bounds for metabolites. If c_min == c_max,  10^-5 is added automatically to ln c_max. This also speeds up the calculation.

* By default, the algorithm starts by running model balancing on an average model state (with metabolic state data given by the geometric mean over the metabolic states). The resulting kinetic constants are used as initial values for the following full calculation with several metabolic states.