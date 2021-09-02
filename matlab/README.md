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

Matlab scripts for generating the example models can be found in the folder `../resources/models`

## Implementation notes

* For numerical reasons, the matlab implementation does not tolerate equal lower and upper bounds for metabolites. If c_{min}==c_{max},  10^{-5} is added automatically to ln c_{max}. This also speeds up the calculation.

