Simple example model for model balancing
========================================

The file contains a model structure + artificial data ("measured" kinetic constants and "omics" data for 1 metabolic state), for a toy example. The flux standard deviations are not used in model balancing (but errors caused by the artificial noise in flux data can be assessed)

Later, model (network structure), kinetic "data", and "omics data" may be in separate files (but always as SBtab tables, ie they can be combined into a single file without problems).

Provenance: matlab script 'cmb_example_for_cvxpy.m'

Files:

Model structure (including reconstructed kinetic constants)
* kinetic_model.tsv

Data:
* kinetic_data.tsv
* state_data.tsv

Prior table:
* see ../../prior-table/cmb_prior.tsv

Reconstructed omics data
* metabolic_states.tsv

Diagnostic files:
* options.tsv
* report.txt
* results.mat
