Example model: E. coli central metabolism
=========================================

Model of E. coli central metabolism from Noor et al. (2016) (see `../README.md` for the reference).

Files defining the model structure:

* `e_coli_noor_2016.xml`: SBML file (network structure)

* `e_coli_noor_2016.tsv`: SBtab file (network structure)

* `e_coli_noor_2016_ConcentrationConstraint.tsv`: SBtab file (metabolite concentration constraints)

* `e_coli_noor_2016_Position.tsv`: SBtab file (layout information for network graphics)

Data files:

* `e_coli_ccm_kinetic_data`: kinetic data (see the README file in that folder)

* `e_coli_ccm_state_data_COPY_FOR_INFO.tsv`: metabolic state data (copy, for information)

* `e_coli_example_reference_state.mat`: reference state for the E. coli model, determined by enzyme cost minimization
  and written by `matlab/cmb_e_coli_c_init.m`. It is only used for the E. coli test model with artificial data,
  as a starting point for generating realistic artificial data.

* `mat-files`: model and data from the enzyme cost minimization toolbox (see the README file in that folder)

Other information

* The matlab script for computing the reference state can be found in the subfolder `matlab`

For running model balancing on this example, please see the demo scripts
`demo_cmb_artificial_data.m` and `demo_cmb_experimental_data.m` in the folder `model-balancing/demo`
(in the main directory of this repository).

## Provenance

* `e_coli_noor_2016.xml`: from the enzyme cost minimization example model for E. coli
  (original filename `ecoli_ccm_ProteinComposition_Network.xml`).
  Note: the information about external metabolites (boundaryCondition=True) in this file is wrong and needs to be fixed.

* `e_coli_noor_2016.tsv`: from the Metabolic Network Toolbox (parameter balancing models, `mnt_parameter_balancing/models`);
  modified: oxaloacetate set to internal.

* `e_coli_ccm_state_data_COPY_FOR_INFO.tsv`: state data used in model balancing runs for this model
  (results folder `e_coli_noor_2016/no_kinetic_data_balanced`).
