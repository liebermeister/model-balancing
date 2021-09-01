Example model: E. coli central metabolism 
==========================================

Files defining the model structure:

* SBML file (network structure): in .xml file

* SBtab files: (network structure, constraints, and layout information in separate .tsv files

Other information

* Matlab scripts for generating the network files can be found in the subfolder `matlab`

* Artificial data can be found in the subfolder `data`

* The network layout can be found in the subfolder `graphics`

For running model balancing on this example, please see `../matlab/demo/demo_cmb_experimental_data.m`


## Provenance:

* file e_coli_noor_2016.xml  from /home/wolfram/projekte/enzyme_cost_minimization/html/data/ecoli_ccm/
  (original filename ecoli_ccm_ProteinComposition_Network.xml)
  -> NOTE WRONG information about external metablites (boundaryCondition=True; needs to be fixed) 

* file e_coli_noor_2016.tsv from MNT repository: mnt/mnt/mnt-extensions/mnt_parameter_balancing/models/
  modified: set oxaloacetate to internal

* file ecoli_ccm_kinetic_data.tsv from /home/wolfram/projekte/model_balancing/matlab/results/e_coli_noor_2016/no_kinetic_data_balanced/data/kinetic_data.tsv

* file ecoli_ccm_state_data.tsv from /home/wolfram/projekte/model_balancing/matlab/results/e_coli_noor_2016/no_kinetic_data_balanced/data/state_data.tsv

* file e_coli_example_reference_state contains a reference state for the e coli model, determined by ECM, and written by cmb_e_coli_c_init.m
  It is only used for the E coli test model with artificial data, as a starting point for generating realistic artificial data.