Simple test model (double branch model)
==========================================

Network structure:

```
X - O       X
     \     /
      O - O 
     /     \
X - O       X
```

Files defining the model structure:

* SBML file (network structure): in .xml file

* SBtab files: (network structure, constraints, and layout information in separate .tsv files

Other information

* Matlab scripts for generating the network files can be found in the subfolder `matlab`

* Artificial data can be found in the subfolder `data`

* The network layout can be found in the subfolder `graphics`

For running model balancing on this example, please see `../matlab/demo/demo_cmb_experimental_data.m`
