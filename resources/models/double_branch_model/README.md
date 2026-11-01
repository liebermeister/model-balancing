Simple test model (double branch model)
=======================================

Network structure:

```
X - O       X
     \     /
      O - O
     /     \
X - O       X
```

Files defining the model structure:

* `double_branch_model.xml`: SBML file (network structure)

* `double_branch_model.tsv`: SBtab file (network structure)

* `double_branch_model_ConcentrationConstraint.tsv`: SBtab file (metabolite concentration constraints)

* `double_branch_model_Position.tsv`: SBtab file (layout information for network graphics)

Other information

* Matlab scripts for generating the network files can be found in the subfolder `matlab`

* Artificial data (true model, kinetic data, and state data, in SBtab format) can be found in the subfolder `data`

* A network graphics (`double_branch_model.pdf`) can be found in the subfolder `graphics`

For running model balancing on this example, please see the demo scripts
`demo_mb_artificial_data.m` and `demo_mb_experimental_data.m` in the folder `model-balancing/demo`
(in the main directory of this repository).
