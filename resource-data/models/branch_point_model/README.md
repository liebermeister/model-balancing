Simple test model (metabolic branch point)
==========================================

Network structure:

```
  X1
   \v1
    \    v3
      S4 --> X3
    /
   /v2
 X2
```

Files defining the model structure:

* `branch_point_model.xml`: SBML file (network structure)

* `branch_point_model.tsv`: SBtab file (network structure)

* `branch_point_model_ConcentrationConstraint.tsv`: SBtab file (metabolite concentration constraints)

* `branch_point_model_Position.tsv`: SBtab file (layout information for network graphics)

Other information

* Matlab scripts for generating the network files can be found in the subfolder `matlab`

* Artificial data (true model, kinetic data, and state data, in SBtab format) can be found in the subfolder `data`

* A network graphics (`branch_point_model.pdf`) can be found in the subfolder `graphics`

For running model balancing on this example, please see the demo scripts
`demo_cmb_artificial_data.m` and `demo_cmb_experimental_data.m` in the folder `matlab/demo`
(in the main directory of this repository).
