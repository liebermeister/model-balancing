Simple test model (chain of three reactions)
============================================

Network structure:

```
X1 <-> S2 <-> S3 <-> X4
```

Files defining the model structure:

* `three_chain_model.xml`: SBML file (network structure)

* `three_chain_model.tsv`: SBtab file (network structure)

* `three_chain_model_ConcentrationConstraint.tsv`: SBtab file (metabolite concentration constraints)

* `three_chain_model_Position.tsv`: SBtab file (layout information for network graphics)

Other information

* Matlab scripts for generating the network files can be found in the subfolder `matlab`

* Artificial data (true model, kinetic data, and state data, in SBtab format) can be found in the subfolder `data`

* A network graphics (`three_chain_model.pdf`) can be found in the subfolder `graphics`

For running model balancing on this example, please see the demo scripts
`demo_cmb_artificial_data.m` and `demo_cmb_experimental_data.m` in the folder `model-balancing/demo`
(in the main directory of this repository).
