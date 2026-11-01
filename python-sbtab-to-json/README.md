# sbtab_to_json

Converts the input files of model balancing (SBtab) into the JSON problem file, so that the
Python version of model balancing can start from the same files as the Matlab version, without
Matlab. It is a stand-alone port of this Matlab code path:

    mb_load_and_convert_model_and_data -> mb_model_and_data
      -> mb_prepare_posterior -> cvxpy_problem_data_structure

It needs Python 3 and numpy. The prior table is bundled in `data/mb_prior.tsv`.
It does not change the Python package, and it is written so that it can be moved into it.

## Usage

    python sbtab_to_json.py network.tsv kinetic_data.tsv state_data.tsv out.json \
        --constraints constraints.tsv --mean S1_Mean S2_Mean --std S1_Std S2_Std \
        --fill-true-values

- `--mean`, `--std`: columns in the state data file, one for each metabolic state. The default
  for `--mean` is a single column `Mean`. Without `--std`, the default data errors of the prior
  table are used.
- `--use-kinetic-data`: `all` (default), `only_Keq`, or `none`
- `--match-data-by`: `ModelElementId` (default), `KeggId`, or `BiggId`
- `--constraints`: SBtab file with the table `ConcentrationConstraint`
- `--prior-file`: use a different prior table
- `--fill-true-values`: needed for now, see below

From Python:

    from sbtab_to_json import build_problem, fill_true_values, write_json
    problem = build_problem("network.tsv", "kinetic_data.tsv", "state_data.tsv",
                            columns_mean=["S1_Mean", "S2_Mean"])
    write_json(fill_true_values(problem), "out.json")

Tables are found by their `TableID`, so one file may contain several of them.

## Input tables

| File | Tables |
|---|---|
| network | `Reaction` (column `ReactionFormula`, optionally `MetabolicRegulation`), `Compound`, optionally `EnzymeRegulation` |
| kinetic data | `ParameterData` or `Parameter` |
| state data | `MetaboliteConcentrationData`, `EnzymeConcentrationData`, `FluxData` |
| constraints | `ConcentrationConstraint` |

## Validation

`../matlab-python-comparison` compares the output with the JSON written by the fixed Matlab code.
All numbers agree within 1e-8 for nine cases: a test model with regulation, missing data and a zero
flux (different options), a model with a closed cycle, the branch point model, and the *E. coli*
model. Only array shapes differ, see below.

The reader of the Python package reads the output, and its solver runs on it without Matlab. On
the branch point model it reaches an objective of 0.08245 (Matlab: 0.08250).

## Differences from the Matlab output

- Arrays keep their shape. Matlab's `jsonencode` turns (N, 1) matrices into flat vectors and
  vectors of length one into numbers, which makes problems with one state unreadable for the
  Python reader.
- Concentration constraints are matched by metabolite. Matlab uses the row order.

The converter reproduces the Matlab code after the fixes in `mb_prepare_posterior` and
`state_data_to_data` of October 2026. It does not reproduce older Matlab versions or the example
JSON files made with them.

## `--fill-true-values`

Real data have no true values, but `read_arguments_from_json` in the Python package fails if they
are missing. This option inserts the preposterior means in their place. They are not used in the
optimisation. The option can go once the reader accepts missing values.
