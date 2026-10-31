# Matlab-Python comparison for model balancing

Permanent test that the Matlab version, the converter `python-sbtab-to-json`, and the Python
version of model balancing describe the same problem. Nothing in the Python package is changed.

    MB_PYTHON_SRC=/path/to/python-package/src ./run_comparison.sh /path/to/python

`MB_PYTHON_SRC` is the `src` folder of the [Python package](https://gitlab.com/elad.noor/model-balancing). The
interpreter needs numpy, scipy, pint, cvxpy and sbtab. The script stops on the first failing step and
writes its files to `results`; set `MB_COMPARISON_OUT` for a different output folder. `matlab_side.m`
sets up the Matlab paths itself; set `MB_MATLAB_PACKAGES` if the Matlab packages (Metabolic Network
Toolbox and others) are not in the folder that contains this repository.

## Steps

1. `matlab_side.m` runs nine cases. For each, it writes the JSON problem (`*_problem.json`) and the
   optimum found by Matlab (`*_matlab.json`).
2. `compare_json.py` runs the converter on the same SBtab files, writes `*_python.json`, and compares
   all numbers with the Matlab JSON (tolerance 1e-8, shape ignored).
3. `compare_objective.py` evaluates the Python objective at the Matlab optimum and compares the
   preposterior distributions, enzyme levels, reverse kcat values, the objective, and the term with
   beta.

## Test cases

| Case | Purpose |
|---|---|
| `test_model_alpha_0.5`, `_1`, `_0` | test model, three values of alpha |
| `test_model_no_std` | no error bars in the data file |
| `test_model_only_Keq`, `test_model_no_kinetic_data` | other settings of `use_kinetic_data` |
| `test_model_cycle` | closed cycle, one Wegscheider condition |
| `branch_point_model`, `e_coli_noor_2016` | models from the Matlab package |

The test model in `test_model` comes from `make_test_model.py` (seeded). It has 6 reactions and 3
states, with an activator and an inhibitor, a stoichiometric coefficient of 2, a reaction that
always runs backwards, a reaction that reverses between states, a reaction with zero flux in two
states, and missing metabolite and enzyme data.

## Reading the output

`compare_objective.py` marks each check as:

- `PASS`: Matlab and Python agree
- `KNOWN`: fails because of a listed issue in the Python package; does not fail the test
- `FAIL`: new difference; the exit code is 1
- `PASS (expected failure ...)`: an issue was fixed, so remove it from `KNOWN_ISSUES` and from the
  expected checks

Current state: all 9 cases identical in `compare_json.py`; in `compare_objective.py`, 63 checks with
33 passed, 30 known issues, 0 failed.

The comparison uses the point found by Matlab, so it does not show differences in the constraints.
In particular, Python has no Wegscheider conditions: on `test_model_cycle`, the solver of the Python
package ends at a lower objective than Matlab (0.331 against 0.348), with Keq values that violate
the condition by 0.50 on the log scale.
