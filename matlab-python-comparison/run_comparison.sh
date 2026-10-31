#!/bin/bash
# Compare the Matlab version, the SBtab-to-JSON converter, and the Python version of model balancing.
# usage: run_comparison.sh [PYTHON]
#   PYTHON: Python interpreter with numpy, scipy, pint, cvxpy, sbtab (default: python3)
#   MB_PYTHON_SRC (required): src folder of the Python package (gitlab.com/elad.noor/model-balancing)
# Output files are written to the folder "results"; the Python package is not changed.
set -e
cd "$(dirname "$0")"
PY=${1:-python3}
export PYTHONDONTWRITEBYTECODE=1
export MB_COMPARISON_OUT=${MB_COMPARISON_OUT:-$PWD/results}
PYTHON_SRC=${MB_PYTHON_SRC:?Please set MB_PYTHON_SRC to the src folder of the Python package}

matlab -batch matlab_side
echo; echo "=== Matlab JSON against python-sbtab-to-json ==="
"$PY" compare_json.py "$MB_COMPARISON_OUT"
echo; echo "=== Matlab against Python, at the optimum found by Matlab ==="
"$PY" compare_objective.py "$MB_COMPARISON_OUT" "$PYTHON_SRC"
