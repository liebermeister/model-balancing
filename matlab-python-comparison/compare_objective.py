#!/usr/bin/env python3
"""Compare the Matlab and Python versions of model balancing at the same point.

usage: compare_objective.py RESULT_DIR PATH_TO_PYTHON_PACKAGE_SRC

The JSON problem files from python-sbtab-to-json are read with the Python package; the Python
objective is then evaluated at the optimum found by Matlab, using the terms that are comparable.
Nothing in the Python package is modified.

Checks that fail because of known issues in the Python package are reported as KNOWN (see
KNOWN_ISSUES). They do not cause a failing exit code, but an unexpected pass is reported, so
that the entry can be removed once the issue is fixed.
"""
import json
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")

LN1000 = np.log(1000.0)
TOLERANCE = 1e-6

# check name -> reason; all of them are issues in the Python package
KNOWN_ISSUES = {
    "read_matlab_json": "reader fails on Matlab files: `true` empty or null; one state gives flat vectors",
    "enzyme_levels_regulated": "ln_eta_regulation in model_balancing_noncvx.py has the wrong sign and inhibitor form",
    "objective_regulated": "follows from enzyme_levels_regulated",
    "beta_definition": "Python uses beta, Matlab and the article use beta^2",
}


def matrix(x):
    a = np.array(x, dtype=float)
    return a.reshape(-1, 1) if a.ndim == 1 else a


def report(results, case, check, passed, detail=""):
    if passed:
        status = "PASS (expected failure: please remove from KNOWN_ISSUES)" if check in EXPECTED.get(case, ()) else "PASS"
    elif check in EXPECTED.get(case, ()):
        status = "KNOWN"
    else:
        status = "FAIL"
    results.append(status)
    print(f"  {check:28s} {status:6s} {detail}")


EXPECTED = {}


def main(result_dir, python_src):
    sys.path.insert(0, python_src)
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python-sbtab-to-json"))
    import sbtab_to_json
    from model_balancing.io import read_arguments_from_json
    from model_balancing.model_balancing_noncvx import ModelBalancing

    results = []
    manifest = json.load(open(os.path.join(result_dir, "manifest.json")))
    for case in manifest:
        name = case["name"]
        regulated = name.startswith("test_model") and "cycle" not in name
        EXPECTED[name] = set()
        EXPECTED[name].add("read_matlab_json")
        EXPECTED[name].add("beta_definition")
        if regulated:
            EXPECTED[name] |= {"enzyme_levels_regulated", "objective_regulated"}
        print(name)

        m = json.load(open(os.path.join(result_dir, name + "_matlab.json")))

        try:
            read_arguments_from_json(json.load(open(os.path.join(result_dir, name + "_problem.json"))))
            report(results, name, "read_matlab_json", True)
        except Exception as e:
            report(results, name, "read_matlab_json", False, repr(e)[:70])

        # the converter writes all arrays in their proper shape; true values are not available
        d = sbtab_to_json.fill_true_values(json.load(open(os.path.join(result_dir, name + "_python.json"))))
        args = read_arguments_from_json(d)

        S = args["S"]
        alpha = float(m["alpha"])
        mb = ModelBalancing(**args, alpha=alpha, beta=0.0)
        point = {
            "ln_Km": np.array(m["lnKM"], dtype=float) - LN1000,
            "ln_Ka": np.array(m["lnKA"], dtype=float) - LN1000 if len(m["lnKA"]) else None,
            "ln_Ki": np.array(m["lnKI"], dtype=float) - LN1000 if len(m["lnKI"]) else None,
            "ln_kcatf": np.array(m["lnKcatf"], dtype=float),
            "ln_Keq": np.array(m["lnKeq"], dtype=float) + S.sum(0) * np.log(1e-3),
            "ln_conc_met": matrix(m["X"]) - LN1000,
        }

        V = matrix(m["V"])
        nonzero = V != 0
        report(results, name, "preposterior_metabolites",
               np.abs(np.log(args["geom_mean_conc_met"].m_as("mM")) - matrix(m["Xmean"])).max() < TOLERANCE
               and np.abs(np.diag(args["precision_ln_conc_met"]) ** -0.5 - matrix(m["Xstd"]).flatten()).max() < TOLERANCE)
        report(results, name, "preposterior_enzymes",
               np.abs(np.log(args["geom_mean_conc_enz"].m_as("mM")) - matrix(m["lnEmean"])).max() < TOLERANCE
               and np.abs(np.diag(args["precision_ln_conc_enz"]) ** -0.5 - matrix(m["lnEstd"]).flatten()).max() < TOLERANCE)

        try:
            mb.extend_var_dict_to_dependent_params(point)
        except Exception as e:
            report(results, name, "enzyme_levels_regulated" if regulated else "enzyme_levels", False, repr(e)[:70])
            continue

        check = "enzyme_levels_regulated" if regulated else "enzyme_levels"
        difference = np.abs(point["ln_conc_enz"][nonzero] - (np.log(matrix(m["E"])[nonzero]) - LN1000)).max()
        report(results, name, check, difference < TOLERANCE, f"max |ln e difference| = {difference:.2e}")

        report(results, name, "kcat_reverse", np.abs(point["ln_kcatr"] - np.array(m["lnKcatr"], dtype=float)).max() < TOLERANCE)

        # Python also counts reactions with zero flux in the enzyme score: this term is constant
        # (ln e = 0 for these reactions) and is removed before comparing the objectives
        z = mb.var_dict_to_z_scores(point)
        gm = np.log(args["geom_mean_conc_enz"].m_as("M"))
        precision = np.diag(args["precision_ln_conc_enz"]).reshape(gm.shape, order="C")
        diff = point["ln_conc_enz"] - gm
        constant = (precision * diff**2 * np.where(diff > 0, 1.0, alpha))[~nonzero].sum()
        python_score = sum(float(v) for v in z.values()) - constant
        matlab_score = -2 * float(m["log_posterior"])
        check = "objective_regulated" if regulated else "objective"
        report(results, name, check, abs(python_score - matlab_score) < 1e-5 * (1 + matlab_score),
               f"Matlab {matlab_score:.6f}, Python {python_score:.6f}")

        mbb = ModelBalancing(**args, alpha=alpha, beta=2.0)
        zb = mbb.var_dict_to_z_scores(point)["c_over_Km"]
        matlab_b = -2 * (float(m["log_posterior_beta_2"]) - float(m["log_posterior"]))
        report(results, name, "beta_definition", abs(zb - matlab_b) < 1e-6 * (1 + matlab_b), f"Matlab {matlab_b:.3f}, Python {zb:.3f}")

    n_fail = results.count("FAIL")
    print(f"\n{len(results)} checks: {results.count('PASS')} passed, {results.count('KNOWN')} known issues, {n_fail} failed")
    for key, reason in KNOWN_ISSUES.items():
        print(f"  known issue {key}: {reason}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
