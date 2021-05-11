import json
import os
import warnings
from typing import Dict, List, Union

import numpy as np
import pandas as pd
import pint

# Disable Pint's old fallback behavior (must come before importing Pint)
os.environ["PINT_ARRAY_PROTOCOL_FALLBACK"] = "0"

ureg = pint.UnitRegistry(system="mks")
Q_ = ureg.Quantity

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Q_([])

RT = Q_(8.31e-3 * 298.15, "kJ/mol")

MIN_DRIVING_FORCE = 1e-3 * RT
MIN_FLUX = Q_(1e-9, "M/s")


def convert_to_irreversible(args: Dict[str, np.array]) -> Dict[str, np.array]:

    # in the python implementation of Model Balancing, only non-negative fluxes
    # are allowed. In order to accomodate that, we add a new reaction for each
    # one that has a negative flux (in at least one condition) which is the
    # opposite reaction of the one in the model, and use that whenever the flux
    # is negative. Note that the Keq, kcatr, and kcatf have to be adapted to
    # the new direction

    Nr = args["S"].shape[1]
    v_pos = abs(args["fluxes"])
    v_neg = abs(args["fluxes"])
    v_pos[args["fluxes"] < 0] = Q_(0.0, "mM/s")
    v_neg[args["fluxes"] > 0] = Q_(0.0, "mM/s")

    new_args = dict()

    for k in ["conc_met_gmean", "conc_met_ln_cov", "metabolite_names", "state_names"]:
        new_args[k] = args[k]

    new_args["S"] = np.block([[args["S"], -args["S"]]])
    new_args["A_act"] = np.block([[args["A_act"], args["A_act"]]])
    new_args["A_inh"] = np.block([[args["A_inh"], args["A_inh"]]])
    new_args["fluxes"] = np.vstack([v_pos, v_neg])
    new_args["Keq_gmean"] = np.block([args["Keq_gmean"], 1.0 / args["Keq_gmean"]])
    new_args["Keq_ln_cov"] = np.block(
        [
            [args["Keq_ln_cov"], args["Keq_ln_cov"].T],
            [args["Keq_ln_cov"].T, args["Keq_ln_cov"]],
        ]
    )
    new_args["kcatf_gmean"] = np.block([args["kcatf_gmean"], args["kcatr_gmean"]])
    new_args["kcatf_ln_cov"] = np.block(
        [
            [args["kcatf_ln_cov"], np.zeros((Nr, Nr))],
            [np.zeros((Nr, Nr)), args["kcatr_ln_cov"]],
        ]
    )
    new_args["kcatr_gmean"] = np.block([args["kcatr_gmean"], args["kcatf_gmean"]])
    new_args["kcatr_ln_cov"] = np.block(
        [
            [args["kcatr_ln_cov"], np.zeros((Nr, Nr))],
            [np.zeros((Nr, Nr)), args["kcatf_ln_cov"]],
        ]
    )

    for p in ["Km", "Ka", "Ki"]:
        new_args[f"{p}_gmean"] = np.block([args[f"{p}_gmean"], args[f"{p}_gmean"]])
        new_args[f"{p}_ln_cov"] = np.block(
            [
                [args[f"{p}_ln_cov"], args[f"{p}_ln_cov"].T],
                [args[f"{p}_ln_cov"].T, args[f"{p}_ln_cov"]],
            ]
        )

    new_args[f"conc_enz_gmean"] = np.vstack([args[f"conc_enz_gmean"]] * 2)
    new_args[f"conc_enz_ln_cov"] = np.block(
        [
            [args[f"conc_enz_ln_cov"], args[f"conc_enz_ln_cov"].T],
            [args[f"conc_enz_ln_cov"].T, args[f"conc_enz_ln_cov"]],
        ]
    )

    new_args["reaction_names"] = args["reaction_names"] + [
        rxn + "_r" for rxn in args["reaction_names"]
    ]

    return new_args


from ._version import get_versions
from .model_balancing import ModelBalancing
from .model_balancing_cvx import ModelBalancingConvex

__version__ = get_versions()["version"]
del get_versions
