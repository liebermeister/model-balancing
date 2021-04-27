import json
import os
import warnings
from typing import Dict, List, Union

import numpy as np
import pandas as pd
import pint
from sbtab import SBtab

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


def standardize_input_matrix(x: Union[List[float], np.array], unit: str) -> Q_:
    """Create a 2D Numpy array of Quantities.

    If the input is only 1D, make sure it becomes a 2D array with a single
    column. This is important because we always assume that all our inputs are
    2D.
    """
    if not x:
        return Q_(x, unit)
    elif type(x) == list:
        return Q_(x, unit).reshape(len(x), -1)
    else:
        return Q_(x, unit).reshape(x.shape[0], -1)


def read_arguments_json(
    json_fname: str,
) -> Dict[str, np.array]:
    with open(json_fname, "rt") as fp:
        data = json.load(fp)

    keq_standard_concentration = Q_(data["standard_concentration"])

    args = {}
    args["S"] = np.array(data["network"]["stoichiometric_matrix"])
    args["A_act"] = np.array(data["network"]["activation_matrix"])
    args["A_inh"] = np.array(data["network"]["inhibition_matrix"])
    args["fluxes"] = standardize_input_matrix(
        data["reaction_fluxes"]["data"]["mean"], data["reaction_fluxes"]["unit"]
    )

    args["Keq_gmean"] = Q_(
        data["kinetic_constants"]["Keq"]["combined"]["geom_mean"], ""
    )

    # "fix" the standard concentration (to a convention where it is 1 M which is used internally)
    args["Keq_gmean"] = (
        np.diag(np.power(keq_standard_concentration.m_as("M"), args["S"].sum(axis=0)))
        @ args["Keq_gmean"]
    )

    args["Keq_ln_cov"] = np.array(
        data["kinetic_constants"]["Keq"]["combined"]["cov_ln"]
    )

    args["kcatf_gmean"] = Q_(
        data["kinetic_constants"]["Kcatf"]["combined"]["geom_mean"],
        data["kinetic_constants"]["Kcatf"]["unit"],
    )
    args["kcatf_ln_cov"] = np.array(
        data["kinetic_constants"]["Kcatf"]["combined"]["cov_ln"]
    )

    args["kcatr_gmean"] = Q_(
        data["kinetic_constants"]["Kcatr"]["combined"]["geom_mean"],
        data["kinetic_constants"]["Kcatr"]["unit"],
    )
    args["kcatr_ln_cov"] = np.array(
        data["kinetic_constants"]["Kcatr"]["combined"]["cov_ln"]
    )

    args["Km_gmean"] = Q_(
        data["kinetic_constants"]["KM"]["combined"]["geom_mean"], "mM"
    )  # data["kinetic_constants"]["KM"]["unit"])
    args["Km_ln_cov"] = np.array(data["kinetic_constants"]["KM"]["combined"]["cov_ln"])

    args["Ka_gmean"] = Q_(
        data["kinetic_constants"]["KA"]["combined"]["geom_mean"], "mM"
    )  # data["kinetic_constants"]["KA"]["unit"])
    args["Ka_ln_cov"] = np.array(data["kinetic_constants"]["KA"]["combined"]["cov_ln"])

    args["Ki_gmean"] = Q_(
        data["kinetic_constants"]["KI"]["combined"]["geom_mean"], "mM"
    )  # data["kinetic_constants"]["KI"]["unit"])
    args["Ki_ln_cov"] = np.array(data["kinetic_constants"]["KI"]["combined"]["cov_ln"])

    args["conc_met_gmean"] = standardize_input_matrix(
        data["metabolite_concentrations"]["combined"]["geom_mean"],
        data["metabolite_concentrations"]["unit"],
    ).reshape(args["S"].shape[0], args["fluxes"].shape[1])
    args["conc_met_ln_cov"] = (
        np.diag(
            np.log(
                np.array(
                    data["metabolite_concentrations"]["combined"]["geom_std"]
                ).T.flatten()
            )
        )
        ** 2.0
    )

    args["conc_enz_gmean"] = standardize_input_matrix(
        data["enzyme_concentrations"]["combined"]["geom_mean"],
        data["metabolite_concentrations"]["unit"],
    ).reshape(args["S"].shape[1], args["fluxes"].shape[1])
    args["conc_enz_ln_cov"] = (
        np.diag(
            np.log(
                np.array(
                    data["enzyme_concentrations"]["combined"]["geom_std"]
                ).T.flatten()
            )
        )
        ** 2.0
    )

    args["metabolite_names"] = data["network"]["metabolite_names"]
    args["reaction_names"] = data["network"]["reaction_names"]
    args["state_names"] = data["state_names"]
    args["rate_law"] = "CM"
    return args


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


def to_state_sbtab(
    v,
    c,
    e,
    delta_g,
    metabolite_names,
    reaction_names,
    state_names,
) -> SBtab.SBtabDocument:

    state_sbtabdoc = SBtab.SBtabDocument(name="MB result")
    flux_df = pd.DataFrame(v.m_as("mM/s"), columns=state_names)
    flux_df.insert(0, "!QuantityType", "rate of reaction")
    flux_df.insert(1, "!Reaction", reaction_names)
    flux_sbtab = SBtab.SBtabTable.from_data_frame(
        flux_df.astype(str),
        table_id="Flux",
        table_name="Metabolic fluxes",
        table_type="QuantityMatrix",
        unit="mM/s",
    )
    state_sbtabdoc.add_sbtab(flux_sbtab)

    conc_met_df = pd.DataFrame(c.m_as("mM"), columns=state_names)
    conc_met_df.insert(0, "!QuantityType", "concentration")
    conc_met_df.insert(1, "!Compound", metabolite_names)
    conc_met_sbtab = SBtab.SBtabTable.from_data_frame(
        conc_met_df.astype(str),
        table_id="MetaboliteConcentration",
        table_name="Metabolite concentration",
        table_type="QuantityMatrix",
        unit="mM",
    )
    state_sbtabdoc.add_sbtab(conc_met_sbtab)

    conc_enz_df = pd.DataFrame(e.m_as("mM"), columns=state_names)
    conc_enz_df.insert(0, "!QuantityType", "concentration of enzyme")
    conc_enz_df.insert(1, "!Reaction", reaction_names)
    conc_enz_sbtab = SBtab.SBtabTable.from_data_frame(
        conc_enz_df.astype(str),
        table_id="EnzymeConcentration",
        table_name="Enzyme concentration",
        table_type="QuantityMatrix",
        unit="mM",
    )
    state_sbtabdoc.add_sbtab(conc_enz_sbtab)

    gibbs_energy_df = pd.DataFrame(delta_g.m_as("kJ/mol"), columns=state_names)
    gibbs_energy_df.insert(0, "!QuantityType", "Gibbs energy of reaction")
    gibbs_energy_df.insert(1, "!Reaction", reaction_names)
    gibbs_energy_sbtab = SBtab.SBtabTable.from_data_frame(
        gibbs_energy_df.astype(str),
        table_id="ReactionGibbsFreeEnergy",
        table_name="Gibbs free energies of reaction",
        table_type="QuantityMatrix",
        unit="kJ/mol",
    )
    state_sbtabdoc.add_sbtab(gibbs_energy_sbtab)
    return state_sbtabdoc


def to_model_sbtab(
    kcatf,
    kcatr,
    Keq,
    Km,
    Ka,
    Ki,
    S,
    A_act,
    A_inh,
    metabolite_names,
    reaction_names,
    state_names,
) -> SBtab.SBtabDocument:

    model_sbtabdoc = SBtab.SBtabDocument(name="MB result")

    parameter_data = []
    parameter_data += [
        ("equilibrium constant", rxn, "", value, "dimensionless")
        for rxn, value in zip(reaction_names, Keq.m_as(""))
    ]
    parameter_data += [
        ("catalytic rate constant geometric mean", rxn, "", value.m_as("1/s"), "1/s")
        for rxn, value in zip(reaction_names, (kcatf * kcatr) ** (0.5))
    ]
    parameter_data += [
        ("forward catalytic rate constant", rxn, "", value.m_as("1/s"), "1/s")
        for rxn, value in zip(reaction_names, kcatf)
    ]
    parameter_data += [
        ("reverse catalytic rate constant", rxn, "", value.m_as("1/s"), "1/s")
        for rxn, value in zip(reaction_names, kcatr)
    ]
    for j, rxn in enumerate(reaction_names):
        for i, met in enumerate(metabolite_names):
            if S[i, j] != 0:
                parameter_data += [
                    ("Michaelis constant", rxn, met, Km[i, j].m_as("mM"), "mM")
                ]
            if A_act[i, j] != 0:
                parameter_data += [
                    ("Activation constant", rxn, met, Ka[i, j].m_as("mM"), "mM")
                ]
            if A_inh[i, j] != 0:
                parameter_data += [
                    ("Inhibition constant", rxn, met, Ki[i, j].m_as("mM"), "mM")
                ]

    parameter_df = pd.DataFrame(
        data=parameter_data,
        columns=["QuantityType", "Reaction", "Compound", "Mode", "Unit"],
    ).sort_values("QuantityType")

    parameter_sbtab = SBtab.SBtabTable.from_data_frame(
        parameter_df.astype(str),
        table_id="Parameter",
        table_name="Parameter",
        table_type="Quantity",
        unit="",
    )
    parameter_sbtab.change_attribute("StandardConcentration", "M")
    model_sbtabdoc.add_sbtab(parameter_sbtab)

    return model_sbtabdoc


from ._version import get_versions
from .model_balancing import ModelBalancing
from .model_balancing_cvx import ModelBalancingConvex

__version__ = get_versions()["version"]
del get_versions
