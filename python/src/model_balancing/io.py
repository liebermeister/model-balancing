import json
from typing import Dict, Union, List

import numpy as np
from sbtab import SBtab
import pandas as pd

from . import Q_


JSON_NAME_MAPPINGS = {
    "Keq": "Keq",
    "Km": "KM",
    "Ka": "KA",
    "Ki": "KI",
    "kcatf": "Kcatf",
    "kcatr": "Kcatr",
}


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

    # Read the kinetic parameters that have 'normal' units:
    for p in ["kcatf", "kcatr", "Km", "Ka", "Ki"]:
        # in the JSON, parameter names have slightly different casing, so we
        # use a dictionary for converting between the conventions.
        p_in_json = JSON_NAME_MAPPINGS[p]
        args[f"{p}_gmean"] = Q_(
            data["kinetic_constants"][p_in_json]["combined"]["geom_mean"],
            data["kinetic_constants"][p_in_json]["unit"],
        )
        args[f"{p}_ln_cov"] = np.array(
            data["kinetic_constants"][p_in_json]["combined"]["cov_ln"]
        )

    # Read the equilibrium constants, which are unitless but if the standard
    # concentration is not 1M, we need to adjust their values based on what
    # it is in the JSON file.
    args["Keq_gmean"] = np.diag(
        np.power(keq_standard_concentration.m_as("M"), args["S"].sum(axis=0))
    ) @ Q_(data["kinetic_constants"]["Keq"]["combined"]["geom_mean"], "")
    args["Keq_ln_cov"] = np.array(
        data["kinetic_constants"]["Keq"]["combined"]["cov_ln"]
    )

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
        columns=["!QuantityType", "!Reaction", "!Compound", "!Mode", "!Unit"],
    ).sort_values("!QuantityType")

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
