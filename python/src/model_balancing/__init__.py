import os
import pint
import warnings
from typing import Dict
from sbtab import SBtab
import numpy as np
import json
import pandas as pd

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
    args["fluxes"] = Q_(
        data["reaction_fluxes"]["data"]["mean"], data["reaction_fluxes"]["unit"]
    )

    args["Keq_gmean"] = Q_(
        data["kinetic_constants"]["Keq"]["combined"]["geom_mean"], ""
    )

    # "fix" the standard concentration (to a convention where it is 1 M which is used internally)
    args["Keq_gmean"] *= np.power(
        keq_standard_concentration.m_as("M"), args["S"].sum(axis=0)
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

    args["conc_met_gmean"] = Q_(
        data["metabolite_concentrations"]["combined"]["geom_mean"],
        data["metabolite_concentrations"]["unit"],
    )
    args["conc_met_gstd"] = np.array(
        data["metabolite_concentrations"]["combined"]["geom_std"]
    )

    args["conc_enz_gmean"] = Q_(
        data["enzyme_concentrations"]["combined"]["geom_mean"],
        data["metabolite_concentrations"]["unit"],
    )
    args["conc_enz_gstd"] = np.array(
        data["enzyme_concentrations"]["combined"]["geom_std"]
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


from .model_balancing import ModelBalancing
from .model_balancing_cvx import ModelBalancingConvex

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
