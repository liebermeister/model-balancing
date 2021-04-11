import os
import warnings
from typing import Union, List
import json
import numpy as np
from model_balancing import ModelBalancing, ModelBalancingConvex, Q_
import cvxpy as cp
import pandas as pd
from sbtab import SBtab

def to_sbtab(
        mb: Union[ModelBalancing, ModelBalancingConvex],
        metabolite_names: List[str],
        reaction_names: List[str],
        state_names: List[str],
) -> SBtab.SBtabDocument:
    condition_names = ["S1", "S2", "S3", "S4"]

    sbtabdoc = SBtab.SBtabDocument(name="balanced_states")

    flux_df = pd.DataFrame(
        mb.fluxes.m_as("mM/s"),
        columns=condition_names
    )
    flux_df.insert(0, "!QuantityType", "rate of reaction")
    flux_df.insert(1, "!Reaction", reaction_names)
    flux_sbtab = SBtab.SBtabTable.from_data_frame(
        flux_df.astype(str),
        table_id="MetabolicFlux",
        table_name="MetabolicFlux",
        table_type="QuantityMatrix",
        unit="mM/s",
    )
    sbtabdoc.add_sbtab(flux_sbtab)

    if isinstance(mb, ModelBalancingConvex):
        c = Q_(np.exp(mb.ln_conc_met.value), "M")
        e = Q_(np.exp(mb.ln_conc_enz.value), "M")
        driving_forces = -mb.driving_forces.value
        kcatf = Q_(np.exp(mb.ln_kcatf.value), "1/s")
        kcatr = Q_(np.exp(mb.ln_kcatr.value), "1/s")
        Keq = Q_(np.exp(mb.ln_Keq.value), "")
        if mb.ln_Km_gmean.size != 0:
            Km = Q_(np.exp(mb._create_dense_matrix(mb.S, mb.ln_Km).value), "M")
        else:
            Km = Q_(np.exp(mb._create_dense_matrix(mb.S, mb.ln_Km)), "M")
        if mb.ln_Ka_gmean.size != 0:
            Ka = Q_(np.exp(mb._create_dense_matrix(mb.A_act, mb.ln_Ka).value), "M")
        else:
            Ka = Q_(np.exp(mb._create_dense_matrix(mb.A_act, mb.ln_Ka)), "M")
        if mb.ln_Ki_gmean.size != 0:
            Ki = Q_(np.exp(mb._create_dense_matrix(mb.A_inh, mb.ln_Ki).value), "M")
        else:
            Ki = Q_(np.exp(mb._create_dense_matrix(mb.A_inh, mb.ln_Ki)), "M")
    else:
        c = Q_(np.exp(mb.ln_conc_met), "M")
        e = Q_(np.exp(mb.ln_conc_enz), "M")
        driving_forces = -mb.driving_forces
        kcatf = Q_(np.exp(mb.ln_kcatf), "1/s")
        kcatr = Q_(np.exp(mb.ln_kcatr), "1/s")
        Keq = Q_(np.exp(mb.ln_Keq), "")
        Km = Q_(np.exp(mb._create_dense_matrix(mb.S, mb.ln_Km)), "M")
        Ka = Q_(np.exp(mb._create_dense_matrix(mb.A_act, mb.ln_Ka)), "M")
        Ki = Q_(np.exp(mb._create_dense_matrix(mb.A_inh, mb.ln_Ki)), "M")

    conc_met_df = pd.DataFrame(
        c.m_as("mM"),
        columns=condition_names
    )
    conc_met_df.insert(0, "!QuantityType", "concentration")
    conc_met_df.insert(1, "!Compound", ["X1", "S2", "X3", "X4"])
    conc_met_sbtab = SBtab.SBtabTable.from_data_frame(
        conc_met_df.astype(str),
        table_id="MetaboliteConcentrations",
        table_name="Metabolite concentrations",
        table_type="QuantityMatrix",
        unit="mM",
    )
    sbtabdoc.add_sbtab(conc_met_sbtab)

    conc_enz_df = pd.DataFrame(
        e.m_as("mM"),
        columns=condition_names
    )
    conc_enz_df.insert(0, "!QuantityType", "concentration of enzyme")
    conc_enz_df.insert(1, "!Reaction", reaction_names)
    conc_enz_sbtab = SBtab.SBtabTable.from_data_frame(
        conc_enz_df.astype(str),
        table_id="EnzymeConcentrations",
        table_name="Enzyme concentrations",
        table_type="QuantityMatrix",
        unit="mM",
    )
    sbtabdoc.add_sbtab(conc_enz_sbtab)

    k_eq_df = pd.DataFrame(
        Keq.m_as(""),
        columns=["Keq"]
    )
    k_eq_df.insert(0, "!QuantityType", "Reaction equilibrium constant")
    k_eq_df.insert(1, "!Reaction", reaction_names)
    gibbs_energy_sbtab = SBtab.SBtabTable.from_data_frame(
        k_eq_df.astype(str),
        table_id="ReactionEquilibriumConstant",
        table_name="Equilibrium constant of reaction",
        table_type="QuantityMatrix",
        unit="",
    )
    sbtabdoc.add_sbtab(gibbs_energy_sbtab)

    gibbs_energy_df = pd.DataFrame(
        driving_forces,
        columns=condition_names
    )
    gibbs_energy_df.insert(0, "!QuantityType", "Gibbs energy of reaction")
    gibbs_energy_df.insert(1, "!Reaction", reaction_names)
    gibbs_energy_sbtab = SBtab.SBtabTable.from_data_frame(
        gibbs_energy_df.astype(str),
        table_id="ReactionGibbsFreeEnergy",
        table_name="Gibbs free energies of reaction",
        table_type="QuantityMatrix",
        unit="kJ/mol",
    )
    sbtabdoc.add_sbtab(gibbs_energy_sbtab)

    kcat_df = pd.DataFrame(
        np.vstack([kcatf.m_as("1/s"), kcatr.m_as("1/s")]).T,
        columns=["kcatf", "kcatr"]
    )
    kcat_df.insert(0, "!QuantityType", "catalytic rate constant")
    kcat_df.insert(1, "!Reaction", reaction_names)
    kcat_sbtab = SBtab.SBtabTable.from_data_frame(
        kcat_df.astype(str),
        table_id="CatalyticRateConstant",
        table_name="Catalytic rate constant",
        table_type="QuantityMatrix",
        unit="1/s",
    )
    sbtabdoc.add_sbtab(kcat_sbtab)

    Km_df = pd.DataFrame(
        Km.m_as("mM").T,
        columns=metabolite_names
    )
    Km_df.insert(0, "!QuantityType", "Michaelis constant")
    Km_df.insert(1, "!Reaction", reaction_names)
    Km_sbtab = SBtab.SBtabTable.from_data_frame(
        Km_df.astype(str),
        table_id="MichaelisConstant",
        table_name="Michaelis-Menten constant",
        table_type="QuantityMatrix",
        unit="mM",
    )
    sbtabdoc.add_sbtab(Km_sbtab)

    Ka_df = pd.DataFrame(
        Ka.m_as("mM").T,
        columns=metabolite_names
    )
    Ka_df.insert(0, "!QuantityType", "activation constant")
    Ka_df.insert(1, "!Reaction", reaction_names)
    Ka_sbtab = SBtab.SBtabTable.from_data_frame(
        Ka_df.astype(str),
        table_id="ActivationConstant",
        table_name="Activation constant",
        table_type="QuantityMatrix",
        unit="mM",
    )
    sbtabdoc.add_sbtab(Ka_sbtab)

    Ki_df = pd.DataFrame(
        Ki.m_as("mM").T,
        columns=metabolite_names
    )
    Ki_df.insert(0, "!QuantityType", "inhibition constant")
    Ki_df.insert(1, "!Reaction", reaction_names)
    Ki_sbtab = SBtab.SBtabTable.from_data_frame(
        Ki_df.astype(str),
        table_id="InhibitionConstant",
        table_name="Inhibition constant",
        table_type="QuantityMatrix",
        unit="mM",
    )
    sbtabdoc.add_sbtab(Ki_sbtab)
    return sbtabdoc

#%%
os.chdir("/home/eladn/git/model-balancing/python")
config_fname = "three_chain_model_artificial_noisy_state_noisy_kinetic"

with open(f"../cvxpy/examples/JSON/{config_fname}.json", "rt") as fp:
    data = json.load(fp)

keq_standard_concentration = Q_(data["standard_concentration"])

args = {}
args["S"] = np.array(data["network"]["stoichiometric_matrix"])
args["A_act"] = np.array(data["network"]["activation_matrix"])
args["A_inh"] = np.array(data["network"]["inhibition_matrix"])
args["fluxes"] = Q_(data["reaction_fluxes"]["data"]["mean"], data["reaction_fluxes"]["unit"])

args["Keq_gmean"] = Q_(data["kinetic_constants"]["Keq"]["combined"]["geom_mean"], "")

# "fix" the standard concentration (to a convention where it is 1 M which is used internally)
args["Keq_gmean"] *= np.power(keq_standard_concentration.m_as("M"), args["S"].sum(axis=0))

args["Keq_ln_cov"] = np.array(data["kinetic_constants"]["Keq"]["combined"]["cov_ln"])

args["kcatf_gmean"] = Q_(data["kinetic_constants"]["Kcatf"]["combined"]["geom_mean"], data["kinetic_constants"]["Kcatf"]["unit"])
args["kcatf_ln_cov"] = np.array(data["kinetic_constants"]["Kcatf"]["combined"]["cov_ln"])

args["kcatr_gmean"] = Q_(data["kinetic_constants"]["Kcatr"]["combined"]["geom_mean"], data["kinetic_constants"]["Kcatr"]["unit"])
args["kcatr_ln_cov"] = np.array(data["kinetic_constants"]["Kcatr"]["combined"]["cov_ln"])

args["Km_gmean"] = Q_(data["kinetic_constants"]["KM"]["combined"]["geom_mean"], "mM") #data["kinetic_constants"]["KM"]["unit"])
args["Km_ln_cov"] = np.array(data["kinetic_constants"]["KM"]["combined"]["cov_ln"])

args["Ka_gmean"] = Q_(data["kinetic_constants"]["KA"]["combined"]["geom_mean"], "mM") #data["kinetic_constants"]["KA"]["unit"])
args["Ka_ln_cov"] = np.array(data["kinetic_constants"]["KA"]["combined"]["cov_ln"])

args["Ki_gmean"] = Q_(data["kinetic_constants"]["KI"]["combined"]["geom_mean"], "mM") #data["kinetic_constants"]["KI"]["unit"])
args["Ki_ln_cov"] = np.array(data["kinetic_constants"]["KI"]["combined"]["cov_ln"])

args["conc_met_gmean"] = Q_(data["metabolite_concentrations"]["combined"]["geom_mean"], data["metabolite_concentrations"]["unit"])
args["conc_met_gstd"] = np.array(data["metabolite_concentrations"]["combined"]["geom_std"])

args["conc_enz_gmean"] = Q_(data["enzyme_concentrations"]["combined"]["geom_mean"], data["metabolite_concentrations"]["unit"])
args["conc_enz_gstd"] = np.array(data["enzyme_concentrations"]["combined"]["geom_std"])

args["rate_law"] = "CM"

metabolite_names = data["network"]["metabolite_names"]
reaction_names = data["network"]["reaction_names"]
state_names = data["state_names"]

#%%
mbc = ModelBalancingConvex(**args)
assert mbc.is_gmean_feasible()
mbc.solve()
print(f"Convex optimization (equivalent to α = 0) ... optimized total squared Z-scores = {mbc.objective_value:.3f}")
to_sbtab(mbc, metabolite_names, reaction_names, state_names).write(f"../res/{config_fname}_convex.tsv")

#%%
mb = ModelBalancing(**args)
assert mb.is_gmean_feasible()
for a in [0., 0.001, 0.01, 0.1, 0.5, 1.0]:
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    print(f"Solving using non-convex solver, α = {a:5.1g} ... ", end="")
    mb.alpha = a
    mb.solve()
    print(f"optimized total squared Z-scores = {mb.objective_value:.3f}")
    to_sbtab(mb, metabolite_names, reaction_names, state_names).write(f"../res/{config_fname}_alpha_{a:.1g}.tsv")