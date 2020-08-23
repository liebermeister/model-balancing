import itertools
import cvxpy as cp
import numpy as np
import pandas as pd
from cvxpy.constraints.constraint import Constraint
from cvxpy.expressions.expression import Expression
from cvxpy.problems.objective import Objective
from equilibrator_api import Q_, R, default_T

# helper function for the CM rate law

def _B_matrix(
    Nc: int, col_subs: np.ndarray, col_prod: np.ndarray
) -> np.ndarray:
    """Build the B matrix for the eta^kin expression.

    row_subs : np.ndarray
        A column from the substrate stoichiometric matrix. We assume
        coefficients represent reactant molecularities so
        only integer values are allowed.

    row_prod : np.ndarray
        A column from the product stoichiometric matrix. We assume
        coefficients represent reactant molecularities so
        only integer values are allowed.
    """

    def K_matrix(n: int) -> np.ndarray:
        """Make the 'K' matrix for the CM rate law."""
        lst = list(itertools.product([0, 1], repeat=n))
        lst.pop(0)  # remove the [0, 0, ..., 0] row
        return np.array(lst)

    def expand_S(coeffs: np.ndarray) -> np.ndarray:
        """Expand a coefficient column into a matrix with duplicates."""
        cs = list(np.cumsum(list(map(int, coeffs.flat))))
        S_tmp = np.zeros((cs[-1], Nc))
        for j, (i_from, i_to) in enumerate(zip([0] + cs, cs)):
            S_tmp[i_from:i_to, j] = 1
        return S_tmp

    S_subs = expand_S(col_subs)
    S_prod = expand_S(col_prod)

    A = np.vstack(
        [
            np.zeros((1, Nc)),
            K_matrix(S_subs.shape[0]) @ S_subs,
            K_matrix(S_prod.shape[0]) @ S_prod,
        ]
    )

    return A - np.ones((A.shape[0], S_subs.shape[0])) @ S_subs

def create_sparse_variable_matrix(S, mu: np.array, ln_cov: np.array) -> Expression:
    Nc, Nr = S.shape

    assert mu.ndim == 1
    ln_K = cp.Variable(mu.size)
    z2_scores = cp.quad_form(ln_K - np.log(mu), np.linalg.pinv(ln_cov))
    
    K_mat = []
    k = 0
    for i in range(Nc):
        row = []
        for j in range(Nr):
            if S[i, j] != 0:
                row.append(ln_K[k])
                k += 1
            else:
                row.append(cp.Constant(0))
        K_mat.append(cp.hstack(row))
    K_mat = cp.vstack(K_mat)
    return K_mat, z2_scores

def solve(
    S,
    fluxes,
    A_act,
    A_inh,
    Keq_gmean, Keq_ln_cov,
    conc_enz_gmean, conc_enz_gsigma,
    conc_met_gmean, conc_met_gsigma,
    kcatf_gmean, kcatf_ln_cov,
    kcatr_gmean, kcatr_ln_cov,
    Km_gmean, Km_ln_cov,
    Ka_gmean, Ka_ln_cov,
    Ki_gmean, Ki_ln_cov,
    rate_law: str = "CM"
):

    Nc, Nr = S.shape
    assert fluxes.shape[0] == Nr
    Ncond = fluxes.shape[1]
    assert A_act.shape == (Nc, Nr)
    assert A_inh.shape == (Nc, Nr)
    assert Keq_gmean.shape == (Nr, )
    assert Keq_ln_cov.shape == (Nr, Nr)
    assert conc_enz_gmean.shape == (Nr, Ncond)
    assert conc_enz_gsigma.shape == (Nr, Ncond)
    assert conc_met_gmean.shape == (Nc, Ncond)
    assert conc_met_gsigma.shape == (Nc, Ncond)

    # helper matrices
    S_subs = abs(S)
    S_prod = abs(S)
    S_subs[S > 0] = 0
    S_prod[S < 0] = 0

    # define the independent variables
    standard_dg = cp.Variable(shape=Nr)
    ln_kcatf = cp.Variable(shape=(Nr))  # in 1/s
    ln_Km, z2_scores_Km = create_sparse_variable_matrix(S, Km_gmean.m_as("M"), Km_ln_cov)
    ln_Ka, z2_scores_Ka = create_sparse_variable_matrix(A_act, Ka_gmean.m_as("M"), Ka_ln_cov)
    ln_Ki, z2_scores_Ki = create_sparse_variable_matrix(A_inh, Ki_gmean.m_as("M"), Ki_ln_cov)
    
    # define the dependencies
    standard_dg_over_rt = (standard_dg / (R * default_T).m_as("kJ/mol"))

    # Haldane relationship constraint
    # diag(S.T @ ln_KMM) + ln_kcatf - ln_kcatr = ln(Keq) = -ΔG'0 / RT
    ln_kcatr = cp.diag(S.T @ ln_Km) + ln_kcatf + standard_dg_over_rt

    # calculate Z-scores for deltaG0, kcatf and kcatr
    z2_scores_dg = cp.quad_form(
        -standard_dg_over_rt - np.log(Keq_gmean),
        np.linalg.pinv(Keq_ln_cov)
    )

    z2_scores_kcatf = cp.quad_form(
        ln_kcatf - np.log(kcatf_gmean.m_as("1/s")),
        np.linalg.pinv(kcatf_ln_cov)
    )

    z2_scores_kcatr = cp.quad_form(
        ln_kcatr - np.log(kcatr_gmean.m_as("1/s")),
        np.linalg.pinv(kcatr_ln_cov)
    )

    ln_conc_met = cp.Variable(shape=(Nc, Ncond))

    driving_forces = -cp.vstack([standard_dg_over_rt] * Ncond).T - S.T @ ln_conc_met
    
    ### capacity term
    ln_capacity = np.log(fluxes.m_as("M/s")) - cp.vstack([ln_kcatf] * Ncond).T

    ### thermodynamic term
    ln_eta_thermodynamic = cp.log(1.0 - cp.exp(-driving_forces))

    ### kinetic term
    ln_eta_kinetic = []
    for i in range(Ncond):
        ln_conc_met_matrix = cp.vstack([ln_conc_met[:, i]] * Nr).T
        ln_D_S = S_subs.T @ (ln_conc_met_matrix - ln_Km)
        ln_D_P = S_prod.T @ (ln_conc_met_matrix - ln_Km)
        
        if rate_law == "S":
            ln_eta_kinetic += [cp.Constant(np.zeros(self.Nr))]
        elif rate_law == "1S":
            ln_eta_kinetic += [-cp.diag(cp.logistic(-ln_D_S))]
        elif rate_law == "SP":
            ln_eta_kinetic += [-cp.diag(cp.logistic(-ln_D_S + ln_D_P))]
        elif rate_law == "1SP":
            ln_eta_kinetic += [-cp.diag(cp.logistic(-ln_D_S + cp.logistic(ln_D_P)))]
        elif rate_law == "CM":
            ln_eta = []
            for j in range(Nr):
                B = _B_matrix(Nc, S_subs[:, j], S_prod[:, j])
                ln_eta += [-cp.log_sum_exp(B @ (ln_conc_met[:, i] - ln_Km[:, j]))]
            ln_eta_kinetic += [cp.reshape(cp.vstack(ln_eta), (Nr,))]
        else:
            raise ValueError(f"unsupported rate law {rate_law}")

    ln_eta_kinetic = cp.vstack(ln_eta_kinetic).T

    ### allosteric term
    ln_eta_allosteric = []
    for i in range(Ncond):
        ln_conc_met_matrix = cp.vstack([ln_conc_met[:, i]] * Nr).T
        ln_act = A_act.T @ cp.logistic(ln_Ka - ln_conc_met_matrix)
        ln_inh = A_inh.T @ cp.logistic(ln_conc_met_matrix - ln_Ki)
        ln_eta_act = -cp.diag(ln_act)
        ln_eta_inh = -cp.diag(ln_inh)
        ln_eta_allosteric += [ln_eta_act + ln_eta_inh]
    ln_eta_allosteric = cp.vstack(ln_eta_allosteric).T
    
    ln_conc_enz = ln_capacity - ln_eta_thermodynamic - ln_eta_kinetic - ln_eta_allosteric
        
    # ln enzyme concentrations are convex functions of the ln metabolite concentrations
    # but, since the z-scores use the square function, we have to take only the positive
    # values (otherwise the result is not convex).
    z2_scores_enz = sum(cp.square(
        cp.pos(ln_conc_enz - np.log(conc_enz_gmean.m_as("M"))) / np.log(conc_enz_gsigma)
    ).flatten())

    z2_scores_met = sum(cp.square(
        (ln_conc_met - np.log(conc_met_gmean.m_as("M"))) / np.log(conc_met_gsigma)
    ).flatten())
        
    objective = cp.Minimize(
        z2_scores_enz +
        z2_scores_met +
        z2_scores_dg +
        z2_scores_kcatf +
        z2_scores_kcatr +
        z2_scores_Km +
        z2_scores_Ka +
        z2_scores_Ki
    )

    prob = cp.Problem(objective, [driving_forces >= 1e-6])
    prob.solve()

    print("\nMetabolite concentrations =\n", np.exp(ln_conc_met.value))
    print("\nEnzyme concentrations =\n", np.exp(ln_conc_enz.value))
    print("\neta_thr =\n", np.exp(ln_eta_thermodynamic.value).round(2))
    print("\neta_kin =\n", np.exp(ln_eta_kinetic.value).round(2))
    print("\neta_reg =\n", np.exp(ln_eta_allosteric.value).round(2))
    print("\n\n\n")
    print("All Z-scores\n" + "-"*50)
    print("enzymes = ", z2_scores_enz.value.round(2))
    print("metabolites = ", z2_scores_met.value.round(2))
    print("ΔG'0 = ", z2_scores_dg.value.round(2))
    print("kcat forward = ", z2_scores_kcatf.value.round(2))
    print("kcat reverse = ", z2_scores_kcatr.value.round(2))
    print("Km = ", z2_scores_Km.value.round(2))
    print("Ka = ", z2_scores_Ka.value.round(2))
    print("Ki = ", z2_scores_Ki.value.round(2))