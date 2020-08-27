from typing import Union

import itertools
import cvxpy as cp
import numpy as np
from cvxpy.constraints.constraint import Constraint
from cvxpy.expressions.expression import Expression
from cvxpy.problems.objective import Objective

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

def _create_dense_matrix(
    S: np.array,
    x: Union[np.array, cp.Variable],
) -> Expression:
    """Converts a sparse list of affinity parameters (e.g. Km) to a matrix."""
    
    Nc, Nr = S.shape
    
    if x.size == 0:
        return np.zeros((Nc, Nr))
    
    K_mat = []
    k = 0
    for i in range(Nc):
        row = []
        for j in range(Nr):
            if S[i, j] != 0:
                row.append(x[k])
                k += 1
            else:
                row.append(cp.Constant(0))
        K_mat.append(cp.hstack(row))
    K_mat = cp.vstack(K_mat)
    return K_mat

def _z_score(
    x: Union[np.array, cp.Expression],
    mu: np.array,
    cov: np.array,
) -> Expression:
    """Calculates the sum of squared Z-scores (with a covariance mat)."""
    return cp.quad_form(x - mu, np.linalg.pinv(cov))

def _driving_forces(
    ln_Keq: Union[np.array, cp.Expression],
    ln_conc_met: Union[np.array, cp.Expression],
    S: np.array,
) -> Expression:
    """Calculates the driving forces of all reactions."""
    Ncond = ln_conc_met.shape[1]
    return cp.vstack([ln_Keq] * Ncond).T - S.T @ ln_conc_met

def _ln_kcatr(
    ln_kcatf: Union[np.array, cp.Expression],
    ln_Km: Union[np.array, cp.Expression],
    ln_Keq: Union[np.array, cp.Expression],
    S: np.array,
) -> Expression:
    """Calculate the kcat-reverse based on Haldane relationship constraint."""
    ln_Km_matrix = _create_dense_matrix(S, ln_Km)
    return cp.diag(S.T @ ln_Km_matrix) + ln_kcatf - ln_Keq

def _ln_capacity(
    fluxes: np.array,
    ln_kcatf: Union[np.array, cp.Expression],
) -> Expression:
    """Calculate the capacity term of the enzyme."""
    Ncond = fluxes.shape[1]
    return np.log(fluxes.m_as("M/s")) - cp.vstack([ln_kcatf] * Ncond).T

def _ln_eta_thermodynamic(
    driving_forces: Union[np.array, cp.Expression]
) -> Expression:
    """Calculate the thermodynamic term of the enzyme."""
    return cp.log(1.0 - cp.exp(-driving_forces))

def _ln_eta_kinetic(
    ln_conc_met: Union[np.array, cp.Expression],
    ln_Km: Union[np.array, cp.Expression],
    S: np.array,
    rate_law: str = "CM",
) -> Expression:
    """Calculate the kinetic (saturation) term of the enzyme."""
    Nc, Nr = S.shape
    Ncond = ln_conc_met.shape[1]
    
    S_subs = abs(S)
    S_prod = abs(S)
    S_subs[S > 0] = 0
    S_prod[S < 0] = 0
    
    ln_Km_matrix = _create_dense_matrix(S, ln_Km)

    ln_eta_kinetic = []
    for i in range(Ncond):
        ln_conc_met_matrix = cp.vstack([ln_conc_met[:, i]] * Nr).T
        ln_D_S = S_subs.T @ (ln_conc_met_matrix - ln_Km_matrix)
        ln_D_P = S_prod.T @ (ln_conc_met_matrix - ln_Km_matrix)
        
        if rate_law == "S":
            ln_eta_kinetic += [cp.Constant(np.zeros(Nr,))]
        elif rate_law == "1S":
            ln_eta_kinetic += [-cp.diag(cp.logistic(-ln_D_S))]
        elif rate_law == "SP":
            ln_eta_kinetic += [-cp.diag(cp.logistic(-ln_D_S + ln_D_P))]
        elif rate_law == "1SP":
            ln_eta_kinetic += [-cp.diag(cp.logistic(-ln_D_S + cp.logistic(ln_D_P)))]
        elif rate_law == "CM":
            _ln_eta_of_reaction = []
            for j in range(Nr):
                B = _B_matrix(Nc, S_subs[:, j], S_prod[:, j])
                _ln_eta_of_reaction += [-cp.log_sum_exp(B @ (ln_conc_met[:, i] - ln_Km_matrix[:, j]))]
            ln_eta_kinetic += [cp.reshape(cp.vstack(_ln_eta_of_reaction), (Nr,))]
        else:
            raise ValueError(f"unsupported rate law {rate_law}")

    return cp.vstack(ln_eta_kinetic).T

def _ln_eta_regulation(
    ln_conc_met: Union[np.array, cp.Expression],
    ln_Ka: Union[np.array, cp.Expression],
    ln_Ki: Union[np.array, cp.Expression],
    A_act: np.array,
    A_inh: np.array
) -> Expression:
    """Calculate the regulation (allosteric) term of the enzyme."""
    Nc, Nr = A_act.shape
    Ncond = ln_conc_met.shape[1]

    ln_Ka_matrix = _create_dense_matrix(A_act, ln_Ka)
    ln_Ki_matrix = _create_dense_matrix(A_inh, ln_Ki)
    
    ln_eta_allosteric = []
    for i in range(Ncond):
        ln_conc_met_matrix = cp.vstack([ln_conc_met[:, i]] * Nr).T
        ln_act = A_act.T @ cp.logistic(ln_Ka_matrix - ln_conc_met_matrix)
        ln_inh = A_inh.T @ cp.logistic(ln_conc_met_matrix - ln_Ki_matrix)
        ln_eta_act = -cp.diag(ln_act)
        ln_eta_inh = -cp.diag(ln_inh)
        ln_eta_allosteric += [ln_eta_act + ln_eta_inh]
    return cp.vstack(ln_eta_allosteric).T

def _ln_conc_enz(
    ln_Keq: Union[np.array, cp.Expression],
    ln_kcatf: Union[np.array, cp.Expression],
    ln_Km: Union[np.array, cp.Expression],
    ln_Ka: Union[np.array, cp.Expression],
    ln_Ki: Union[np.array, cp.Expression],
    ln_conc_met: Union[np.array, cp.Expression],
    fluxes: np.array,
    S: np.array,
    A_act: np.array,
    A_inh: np.array,
    rate_law: str = "CM",
) -> Expression:
    """Calculate the required enzyme levels based on fluxes and rate laws."""
    driving_forces = _driving_forces(ln_Keq, ln_conc_met, S)
    ln_capacity = _ln_capacity(fluxes, ln_kcatf)
    ln_eta_thermodynamic = _ln_eta_thermodynamic(driving_forces)
    ln_eta_kinetic = _ln_eta_kinetic(ln_conc_met, ln_Km, S, rate_law=rate_law)
    ln_eta_regulation = _ln_eta_regulation(ln_conc_met, ln_Ka, ln_Ki, A_act, A_inh)
    return ln_capacity - ln_eta_thermodynamic - ln_eta_kinetic - ln_eta_regulation

def solve(
    S: np.array,
    fluxes: np.array,
    A_act: np.array,
    A_inh: np.array,
    Keq_gmean: np.array,
    Keq_ln_cov: np.array,
    conc_enz_gmean: np.array,
    conc_enz_gstd: np.array,
    conc_met_gmean: np.array,
    conc_met_gstd: np.array,
    kcatf_gmean: np.array,
    kcatf_ln_cov: np.array,
    kcatr_gmean: np.array,
    kcatr_ln_cov: np.array,
    Km_gmean: np.array,
    Km_ln_cov: np.array,
    Ka_gmean: np.array,
    Ka_ln_cov: np.array,
    Ki_gmean: np.array,
    Ki_ln_cov: np.array,
    rate_law: str = "CM",
    solver: str = "SCS",
):

    Nc, Nr = S.shape
    assert fluxes.shape[0] == Nr
    Ncond = fluxes.shape[1]
    assert A_act.shape == (Nc, Nr)
    assert A_inh.shape == (Nc, Nr)
    assert Keq_gmean.shape == (Nr, )
    assert Keq_ln_cov.shape == (Nr, Nr)
    assert conc_enz_gmean.shape == (Nr, Ncond)
    assert conc_enz_gstd.shape == (Nr, Ncond)
    assert conc_met_gmean.shape == (Nc, Ncond)
    assert conc_met_gstd.shape == (Nc, Ncond)

    # define the independent variables and their z-scores
    ln_conc_met = cp.Variable(shape=(Nc, Ncond), value=np.log(conc_met_gmean.m_as("M")))

    # conc_met is given as a matrix (with conditions as columns) and therefore
    # we assume a diagonal covariance matrix (for simplicity). Instead of a
    # ln_cov matrix, we simply have the geometric means and stds arranged in 
    # the same shape as the variables.
    z2_scores_met = sum(cp.square(
        (ln_conc_met - np.log(conc_met_gmean.m_as("M"))) / np.log(conc_met_gstd)
    ).flatten())

    ln_Keq = cp.Variable(shape=Nr, value=np.log(Keq_gmean.m_as("")))
    z2_scores_Keq = _z_score(ln_Keq, np.log(Keq_gmean.m_as("")), Keq_ln_cov)

    ln_kcatf = cp.Variable(shape=(Nr), value=np.log(kcatf_gmean.m_as("1/s")))
    z2_scores_kcatf = _z_score(ln_kcatf, np.log(kcatf_gmean.m_as("1/s")), kcatf_ln_cov)

    if Km_gmean.size == 0:
        ln_Km = np.array([])
        z2_scores_Km = cp.Constant(0)
    else:
        ln_Km = cp.Variable(shape=Km_gmean.size, value=np.log(Km_gmean.m_as("M")))
        z2_scores_Km = _z_score(ln_Km, np.log(Km_gmean.m_as("M")), Km_ln_cov)

    if Ka_gmean.size == 0:
        ln_Ka = np.array([])
        z2_scores_Ka = cp.Constant(0)
    else:
        ln_Ka = cp.Variable(shape=Ka_gmean.size, value=np.log(Ka_gmean.m_as("M")))
        z2_scores_Ka = _z_score(ln_Ka, np.log(Ka_gmean.m_as("M")), Ka_ln_cov)

    if Ki_gmean.size == 0:
        ln_Ki = np.array([])
        z2_scores_Ki = cp.Constant(0)
    else:
        ln_Ki = cp.Variable(shape=Ki_gmean.size, value=np.log(Ki_gmean.m_as("M")))
        z2_scores_Ki = _z_score(ln_Ki, np.log(Ki_gmean.m_as("M")), Ki_ln_cov)

    # the dependent parameters are:
    # - reverse kcats
    # - enzyme concentrations
    
    ln_kcatr = _ln_kcatr(ln_kcatf, ln_Km, ln_Keq, S)
    z2_scores_kcatr = _z_score(ln_kcatr, np.log(kcatr_gmean.m_as("1/s")), kcatr_ln_cov)

    driving_forces = _driving_forces(ln_Keq, ln_conc_met, S)
    
    ### capacity term
    ln_capacity = _ln_capacity(fluxes, ln_kcatf)

    ### thermodynamic term
    ln_eta_thermodynamic = _ln_eta_thermodynamic(driving_forces)

    ### kinetic term
    ln_eta_kinetic = _ln_eta_kinetic(ln_conc_met, ln_Km, S, rate_law=rate_law)

    ### allosteric term
    ln_eta_regulation = _ln_eta_regulation(ln_conc_met, ln_Ka, ln_Ki, A_act, A_inh)
    
    ln_conc_enz = ln_capacity - ln_eta_thermodynamic - ln_eta_kinetic - ln_eta_regulation
        
    # ln enzyme concentrations are convex functions of the ln metabolite concentrations
    # but, since the z-scores use the square function, we have to take only the positive
    # values (otherwise the result is not convex).
    z2_scores_enz = sum(cp.square(
        cp.pos(ln_conc_enz - np.log(conc_enz_gmean.m_as("M"))) / np.log(conc_enz_gstd)
    ).flatten())
        
    objective = cp.Minimize(
        z2_scores_met +
        z2_scores_Keq +
        z2_scores_kcatf +
        z2_scores_Km +
        z2_scores_Ka +
        z2_scores_Ki +
        z2_scores_kcatr +
        z2_scores_enz
    )
    
    print(f"initial Z-score = {objective.value:.2e}")

    prob = cp.Problem(objective, [driving_forces >= 1e-6])
    prob.solve(solver=solver)

    print(f"optimal Z-score = {objective.value:.2e}")

    print("\nMetabolite concentrations (M) =\n", np.exp(ln_conc_met.value))
    print("\nEnzyme concentrations (M) =\n", np.exp(ln_conc_enz.value))
    print("\nDriving forces (RT) =\n", driving_forces.value)
    print("\nη(thr) =\n", np.exp(ln_eta_thermodynamic.value).round(2))
    print("\nη(kin) =\n", np.exp(ln_eta_kinetic.value).round(2))
    print("\nη(reg) =\n", np.exp(ln_eta_regulation.value).round(2))
    print("\n\n\n")
    print("All Z-scores\n" + "-"*50)
    print("enzymes = ", z2_scores_enz.value.round(2))
    print("metabolites = ", z2_scores_met.value.round(2))
    print("Keq = ", z2_scores_Keq.value.round(2))
    print("kcat forward = ", z2_scores_kcatf.value.round(2))
    print("kcat reverse = ", z2_scores_kcatr.value.round(2))
    print("Km = ", z2_scores_Km.value.round(2))
    print("Ka = ", z2_scores_Ka.value.round(2))
    print("Ki = ", z2_scores_Ki.value.round(2))