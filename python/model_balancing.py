from typing import Union
import os
import warnings
import itertools
import pint
import cvxpy as cp
import numpy as np
from cvxpy.constraints.constraint import Constraint
from cvxpy.expressions.expression import Expression
from cvxpy.problems.objective import Objective

# Disable Pint's old fallback behavior (must come before importing Pint)
os.environ["PINT_ARRAY_PROTOCOL_FALLBACK"] = "0"

ureg = pint.UnitRegistry(system="mks")
Q_ = ureg.Quantity

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Q_([])

class ModelBalancing(object):
    
    MIN_DRIVING_FORCE = 1e-6  # in units of RT

    def __init__(
        self,
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
        solver: str = "MOSEK",
    ) -> None:
        self.S = S
        self.fluxes = fluxes
        self.A_act = A_act
        self.A_inh = A_inh

        self.Nc, self.Nr = S.shape
        assert self.fluxes.shape[0] == self.Nr
        self.Ncond = self.fluxes.shape[1]
        assert self.A_act.shape == (self.Nc, self.Nr)
        assert self.A_inh.shape == (self.Nc, self.Nr)
        
        self.ln_Keq_gmean = cp.Parameter(shape=(self.Nr, ), value=np.log(Keq_gmean.m_as("")))
        self.Keq_ln_cov = Keq_ln_cov
        self.ln_conc_enz_gmean = cp.Parameter(shape=(self.Nr, self.Ncond), value=np.log(conc_enz_gmean.m_as("M")))
        self.ln_conc_enz_gstd = np.log(conc_enz_gstd)
        self.ln_conc_met_gmean = cp.Parameter(shape=(self.Nc, self.Ncond), value=np.log(conc_met_gmean.m_as("M")))
        self.ln_conc_met_gstd = np.log(conc_met_gstd)
        self.ln_kcatf_gmean = cp.Parameter(shape=(self.Nr, ), value=np.log(kcatf_gmean.m_as("1/s")))
        self.kcatf_ln_cov = kcatf_ln_cov
        self.ln_kcatr_gmean = cp.Parameter(shape=(self.Nr, ), value=np.log(kcatr_gmean.m_as("1/s")))
        self.kcatr_ln_cov = kcatr_ln_cov
        self.ln_Km_gmean = np.log(Km_gmean.m_as("M"))
        self.Km_ln_cov = Km_ln_cov
        self.ln_Ka_gmean = np.log(Ka_gmean.m_as("M"))
        self.Ka_ln_cov = Ka_ln_cov
        self.ln_Ki_gmean = np.log(Ki_gmean.m_as("M"))
        self.Ki_ln_cov = Ki_ln_cov
        self.rate_law = "CM"
        self.solver = "MOSEK"

        assert self.Keq_ln_cov.shape == (self.Nr, self.Nr)
        assert self.ln_conc_enz_gstd.shape == (self.Nr, self.Ncond)
        assert self.ln_conc_met_gstd.shape == (self.Nc, self.Ncond)
        
        self.ln_conc_met = cp.Variable(shape=(self.Nc, self.Ncond))
        self.ln_Keq = cp.Variable(shape=(self.Nr, ))
        self.ln_kcatf = cp.Variable(shape=(self.Nr, ))
        
        for p in ["Km", "Ka", "Ki"]:
            ln_gmean = self.__getattribute__(f"ln_{p}_gmean")
            ln_cov = self.__getattribute__(f"{p}_ln_cov")
            if ln_gmean.size == 0:
                self.__setattr__(f"ln_{p}", np.array([]))
                self.__setattr__(f"z2_scores_{p}", cp.Constant(0))
            else:
                ln_p = cp.Variable(shape=ln_gmean.size)
                self.__setattr__(f"z2_scores_{p}", ModelBalancing._z_score(ln_p, ln_gmean, ln_cov))
                self.__setattr__(f"ln_{p}", ln_p)
        
        for p in ["Keq", "kcatf", "kcatr"]:
            ln_p = self.__getattribute__(f"ln_{p}")
            ln_gmean = self.__getattribute__(f"ln_{p}_gmean")
            ln_cov = self.__getattribute__(f"{p}_ln_cov")
            self.__setattr__(f"z2_scores_{p}", ModelBalancing._z_score(ln_p, ln_gmean, ln_cov))

        # conc_met is given as a matrix (with conditions as columns) and therefore
        # we assume a diagonal covariance matrix (for simplicity). Instead of a
        # ln_cov matrix, we simply have the geometric means and stds arranged in 
        # the same shape as the variables.
        self.z2_scores_met = sum(cp.square(
            (self.ln_conc_met - self.ln_conc_met_gmean) / self.ln_conc_met_gstd
        ).flatten())

        # ln enzyme concentrations are convex functions of the ln metabolite concentrations
        # but, since the z-scores use the square function, we have to take only the positive
        # values (otherwise the result is not convex).
        self.z2_scores_enz = sum(cp.square(
            cp.pos(self.ln_conc_enz - self.ln_conc_enz_gmean) / self.ln_conc_enz_gstd
        ).flatten())

        total_z2_scores = sum(
            [
                self.__getattribute__(f"z2_scores_{p}")
                for p in ["Km", "Ka", "Ki", "Keq", "kcatf", "kcatr", "met", "enz"]
            ]
        )
        self.objective = cp.Minimize(total_z2_scores)

    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def _z_score(
        x: Union[np.array, cp.Expression],
        mu: np.array,
        cov: np.array,
    ) -> Expression:
        """Calculates the sum of squared Z-scores (with a covariance mat)."""
        return cp.quad_form(x - mu, np.linalg.pinv(cov))

    def _driving_forces(
        self,
        ln_Keq: Union[np.array, cp.Expression],
        ln_conc_met: Union[np.array, cp.Expression],
    ) -> Expression:
        """Calculates the driving forces of all reactions."""
        return cp.vstack([ln_Keq] * self.Ncond).T - self.S.T @ ln_conc_met

    @property
    def driving_forces(self) -> Expression:
        return self._driving_forces(self.ln_Keq, self.ln_conc_met)
    
    @staticmethod
    def _ln_kcatr(
        S: np.array,
        ln_kcatf: Union[np.array, cp.Expression],
        ln_Km: Union[np.array, cp.Expression],
        ln_Keq: Union[np.array, cp.Expression],
    ) -> Expression:
        """Calculate the kcat-reverse based on Haldane relationship constraint."""
        ln_Km_matrix = ModelBalancing._create_dense_matrix(S, ln_Km)
        return cp.diag(S.T @ ln_Km_matrix) + ln_kcatf - ln_Keq

    @property
    def ln_kcatr(self) -> Expression:
        return ModelBalancing._ln_kcatr(self.S, self.ln_kcatf, self.ln_Km, self.ln_Keq)
    
    def _ln_capacity(
        self,
        ln_kcatf: Union[np.array, cp.Expression],
    ) -> Expression:
        """Calculate the capacity term of the enzyme."""
        return np.log(self.fluxes.m_as("M/s")) - cp.vstack([ln_kcatf] * self.Ncond).T

    def _ln_eta_thermodynamic(
        self,
        driving_forces: Union[np.array, cp.Expression]
    ) -> Expression:
        """Calculate the thermodynamic term of the enzyme."""
        return cp.log(1.0 - cp.exp(-driving_forces))

    @property
    def ln_eta_thermodynamic(self) -> Expression:
        return self._ln_eta_thermodynamic(self.driving_forces)
    
    def _ln_eta_kinetic(
        self,
        ln_conc_met: Union[np.array, cp.Expression],
        ln_Km: Union[np.array, cp.Expression],
    ) -> Expression:
        """Calculate the kinetic (saturation) term of the enzyme."""
        S_subs = abs(self.S)
        S_prod = abs(self.S)
        S_subs[self.S > 0] = 0
        S_prod[self.S < 0] = 0

        ln_Km_matrix = ModelBalancing._create_dense_matrix(self.S, ln_Km)

        ln_eta_kinetic = []
        for i in range(self.Ncond):
            ln_conc_met_matrix = cp.vstack([ln_conc_met[:, i]] * self.Nr).T
            ln_D_S = S_subs.T @ (ln_conc_met_matrix - ln_Km_matrix)
            ln_D_P = S_prod.T @ (ln_conc_met_matrix - ln_Km_matrix)

            if self.rate_law == "S":
                ln_eta_kinetic += [cp.Constant(np.zeros(self.Nr,))]
            elif self.rate_law == "1S":
                ln_eta_kinetic += [-cp.diag(cp.logistic(-ln_D_S))]
            elif self.rate_law == "SP":
                ln_eta_kinetic += [-cp.diag(cp.logistic(-ln_D_S + ln_D_P))]
            elif self.rate_law == "1SP":
                ln_eta_kinetic += [-cp.diag(cp.logistic(-ln_D_S + cp.logistic(ln_D_P)))]
            elif self.rate_law == "CM":
                ln_eta_of_reaction = []
                for j in range(self.Nr):
                    B = ModelBalancing._B_matrix(self.Nc, S_subs[:, j], S_prod[:, j])
                    ln_eta_of_reaction += [-cp.log_sum_exp(B @ (ln_conc_met[:, i] - ln_Km_matrix[:, j]))]
                ln_eta_kinetic += [cp.reshape(cp.vstack(ln_eta_of_reaction), (self.Nr,))]
            else:
                raise ValueError(f"unsupported rate law {self.rate_law}")

        return cp.vstack(ln_eta_kinetic).T

    @property
    def ln_eta_kinetic(self) -> Expression:
        return self._ln_eta_kinetic(self.ln_conc_met, self.ln_Km)
    
    def _ln_eta_regulation(
        self,
        ln_conc_met: Union[np.array, cp.Expression],
        ln_Ka: Union[np.array, cp.Expression],
        ln_Ki: Union[np.array, cp.Expression],
    ) -> Expression:
        """Calculate the regulation (allosteric) term of the enzyme."""
        ln_Ka_matrix = ModelBalancing._create_dense_matrix(self.A_act, ln_Ka)
        ln_Ki_matrix = ModelBalancing._create_dense_matrix(self.A_inh, ln_Ki)

        ln_eta_allosteric = []
        for i in range(self.Ncond):
            ln_conc_met_matrix = cp.vstack([ln_conc_met[:, i]] * self.Nr).T
            ln_act = self.A_act.T @ cp.logistic(ln_Ka_matrix - ln_conc_met_matrix)
            ln_inh = self.A_inh.T @ cp.logistic(ln_conc_met_matrix - ln_Ki_matrix)
            ln_eta_act = -cp.diag(ln_act)
            ln_eta_inh = -cp.diag(ln_inh)
            ln_eta_allosteric += [ln_eta_act + ln_eta_inh]
        return cp.vstack(ln_eta_allosteric).T

    @property
    def ln_eta_regulation(self) -> Expression:
        return self._ln_eta_regulation(self.ln_conc_met, self.ln_Ka, self.ln_Ki)
    
    def _ln_conc_enz(
        self,
        ln_Keq: Union[np.array, cp.Expression],
        ln_kcatf: Union[np.array, cp.Expression],
        ln_Km: Union[np.array, cp.Expression],
        ln_Ka: Union[np.array, cp.Expression],
        ln_Ki: Union[np.array, cp.Expression],
        ln_conc_met: Union[np.array, cp.Expression],
    ) -> Expression:
        """Calculate the required enzyme levels based on fluxes and rate laws."""
        driving_forces = self._driving_forces(ln_Keq, ln_conc_met)
        ln_capacity = self._ln_capacity(ln_kcatf)
        ln_eta_thermodynamic = self._ln_eta_thermodynamic(driving_forces)
        ln_eta_kinetic = self._ln_eta_kinetic(ln_conc_met, ln_Km)
        ln_eta_regulation = self._ln_eta_regulation(ln_conc_met, ln_Ka, ln_Ki)
        return ln_capacity - ln_eta_thermodynamic - ln_eta_kinetic - ln_eta_regulation

    @property
    def ln_conc_enz(self) -> Expression:
        return self._ln_conc_enz(self.ln_Keq, self.ln_kcatf, self.ln_Km, self.ln_Ka, self.ln_Ki, self.ln_conc_met)
        
    def is_gmean_feasible(self) -> bool:
        return (self._driving_forces(self.ln_Keq_gmean, self.ln_conc_met_gmean).value >= self.MIN_DRIVING_FORCE).all()
    
    def initialize_with_gmeans(self) -> None:
        # define the independent variables and their z-scores
        self.ln_conc_met.value = self.ln_conc_met_gmean.value
        self.ln_Keq.value = self.ln_Keq_gmean.value
        self.ln_kcatf.value = self.ln_kcatf_gmean.value

        if self.ln_Km_gmean.size != 0:
            self.ln_Km.value = self.ln_Km_gmean

        if self.ln_Ka_gmean.size != 0:
            self.ln_Ka.value = self.ln_Ka_gmean

        if self.ln_Ki_gmean.size != 0:
            self.ln_Ki.value = self.ln_Ki_gmean
           
        self.ln_kcatr_gmean.value = self.ln_kcatr.value
        self.ln_conc_enz_gmean.value = self.ln_conc_enz.value
    
    def solve(self) -> None:
        prob = cp.Problem(self.objective, [self.driving_forces >= self.MIN_DRIVING_FORCE])
        prob.solve(solver=self.solver)

    @property
    def objective_value(self) -> float:
        return self.objective.value

    def print_z_scores(self, percision: int = 2) -> None:
        for p in ["enz", "met", "Keq", "kcatf", "kcatr", "Km", "Ka", "Ki"]:
            z = self.__getattribute__(f"z2_scores_{p}").value
            print(f"{p} = {z.round(percision)}")
        
    def print_status(self) -> None:
        print("\nMetabolite concentrations (M) =\n", np.exp(self.ln_conc_met.value))
        print("\nEnzyme concentrations (M) =\n", np.exp(self.ln_conc_enz.value))
        print("\nDriving forces (RT) =\n", self.driving_forces.value)
        print("\nη(thr) =\n", np.exp(self.ln_eta_thermodynamic.value).round(2))
        print("\nη(kin) =\n", np.exp(self.ln_eta_kinetic.value).round(2))
        print("\nη(reg) =\n", np.exp(self.ln_eta_regulation.value).round(2))
        print("\n\n\n")
