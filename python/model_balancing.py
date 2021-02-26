import itertools
import os
import warnings
from scipy.optimize import minimize

import numpy as np
import scipy.special
import pint

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
        assert fluxes.shape[0] == S.shape[1]
        assert A_act.shape == S.shape
        assert A_inh.shape == S.shape
        assert Keq_ln_cov.shape == (S.shape[1], S.shape[1])
        assert conc_enz_gstd.shape == (S.shape[1], fluxes.shape[1])
        assert conc_met_gstd.shape == (S.shape[0], fluxes.shape[1])

        self.S = S
        self.fluxes = fluxes
        self.A_act = A_act
        self.A_inh = A_inh

        self.Nc, self.Nr = self.S.shape
        self.Ncond = self.fluxes.shape[1]

        self.conc_enz_ln_gmean = np.log(conc_enz_gmean.m_as("M").flatten())
        self.conc_enz_inv_ln_cov = np.diag(1.0 / np.log(conc_enz_gstd.flatten()))
        self.conc_met_ln_gmean = np.log(conc_met_gmean.m_as("M").flatten())
        self.conc_met_inv_ln_cov = np.diag(1.0 / np.log(conc_met_gstd.flatten()))

        self.Keq_ln_gmean = np.log(Keq_gmean.m_as(""))
        self.Keq_inv_ln_cov = np.linalg.pinv(Keq_ln_cov)
        self.kcatf_ln_gmean = np.log(kcatf_gmean.m_as("1/s"))
        self.kcatf_inv_ln_cov = np.linalg.pinv(kcatf_ln_cov)
        self.kcatr_ln_gmean = np.log(kcatr_gmean.m_as("1/s"))
        self.kcatr_inv_ln_cov = np.linalg.pinv(kcatr_ln_cov)
        self.Km_ln_gmean = np.log(Km_gmean.m_as("M"))
        self.Km_inv_ln_cov = np.linalg.pinv(Km_ln_cov)
        self.Ka_ln_gmean = np.log(Ka_gmean.m_as("M"))
        self.Ka_inv_ln_cov = np.linalg.pinv(Ka_ln_cov)
        self.Ki_ln_gmean = np.log(Ki_gmean.m_as("M"))
        self.Ki_inv_ln_cov = np.linalg.pinv(Ki_ln_cov)

        self.rate_law = "CM"
        self.solver = "MOSEK"

        self.ln_conc_met = np.zeros(self.Nc * self.Ncond)
        self.ln_Keq = np.zeros(self.Nr)
        self.ln_kcatf = np.zeros(self.Nr)
        self.ln_Ka = np.zeros(Ka_gmean.shape)
        self.ln_Ki = np.zeros(Ki_gmean.shape)
        self.ln_Km = np.zeros(Km_gmean.shape)

    def objective_function_old(
        self,
    ) -> float:
        """Calculate the sum of squares of all Z-scores.

        This includes the 6 independent variables, and the two dependent ones
        (kcatr and conc_enz).
        """
        total_z2_scores = 0.0

        for p in ["Km", "Ka", "Ki", "Keq", "kcatf", "kcatr", "conc_met", "conc_enz"]:
            p_ln_gmean = self.__getattribute__(f"{p}_ln_gmean")
            p_inv_ln_cov = self.__getattribute__(f"{p}_inv_ln_cov")
            ln_p = self.__getattribute__(f"ln_{p}")
            if p_ln_gmean.size == 0:
                continue
            else:
                total_z2_scores += ModelBalancing._z_score(ln_p, p_ln_gmean, p_inv_ln_cov)

        return total_z2_scores

    def objective_function(
        self,
        x: np.ndarray
    ) -> float:
        """Calculate the sum of squares of all Z-scores.

        The input (x) is a stacked version of all the independent variables, assuming
        the following order: Km, Ka, Ki, Keq, kcatf, conc_met
        """
        total_z2_scores = 0.0

        var_dict = {}
        i = 0
        for p in ["Km", "Ka", "Ki", "Keq", "kcatf", "conc_met"]:
            ln_p = self.__getattribute__(f"ln_{p}")
            var_dict[f"ln_{p}"] = x[i:i+ln_p.size].reshape(ln_p.shape)
            i += ln_p.size

        ln_kcatr = ModelBalancing._ln_kcatr(
            self.S, var_dict["ln_kcatf"], var_dict["ln_Km"], var_dict["ln_Keq"]
        )
        ln_conc_enz = self._ln_conc_enz(**var_dict)

        var_dict["ln_kcatr"] = ln_kcatr
        var_dict["ln_conc_enz"] = ln_conc_enz

        # first sum the Z-scores of the variables whose prior is a multivariate
        # log-normal distribution (with a full covariance matrix).
        # note that kcatr is a dependent parameter, while all the others are independent.
        for p in ["Km", "Ka", "Ki", "Keq", "kcatf", "kcatr", "conc_met", "conc_enz"]:
            ln_p_gmean = self.__getattribute__(f"{p}_ln_gmean")
            p_inv_ln_cov = self.__getattribute__(f"{p}_inv_ln_cov")
            ln_p = var_dict[f"ln_{p}"]
            if ln_p_gmean.size == 0:
                continue
            else:
                total_z2_scores += ModelBalancing._z_score(ln_p, ln_p_gmean, p_inv_ln_cov)

        return total_z2_scores

    @property
    def objective_value(self) -> float:
        x0 = []
        for p in ["Km", "Ka", "Ki", "Keq", "kcatf", "conc_met"]:
            x0.append(self.__getattribute__(f"ln_{p}").flatten())
        x0 = np.hstack(x0).flatten()
        return self.objective_function(x0)

    @staticmethod
    def _z_score(
        x: np.array,
        mu: np.array,
        inv_cov: np.array,
    ) -> float:
        """Calculates the sum of squared Z-scores (with a covariance mat)."""
        normed = (x.flatten() - mu.flatten()) @ inv_cov
        return sum(map(np.square, normed.flat))

    @staticmethod
    def _B_matrix(Nc: int, col_subs: np.ndarray, col_prod: np.ndarray) -> np.ndarray:
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
    def _logistic(x: np.ndarray) -> np.ndarray:
        """elementwise calculation of: log(1 + e ^ x)"""
        return np.log(1.0 + np.exp(x))

    @staticmethod
    def _create_dense_matrix(
        S: np.ndarray,
        x: np.ndarray,
    ) -> np.ndarray:
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
                    row.append(0)
            K_mat.append(np.hstack(row))
        K_mat = np.vstack(K_mat)
        return K_mat

    def _driving_forces(
        self,
        ln_Keq: np.ndarray,
        ln_conc_met: np.ndarray,
    ) -> np.ndarray:
        """Calculates the driving forces of all reactions."""
        return np.vstack([ln_Keq] * self.Ncond).T - self.S.T @ ln_conc_met.reshape((self.Nc, self.Ncond))

    @property
    def driving_forces(self) -> np.ndarray:
        return self._driving_forces(self.ln_Keq, self.ln_conc_met)

    @staticmethod
    def _ln_kcatr(
        S: np.array,
        ln_kcatf: np.ndarray,
        ln_Km: np.ndarray,
        ln_Keq: np.ndarray,
    ) -> np.ndarray:
        """Calculate the kcat-reverse based on Haldane relationship constraint."""
        ln_Km_matrix = ModelBalancing._create_dense_matrix(S, ln_Km)
        return np.diag(S.T @ ln_Km_matrix) + ln_kcatf - ln_Keq

    @property
    def ln_kcatr(self) -> np.ndarray:
        return ModelBalancing._ln_kcatr(self.S, self.ln_kcatf, self.ln_Km, self.ln_Keq)

    def _ln_capacity(
        self,
        ln_kcatf: np.ndarray,
    ) -> np.ndarray:
        """Calculate the capacity term of the enzyme."""
        return np.log(self.fluxes.m_as("M/s")) - np.vstack([ln_kcatf] * self.Ncond).T

    def _ln_eta_thermodynamic(self, driving_forces: np.ndarray) -> np.ndarray:
        """Calculate the thermodynamic term of the enzyme."""
        return np.log(1.0 - np.exp(-driving_forces))

    @property
    def ln_eta_thermodynamic(self) -> np.ndarray:
        return self._ln_eta_thermodynamic(self.driving_forces)

    def _ln_eta_kinetic(
        self,
        ln_conc_met: np.ndarray,
        ln_Km: np.ndarray,
    ) -> np.ndarray:
        """Calculate the kinetic (saturation) term of the enzyme."""
        S_subs = abs(self.S)
        S_prod = abs(self.S)
        S_subs[self.S > 0] = 0
        S_prod[self.S < 0] = 0

        ln_Km_matrix = ModelBalancing._create_dense_matrix(self.S, ln_Km)

        ln_eta_kinetic = []
        for i in range(self.Ncond):
            condition_ln_conc_met = ln_conc_met[self.Nc*i:self.Nc*(i+1)]
            ln_conc_met_matrix = np.vstack([condition_ln_conc_met] * self.Nr).T
            ln_D_S = S_subs.T @ (ln_conc_met_matrix - ln_Km_matrix)
            ln_D_P = S_prod.T @ (ln_conc_met_matrix - ln_Km_matrix)

            if self.rate_law == "S":
                ln_eta_kinetic += [
                    np.zeros(
                        self.Nr,
                    )
                ]
            elif self.rate_law == "1S":
                ln_eta_kinetic += [-np.diag(self._logistic(-ln_D_S))]
            elif self.rate_law == "SP":
                ln_eta_kinetic += [-np.diag(self._logistic(-ln_D_S + ln_D_P))]
            elif self.rate_law == "1SP":
                ln_eta_kinetic += [
                    -np.diag(self._logistic(-ln_D_S + self._logistic(ln_D_P)))
                ]
            elif self.rate_law == "CM":
                ln_eta_of_reaction = []
                for j in range(self.Nr):
                    B = ModelBalancing._B_matrix(self.Nc, S_subs[:, j], S_prod[:, j])
                    ln_eta_of_reaction += [
                        -scipy.special.logsumexp(
                            B @ (condition_ln_conc_met - ln_Km_matrix[:, j])
                        )
                    ]
                ln_eta_kinetic += [
                    np.reshape(np.vstack(ln_eta_of_reaction), (self.Nr,))
                ]
            else:
                raise ValueError(f"unsupported rate law {self.rate_law}")

        return np.vstack(ln_eta_kinetic).T

    @property
    def ln_eta_kinetic(self) -> np.ndarray:
        return self._ln_eta_kinetic(self.ln_conc_met, self.ln_Km)

    def _ln_eta_regulation(
        self,
        ln_conc_met: np.ndarray,
        ln_Ka: np.ndarray,
        ln_Ki: np.ndarray,
    ) -> np.ndarray:
        """Calculate the regulation (allosteric) term of the enzyme."""
        ln_Ka_matrix = ModelBalancing._create_dense_matrix(self.A_act, ln_Ka)
        ln_Ki_matrix = ModelBalancing._create_dense_matrix(self.A_inh, ln_Ki)

        ln_eta_allosteric = []
        for i in range(self.Ncond):
            condition_ln_conc_met = ln_conc_met[self.Nc*i:self.Nc*(i+1)]
            ln_conc_met_matrix = np.vstack([condition_ln_conc_met] * self.Nr).T
            ln_act = self.A_act.T @ self._logistic(ln_Ka_matrix - ln_conc_met_matrix)
            ln_inh = self.A_inh.T @ self._logistic(ln_conc_met_matrix - ln_Ki_matrix)
            ln_eta_act = -np.diag(ln_act)
            ln_eta_inh = -np.diag(ln_inh)
            ln_eta_allosteric += [ln_eta_act + ln_eta_inh]
        return np.vstack(ln_eta_allosteric).T

    @property
    def ln_eta_regulation(self) -> np.ndarray:
        return self._ln_eta_regulation(self.ln_conc_met, self.ln_Ka, self.ln_Ki)

    def _ln_conc_enz(
        self,
        ln_Keq: np.ndarray,
        ln_kcatf: np.ndarray,
        ln_Km: np.ndarray,
        ln_Ka: np.ndarray,
        ln_Ki: np.ndarray,
        ln_conc_met: np.ndarray,
    ) -> np.ndarray:
        """Calculate the required enzyme levels based on fluxes and rate laws."""
        driving_forces = self._driving_forces(ln_Keq, ln_conc_met)
        ln_capacity = self._ln_capacity(ln_kcatf)
        ln_eta_thermodynamic = self._ln_eta_thermodynamic(driving_forces)
        ln_eta_kinetic = self._ln_eta_kinetic(ln_conc_met, ln_Km)
        ln_eta_regulation = self._ln_eta_regulation(ln_conc_met, ln_Ka, ln_Ki)
        return ln_capacity - ln_eta_thermodynamic - ln_eta_kinetic - ln_eta_regulation

    @property
    def ln_conc_enz(self) -> np.ndarray:
        return self._ln_conc_enz(
            self.ln_Keq,
            self.ln_kcatf,
            self.ln_Km,
            self.ln_Ka,
            self.ln_Ki,
            self.ln_conc_met,
        )

    def is_gmean_feasible(self) -> bool:
        return (
            self._driving_forces(self.Keq_ln_gmean, self.conc_met_ln_gmean)
            >= self.MIN_DRIVING_FORCE
        ).all()

    def initialize_with_gmeans(self) -> None:
        # set the independent parameters values to the geometric means
        for p in ["Km", "Ka", "Ki", "Keq", "kcatf", "conc_met"]:
            self.__setattr__(f"ln_{p}", self.__getattribute__(f"{p}_ln_gmean"))

        # set the geometric means of the dependent parameters (kcatr and conc_enz)
        # to the values calculated using all the independent parameters
        for p in ["kcatr", "conc_enz"]:
            self.__setattr__(f"{p}_ln_gmean", self.__getattribute__(f"ln_{p}"))

    def solve(self) -> None:
        # in order to use the scipy solver, we need to stack all the independent variables
        # into one 1-D array (denoted 'x').

        x0 = []
        for p in ["Km", "Ka", "Ki", "Keq", "kcatf", "conc_met"]:
            x0.append(self.__getattribute__(f"ln_{p}"))
        x0 = np.hstack(x0)

        r = minimize(lambda x: self.objective_function(x), x0=x0)
        if not r.success:
            raise Exception(f"optimization unsuccessful because of {r.message}")

        i = 0
        for p in ["Km", "Ka", "Ki", "Keq", "kcatf", "conc_met"]:
            p_length = self.__getattribute__(f"ln_{p}_gmean").shape[0]
            self.__setattr__(f"ln_{p}", r.x[i:i+p_length])
            i += p_length
