import itertools
import os
import warnings
from scipy.optimize import minimize

from typing import Dict, Tuple, Optional, List
import numpy as np
import scipy.special
from . import (
    read_arguments_json,
    to_state_sbtab,
    to_model_sbtab,
    SBtab,
    RT,
    Q_,
    MIN_DRIVING_FORCE,
    MIN_FLUX,
)


class ModelBalancing(object):

    INDEPENDENT_VARIABLES = ["Km", "Ka", "Ki", "Keq", "kcatf", "conc_met"]
    DEPENDENT_VARIABLES = ["kcatr", "conc_enz"]

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
        metabolite_names: List[str],
        reaction_names: List[str],
        state_names: List[str],
        rate_law: str = "CM",
        alpha: float = 1.0,
    ) -> None:
        self.S = S.copy()
        self.fluxes = fluxes.copy()
        self.A_act = A_act.copy()
        self.A_inh = A_inh.copy()

        self.Nc, self.Nr = S.shape
        assert self.fluxes.shape[0] == self.Nr
        if self.fluxes.ndim == 1:
            self.fluxes = self.fluxes.reshape(self.Nr, 1)
        self.Ncond = self.fluxes.shape[1]
        assert self.A_act.shape == (self.Nc, self.Nr)
        assert self.A_inh.shape == (self.Nc, self.Nr)

        # # since we use a log transform on the fluxes, we cannot allow any
        # # zero or negative values.
        # # First, we reverse the reactions that have negative fluxes
        # for i, j in zip(*np.where(self.fluxes < 0)):
        #     self.fluxes[i, j] = -self.fluxes[i, j]
        #     self.S[:, i] = -self.S[:, i]
        #     Keq_gmean[i] = 1.0 / Keq_gmean[i]
        #     kcatf_gmean[i], kcatr_gmean[i] = (
        #         kcatr_gmean[i],
        #         kcatf_gmean[i],
        #     )

        # For any flux that is too close to 0, we replace with a minimal value
        # of the same sign
        # for i, j in zip(*np.where(np.abs(self.fluxes) < MIN_FLUX)):
        #     self.fluxes[i, j] = np.sign(self.fluxes[i, j]) * MIN_FLUX

        self.ln_Keq_gmean = np.log(Keq_gmean.m_as(""))
        self.ln_Keq_precision = np.linalg.pinv(Keq_ln_cov)
        self.ln_kcatf_gmean = np.log(kcatf_gmean.m_as("1/s"))
        self.ln_kcatf_precision = np.linalg.pinv(kcatf_ln_cov)
        self.ln_kcatr_gmean = np.log(kcatr_gmean.m_as("1/s"))
        self.ln_kcatr_precision = np.linalg.pinv(kcatr_ln_cov)

        self.ln_Km_gmean = np.log(Km_gmean.m_as("M"))
        self.ln_Ka_gmean = np.log(Ka_gmean.m_as("M"))
        self.ln_Ki_gmean = np.log(Ki_gmean.m_as("M"))

        self.ln_Km_precision = np.linalg.pinv(Km_ln_cov)
        if Ka_ln_cov.size == 0:
            self.ln_Ka_precision = None
        else:
            self.ln_Ka_precision = np.linalg.pinv(Ka_ln_cov)
        if Ki_ln_cov.size == 0:
            self.ln_Ki_precision = None
        else:
            self.ln_Ki_precision = np.linalg.pinv(Ki_ln_cov)

        self.ln_conc_enz_gmean = np.log(conc_enz_gmean.m_as("M"))
        self.ln_conc_enz_precision = np.diag(
            np.log(conc_enz_gstd.T.flatten()) ** (-2.0)
        )
        self.ln_conc_met_gmean = np.log(conc_met_gmean.m_as("M"))
        self.ln_conc_met_precision = np.diag(
            np.log(conc_met_gstd.T.flatten()) ** (-2.0)
        )

        self.metabolite_names = metabolite_names
        self.reaction_names = reaction_names
        self.state_names = state_names

        self.rate_law = rate_law
        self.alpha = alpha

        for p in self.INDEPENDENT_VARIABLES:
            self.__setattr__(f"ln_{p}", self.__getattribute__(f"ln_{p}_gmean"))

    @staticmethod
    def from_json(fname: str) -> "ModelBalancing":
        args = read_arguments_json(fname)
        return ModelBalancing(**args)

    def _get_variable_shape(self, p: str) -> int:
        return self.__getattribute__(f"ln_{p}").shape

    def _get_variable_size(self, p: str) -> int:
        return self.__getattribute__(f"ln_{p}").size

    def _variable_vector_to_dict(
        self, x: Optional[np.ndarray] = None
    ) -> Dict[str, np.ndarray]:
        """Convert the variable vector into a dictionary."""
        if x is None:
            x = self._get_variable_vector()
        var_dict = {}
        i = 0
        for p in self.INDEPENDENT_VARIABLES:
            size = self._get_variable_size(p)
            shape = self._get_variable_shape(p)
            var_dict[f"ln_{p}"] = x[i : i + size].reshape(shape, order="F")
            i += size
        return var_dict

    def _get_variable_vector(self) -> np.ndarray:
        """Get the variable vector (x)."""
        # in order to use the scipy solver, we need to stack all the independent variables
        # into one 1-D array (denoted 'x').
        x0 = []
        for p in self.INDEPENDENT_VARIABLES:
            try:
                x0.append(self.__getattribute__(f"ln_{p}").T.flatten())
            except AttributeError:
                raise KeyError(f"Cannot find variable ln_{p}")
        return np.hstack(x0).flatten()

    def _get_full_variable_dictionary(
        self, x: Optional[np.ndarray] = None
    ) -> Dict[str, np.ndarray]:
        """Get a dictionary with all dependent and independent variable values."""
        var_dict = self._variable_vector_to_dict(x)
        var_dict["ln_conc_enz"] = self._ln_conc_enz(**var_dict).flatten()
        var_dict["ln_kcatr"] = ModelBalancing._ln_kcatr(
            self.S, var_dict["ln_kcatf"], var_dict["ln_Km"], var_dict["ln_Keq"]
        )
        return var_dict

    def objective_function(self, x: Optional[np.ndarray] = None) -> float:
        """Calculate the sum of squares of all Z-scores.

        The input (x) is a stacked version of all the independent variables, assuming
        the following order: Km, Ka, Ki, Keq, kcatf, conc_met
        """
        var_dict = self._get_full_variable_dictionary(x)

        total_z2_scores = 0.0
        for p in self.INDEPENDENT_VARIABLES + self.DEPENDENT_VARIABLES:
            ln_p_gmean = self.__getattribute__(f"ln_{p}_gmean")
            ln_p_precision = self.__getattribute__(f"ln_{p}_precision")
            ln_p = var_dict[f"ln_{p}"]

            if p == "conc_enz":
                # take a scaled version of the negative part of
                # the z-score of ln_conc_enz. (alpha = 0 would be convex, and alpha = 1
                # would be the true cost function)
                z2_score = ModelBalancing._z_score(
                    ln_p, ln_p_gmean, ln_p_precision, self.alpha
                )
            else:
                z2_score = ModelBalancing._z_score(ln_p, ln_p_gmean, ln_p_precision)
            total_z2_scores += z2_score

        return total_z2_scores

    @property
    def objective_value(self) -> float:
        return self.objective_function()

    @staticmethod
    def _z_score(
        x: np.array,
        mu: np.array,
        precision: np.array,
        alpha: float = 1,
    ) -> float:
        """Calculates the sum of squared Z-scores (with a covariance mat)."""
        if x.size == 0:
            return 0.0

        diff = x.flatten() - mu.flatten()

        full_z_score = diff.T @ precision @ diff

        if alpha == 1:
            return full_z_score
        else:
            pos_diff = np.array(diff)
            pos_diff[pos_diff < 0.0] = 0.0
            pos_z_score = pos_diff.T @ precision @ pos_diff
            return (1 - alpha) * pos_z_score + alpha * full_z_score

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
        return np.vstack([ln_Keq] * self.Ncond).T - self.S.T @ ln_conc_met

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
        ln_kcatr: np.ndarray,
    ) -> np.ndarray:
        """Calculate the capacity term of the enzyme."""

        # for positive fluxes we take the kcatf, and for negative ones we take
        # the kcatr. reactions with zero flux, are not counted at all.

        ln_cap = np.zeros((self.Nr, self.Ncond))
        for i in range(self.Nr):
            for j in range(self.Ncond):
                f = self.fluxes[i, j].m_as("M/s")
                if f > 0:
                    ln_cap[i, j] = np.log(f) - ln_kcatf[i]
                elif f < 0:
                    ln_cap[i, j] = np.log(-f) - ln_kcatr[i]

        return ln_cap

    def _ln_eta_thermodynamic(self, driving_forces: np.ndarray) -> np.ndarray:
        """Calculate the thermodynamic term of the enzyme."""
        eta_thermo = 1.0 - np.exp(-np.abs(driving_forces))

        # we ignore reactions that have exactly 0 flux (since we can assume
        # that their catalyzing enzyme is not expressed and therefore the
        # driving force can have any value (we do not need to impose any
        # probability distribution on it based on this reaction).
        eta_thermo[self.fluxes == 0] = 1.0
        return np.log(eta_thermo)

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

        ln_eta_kinetic = np.zeros((self.Nr, self.Ncond))
        for i in range(self.Nr):
            for j in range(self.Ncond):
                ln_c_over_Km = ln_conc_met[:, j] - ln_Km_matrix[:, i]
                try:
                    ln_1_plus_c_over_Km = np.log(1.0 + np.exp(ln_c_over_Km))
                except RuntimeWarning:
                    print(ln_c_over_Km)

                if self.rate_law == "S":
                    # numerator and denominator are identical
                    ln_eta_kinetic[i, j] = 0.0
                elif self.rate_law == "1S":
                    # ln(S / (1+S)) = -ln(1 + e^(-lnS))
                    ln_D_S = S_subs[:, i].T @ ln_c_over_Km
                    ln_eta_kinetic[i, j] = -np.log(1.0 + np.exp(-ln_D_S))
                elif self.rate_law == "SP":
                    # ln(S / (S+P)) = -ln(1 + e^(-lnS+lnP))
                    ln_D_SP = self.S[:, i].T @ ln_c_over_Km
                    ln_eta_kinetic[i, j] = -np.log(1.0 + np.exp(ln_D_SP))
                elif self.rate_law == "1SP":
                    # ln(S / (1+S+P)) = -ln(1 + e^(-lnS) + e^(-lnS+lnP))
                    ln_D_S = S_subs[:, i].T @ ln_c_over_Km
                    ln_D_SP = self.S[:, i].T @ ln_c_over_Km
                    ln_eta_kinetic[i, j] = -np.log(
                        1.0 + np.exp(-ln_D_S) + np.exp(ln_D_SP)
                    )
                elif self.rate_law == "CM":
                    # the common modular (CM) rate law
                    ln_D_S = S_subs[:, i].T @ ln_c_over_Km
                    ln_D_CM_S = S_subs[:, i].T @ ln_1_plus_c_over_Km
                    ln_D_CM_P = S_prod[:, i].T @ ln_1_plus_c_over_Km
                    ln_D_CM = np.log(np.exp(ln_D_CM_S) + np.exp(ln_D_CM_P) - 1.0)
                    ln_eta_kinetic[i, j] = ln_D_S - ln_D_CM
                else:
                    raise ValueError(f"unsupported rate law {self.rate_law}")

        return ln_eta_kinetic

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

        ln_eta_regulation = np.zeros((self.Nr, self.Ncond))
        for i in range(self.Nr):
            for j in range(self.Ncond):
                ln_c_over_Ka = ln_conc_met[:, j] - ln_Ka_matrix[:, i]
                ln_c_over_Ki = ln_conc_met[:, j] - ln_Ki_matrix[:, i]
                try:
                    ln_act = self.A_act[:, i].T @ np.log(1.0 + np.exp(-ln_c_over_Ka))
                    ln_inh = self.A_inh[:, i].T @ np.log(1.0 + np.exp(-ln_c_over_Ki))
                except RuntimeWarning:
                    print(ln_c_over_Ka)
                    print(ln_c_over_Ki)
                ln_eta_regulation[i, j] = ln_act + ln_inh
        return ln_eta_regulation

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
        ln_kcatr = self._ln_kcatr(self.S, ln_kcatf, ln_Km, ln_Keq)
        ln_capacity = self._ln_capacity(ln_kcatf, ln_kcatr)
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
            self._driving_forces(self.ln_Keq_gmean, self.ln_conc_met_gmean)
            >= (MIN_DRIVING_FORCE / RT).m_as("")
        ).all()

    def _thermodynamic_constraints(self) -> scipy.optimize.LinearConstraint:
        """Construct the thermodynamic constraints for the variable vector."""

        # given a constraint matrix (A), lower bound (lb), and a variable vector (x)
        # we want the following to two equations to hold:
        # 1) A @ x = np.vstack([ln_Keq] * self.Ncond).T - self.S.T @ ln_conc_met
        # 2) lb = MIN_DRIVING_FORCE
        # so that the thermodynamic constraint would be: A @ x >= lb

        # first, we find the indices of ln_Keq and ln_conc_met in x:
        i = 0
        i_Keq = None
        i_conc_met = None
        for p in self.INDEPENDENT_VARIABLES:
            if p == "Keq":
                i_Keq = i
            elif p == "conc_met":
                i_conc_met = i
            i += self._get_variable_size(p)

        A = np.zeros((self.Nr * self.Ncond, i))
        for j in range(self.Ncond):
            A[self.Nr * j : self.Nr * (j + 1), i_Keq : i_Keq + self.Nr] = np.eye(
                self.Nr
            )
            A[
                self.Nr * j : self.Nr * (j + 1),
                (i_conc_met + self.Nc * j) : (i_conc_met + self.Nc * (j + 1)),
            ] = -self.S.T

        lb = np.ones((self.Nr, self.Ncond)) * -np.inf
        ub = np.ones((self.Nr, self.Ncond)) * np.inf
        lb[self.fluxes > 0] = (MIN_DRIVING_FORCE / RT).m_as("")
        ub[self.fluxes < 0] = -(MIN_DRIVING_FORCE / RT).m_as("")
        lb = lb.T.reshape((self.Nr * self.Ncond,))
        ub = ub.T.reshape((self.Nr * self.Ncond,))
        # for i in range(self.Nr):
        #     for j in range(self.Ncond):
        #         f = self.fluxes[i, j].m_as("M/s")
        #         if f > 0:
        #             lb[i, j] = (MIN_DRIVING_FORCE / RT).m_as("")
        #         elif f < 0:
        #             ub[i, j] = -(MIN_DRIVING_FORCE / RT).m_as("")

        return scipy.optimize.LinearConstraint(A, lb, ub)

    def initialize_with_gmeans(self) -> None:
        # set the independent parameters values to the geometric means
        for p in self.INDEPENDENT_VARIABLES:
            self.__setattr__(f"ln_{p}", self.__getattribute__(f"ln_{p}_gmean"))

        # set the geometric means of the dependent parameters (kcatr and conc_enz)
        # to the values calculated using all the independent parameters
        for p in self.DEPENDENT_VARIABLES:
            self.__setattr__(f"ln_{p}_gmean", self.__getattribute__(f"ln_{p}"))

    def solve(self, solver: str = "SLSQP", options: Optional[dict] = None) -> None:

        x0 = self._get_variable_vector()

        # calculate the thermodynamic constraints (driving_forces >= MIN_DRIVING_FORCE)
        constraints = self._thermodynamic_constraints()

        r = minimize(
            fun=self.objective_function,
            x0=x0,
            constraints=constraints,
            method=solver,
            options=options,
        )
        if not r.success:
            raise Exception(f"optimization unsuccessful because of {r.message}")

        # copy the values from the solution to the class members
        for key, val in self._variable_vector_to_dict(r.x).items():
            self.__setattr__(key, val)

    def print_z_scores(self) -> None:
        for p in self.INDEPENDENT_VARIABLES + self.DEPENDENT_VARIABLES:
            ln_p_gmean = self.__getattribute__(f"ln_{p}_gmean")
            ln_p_precision = self.__getattribute__(f"ln_{p}_precision")
            ln_p = self.__getattribute__(f"ln_{p}")

            if p == "conc_enz":
                # take a scaled version of the negative part of
                # the z-score of ln_conc_enz. (alpha = 0 would be convex, and alpha = 1
                # would be the true cost function)
                z = ModelBalancing._z_score(
                    ln_p, ln_p_gmean, ln_p_precision, self.alpha
                )
            else:
                z = ModelBalancing._z_score(ln_p, ln_p_gmean, ln_p_precision)

            print(f"{p} = {z:.2f}")

    def print_status(self) -> None:
        print("\nMetabolite concentrations (M) =\n", np.exp(self.ln_conc_met))
        print("\nEnzyme concentrations (M) =\n", np.exp(self.ln_conc_enz))
        print("\nDriving forces (RT) =\n", self.driving_forces)
        print("\nη(thr) =\n", np.exp(self.ln_eta_thermodynamic).round(2))
        print("\nη(kin) =\n", np.exp(self.ln_eta_kinetic).round(2))
        print("\nη(reg) =\n", np.exp(self.ln_eta_regulation).round(2))
        print("\n\n\n")

    def to_state_sbtab(self) -> SBtab.SBtabDocument:
        v = self.fluxes
        c = Q_(np.exp(self.ln_conc_met), "M")
        e = Q_(np.exp(self.ln_conc_enz), "M")
        delta_g = -self.driving_forces * RT
        state_sbtabdoc = to_state_sbtab(
            v,
            c,
            e,
            delta_g,
            self.metabolite_names,
            self.reaction_names,
            self.state_names,
        )
        state_sbtabdoc.set_name("MB result")
        state_sbtabdoc.change_attribute("RelaxationAlpha", f"{self.alpha}")
        return state_sbtabdoc

    def to_model_sbtab(self) -> SBtab.SBtabDocument:
        kcatf = Q_(np.exp(self.ln_kcatf), "1/s")
        kcatr = Q_(np.exp(self.ln_kcatr), "1/s")
        Keq = Q_(np.exp(self.ln_Keq), "")
        Km = Q_(np.exp(self._create_dense_matrix(self.S, self.ln_Km)), "M")
        Ka = Q_(np.exp(self._create_dense_matrix(self.A_act, self.ln_Ka)), "M")
        Ki = Q_(np.exp(self._create_dense_matrix(self.A_inh, self.ln_Ki)), "M")
        model_sbtabdoc = to_model_sbtab(
            kcatf,
            kcatr,
            Keq,
            Km,
            Ka,
            Ki,
            self.S,
            self.A_act,
            self.A_inh,
            self.metabolite_names,
            self.reaction_names,
            self.state_names,
        )
        model_sbtabdoc.set_name("MB result")
        model_sbtabdoc.change_attribute("RelaxationAlpha", f"{self.alpha}")
        return model_sbtabdoc
