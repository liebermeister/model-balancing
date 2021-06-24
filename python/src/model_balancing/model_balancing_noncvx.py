"""
A module for performing full model balancing (using SciPy).

"""

import itertools
import os
import warnings
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import scipy.special
from sbtab import SBtab
from scipy.optimize import minimize

from . import (
    DEFAULT_UNITS,
    DEPENDENT_VARIABLES,
    INDEPENDENT_VARIABLES,
    MIN_DRIVING_FORCE,
    Q_,
    RT,
)
from .io import read_arguments_json, to_model_sbtab, to_state_sbtab


class ModelBalancing(object):
    """A class for performing Model Balancing (non-convex version).

    This version of model balancing solves the exact non-convex problem. The
    α parameter can be used to tune between a convex version (α = 0) and the
    full version (α = 1).
    """

    def __init__(
        self,
        S: np.array,
        fluxes: np.array,
        A_act: np.array,
        A_inh: np.array,
        metabolite_names: List[str],
        reaction_names: List[str],
        state_names: List[str],
        rate_law: str = "CM",
        alpha: float = 1.0,
        beta: float = 0.0,
        **kwargs,
    ) -> None:
        self.S = S.copy()
        self.fluxes = fluxes.copy()
        self.A_act = A_act.copy()
        self.A_inh = A_inh.copy()
        self.metabolite_names = metabolite_names
        self.reaction_names = reaction_names
        self.state_names = state_names
        self.rate_law = rate_law
        self.alpha = alpha
        self.beta = beta

        self.Nc, self.Nr = S.shape
        assert self.fluxes.shape[0] == self.Nr
        if self.fluxes.ndim == 1:
            self.fluxes = self.fluxes.reshape(self.Nr, 1)
        self.Ncond = self.fluxes.shape[1]
        assert self.A_act.shape == (self.Nc, self.Nr)
        assert self.A_inh.shape == (self.Nc, self.Nr)
        assert len(self.metabolite_names) == self.Nc
        assert len(self.reaction_names) == self.Nr
        assert len(self.state_names) == self.Ncond
        assert self.rate_law in [
            "S",
            "1S",
            "SP",
            "1SP",
            "CM",
        ], f"unsupported rate law {self.rate_law}"

        for p in DEFAULT_UNITS.keys():
            assert f"{p}_gmean" in kwargs
            assert f"{p}_ln_cov" in kwargs

        self.ln_gmean = {}
        self.ln_precision = {}
        self.ln_lower_bound = {}
        self.ln_upper_bound = {}

        for p, unit in DEFAULT_UNITS.items():
            # geometric means (in log-scale)
            self.ln_gmean[p] = np.log(kwargs[f"{p}_gmean"].m_as(unit))
            if self.ln_gmean[p].size > 0:
                self.ln_precision[p] = np.linalg.pinv(kwargs[f"{p}_ln_cov"])
                if np.all(kwargs[f"{p}_lower_bound"] != None):
                    self.ln_lower_bound[p] = np.log(
                        kwargs[f"{p}_lower_bound"].m_as(unit)
                    )
                else:
                    self.ln_lower_bound[p] = self.ln_gmean[
                        p
                    ].T.flatten() - 20.0 * np.diag(self.ln_precision[p]) ** (-1.0)
                if np.all(kwargs[f"{p}_upper_bound"] != None):
                    self.ln_upper_bound[p] = np.log(
                        kwargs[f"{p}_upper_bound"].m_as(unit)
                    )
                else:
                    self.ln_upper_bound[p] = self.ln_gmean[
                        p
                    ].T.flatten() + 20.0 * np.diag(self.ln_precision[p]) ** (-1.0)
            else:
                self.ln_precision[p] = None
                self.ln_lower_bound[p] = None
                self.ln_upper_bound[p] = None

        assert self.ln_gmean["conc_met"].shape == (self.Nc, self.Ncond)
        assert self.ln_precision["conc_met"].shape == (
            self.Nc * self.Ncond,
            self.Nc * self.Ncond,
        )
        assert self.ln_gmean["conc_enz"].shape == (self.Nr, self.Ncond)
        assert self.ln_precision["conc_enz"].shape == (
            self.Nr * self.Ncond,
            self.Nr * self.Ncond,
        )

        # initialize the independent variables with their geometric means
        for p in INDEPENDENT_VARIABLES:
            self.__setattr__(f"ln_{p}", self.ln_gmean[p])

    @staticmethod
    def from_json(fname: str) -> "ModelBalancing":
        """Initialize a ModelBalancing object using a JSON file.

        See our page about the :ref:`JSON specification sheet <json>`.
        """
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
            x = np.array(list(self.ln_x))
        var_dict = {}
        i = 0
        for p in INDEPENDENT_VARIABLES:
            size = self._get_variable_size(p)
            shape = self._get_variable_shape(p)
            var_dict[f"ln_{p}"] = x[i : i + size].reshape(shape, order="F")
            i += size
        return var_dict

    @property
    def ln_x(self) -> Iterable[float]:
        """Get the variable vector (x)."""
        # in order to use the scipy solver, we need to stack all the independent variables
        # into one 1-D array (denoted 'x').
        for p in INDEPENDENT_VARIABLES:
            for x in self.__getattribute__(f"ln_{p}").T.flat:
                yield x

    def _get_full_variable_dictionary(
        self, x: Optional[np.ndarray] = None
    ) -> Dict[str, np.ndarray]:
        """Get a dictionary with all dependent and independent variable values."""
        var_dict = self._variable_vector_to_dict(x)
        var_dict["ln_conc_enz"] = self._ln_conc_enz(**var_dict)
        var_dict["ln_kcatr"] = ModelBalancing._ln_kcatr(
            self.S, var_dict["ln_kcatf"], var_dict["ln_Km"], var_dict["ln_Keq"]
        )
        return var_dict

    def objective_function(self, x: Optional[np.ndarray] = None) -> float:
        """Calculate the sum of squares of all Z-scores for a given point (x).

        The input (x) is a stacked version of all the independent variables, assuming
        the following order: Km, Ka, Ki, Keq, kcatf, conc_met
        By default, if x is None, the objective value for (the exponent of)
        self.ln_x is returned.
        """
        var_dict = self._get_full_variable_dictionary(x)

        all_z2_scores = []
        for p in INDEPENDENT_VARIABLES + DEPENDENT_VARIABLES:
            ln_p_gmean = self.ln_gmean[p]
            ln_p_precision = self.ln_precision[p]
            ln_p = var_dict[f"ln_{p}"]

            # take a scaled version of the negative part of
            # the z-score of ln_conc_enz. (α = 0 would be convex, and α = 1
            # would be the true cost function)
            all_z2_scores.append(
                ModelBalancing._z_score(
                    ln_p,
                    ln_p_gmean,
                    ln_p_precision,
                    self.alpha if p == "conc_enz" else None,
                )
            )

        # add an extra term for the c/Km pseudo-parameters. we assume that
        # they are log-normal with a mean of 1 and stdev of 1/β
        if self.beta > 0:
            ln_Km_matrix = ModelBalancing._create_dense_matrix(
                self.S, var_dict["ln_Km"]
            )
            ln_c_over_Km = (
                np.repeat(var_dict["ln_conc_met"], (1, self.Nr)) - ln_Km_matrix
            )

            all_z2_scores.append(
                ModelBalancing._z_score(
                    ln_c_over_Km,
                    np.ones((self.Nc, self.Nr)),
                    1.0 / self.beta,
                )
            )

        return sum(all_z2_scores)

    @property
    def objective_value(self) -> float:
        """Get the objective value (i.e. the sum of squares of all z-scores)."""
        return self.objective_function()

    @staticmethod
    def _z_score(
        x: np.array,
        mu: np.array,
        precision: np.array,
        alpha: Optional[float] = None,
    ) -> float:
        """Calculates the sum of squared Z-scores (with a covariance mat).

        alpha - if given, scale the negative part of the z-score by alpha.
        if alpha is None, do not scale, which is equivalent to alpha = 1.
        (Default value is None)
        """
        if x.size == 0:
            return 0.0

        diff = x.flatten() - mu.flatten()

        full_z_score = diff.T @ precision @ diff

        if alpha is None:
            return full_z_score
        else:
            pos_diff = np.array(diff)
            pos_diff[pos_diff < 0.0] = 0.0
            pos_z_score = pos_diff.T @ precision @ pos_diff
            return (1.0 - alpha) * pos_z_score + alpha * full_z_score

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
        return np.vstack([ln_Keq] * self.Ncond).T - self.S.T @ ln_conc_met

    @property
    def driving_forces(self) -> np.ndarray:
        """Calculates the driving forces of all reactions."""
        return self._driving_forces(self.ln_Keq, self.ln_conc_met)

    @staticmethod
    def _ln_kcatr(
        S: np.array,
        ln_kcatf: np.ndarray,
        ln_Km: np.ndarray,
        ln_Keq: np.ndarray,
    ) -> np.ndarray:
        ln_Km_matrix = ModelBalancing._create_dense_matrix(S, ln_Km)
        return np.diag(S.T @ ln_Km_matrix) + ln_kcatf - ln_Keq

    @property
    def ln_kcatr(self) -> np.ndarray:
        """Calculate the kcat-reverse based on Haldane relationship constraint."""
        return ModelBalancing._ln_kcatr(self.S, self.ln_kcatf, self.ln_Km, self.ln_Keq)

    def _ln_capacity(
        self,
        ln_kcatf: np.ndarray,
        ln_kcatr: np.ndarray,
    ) -> np.ndarray:
        """Calculate the capacity term of the enzyme.

        for positive fluxes we take the kcatf, and for negative ones we take
        the kcatr. reactions with zero flux, are not counted at all.
        """

        ln_cap = np.zeros((self.Nr, self.Ncond))
        for i in range(self.Nr):
            for j in range(self.Ncond):
                f = self.fluxes[i, j].m_as("M/s")
                if f > 0.0:
                    ln_cap[i, j] = np.log(f) - ln_kcatf[i]
                elif f < 0.0:
                    ln_cap[i, j] = np.log(-f) - ln_kcatr[i]

        return ln_cap

    def _ln_eta_thermodynamic(self, driving_forces: np.ndarray) -> np.ndarray:
        eta_thermo = 1.0 - np.exp(-np.abs(driving_forces))

        # we ignore reactions that have exactly 0 flux (since we can assume
        # that their catalyzing enzyme is not expressed and therefore the
        # driving force can have any value (we do not need to impose any
        # probability distribution on it based on this reaction).
        for i, j in zip(*np.where(self.fluxes == 0)):
            eta_thermo[i, j] = 1.0
        return np.log(eta_thermo)

    @property
    def ln_eta_thermodynamic(self) -> np.ndarray:
        """Calculate the thermodynamic term of the enzyme."""
        return self._ln_eta_thermodynamic(self.driving_forces)

    def _ln_eta_kinetic(
        self,
        ln_conc_met: np.ndarray,
        ln_Km: np.ndarray,
    ) -> np.ndarray:
        S_neg = abs(self.S)
        S_pos = abs(self.S)
        S_neg[self.S > 0] = 0.0
        S_pos[self.S < 0] = 0.0

        ln_Km_matrix = ModelBalancing._create_dense_matrix(self.S, ln_Km)

        ln_eta_kinetic = np.zeros((self.Nr, self.Ncond))
        for i in range(self.Nr):
            for j in range(self.Ncond):
                if self.fluxes[i, j] > 0:
                    s_subs, s_prod = S_neg[:, i].T, S_pos[:, i].T
                elif self.fluxes[i, j] < 0:
                    s_subs, s_prod = S_pos[:, i].T, S_neg[:, i].T
                else:
                    continue

                ln_c_over_Km = ln_conc_met[:, j] - ln_Km_matrix[:, i]
                ln_1_plus_c_over_Km = np.log(1.0 + np.exp(ln_c_over_Km))

                if self.rate_law == "S":
                    # numerator and denominator are identical
                    ln_eta_kinetic[i, j] = 0.0
                elif self.rate_law == "1S":
                    # ln(S / (1+S)) = -ln(1 + e^(-lnS))
                    ln_D_S = s_subs @ ln_c_over_Km
                    ln_eta_kinetic[i, j] = -np.log(1.0 + np.exp(-ln_D_S))
                elif self.rate_law == "SP":
                    # ln(S / (S+P)) = -ln(1 + e^(-lnS+lnP))
                    ln_D_SP = (s_prod - s_subs) @ ln_c_over_Km
                    ln_eta_kinetic[i, j] = -np.log(1.0 + np.exp(ln_D_SP))
                elif self.rate_law == "1SP":
                    # ln(S / (1+S+P)) = -ln(1 + e^(-lnS) + e^(-lnS+lnP))
                    ln_D_S = s_subs @ ln_c_over_Km
                    ln_D_SP = (s_prod - s_subs) @ ln_c_over_Km
                    ln_eta_kinetic[i, j] = -np.log(
                        1.0 + np.exp(-ln_D_S) + np.exp(ln_D_SP)
                    )
                elif self.rate_law == "CM":
                    # the common modular (CM) rate law
                    ln_D_S = s_subs @ ln_c_over_Km
                    ln_D_CM_S = s_subs @ ln_1_plus_c_over_Km
                    ln_D_CM_P = s_prod @ ln_1_plus_c_over_Km
                    ln_D_CM = np.log(np.exp(ln_D_CM_S) + np.exp(ln_D_CM_P) - 1.0)
                    ln_eta_kinetic[i, j] = ln_D_S - ln_D_CM
                else:
                    raise ValueError(f"unsupported rate law {self.rate_law}")

        return ln_eta_kinetic

    @property
    def ln_eta_kinetic(self) -> np.ndarray:
        """Calculate the kinetic (saturation) term of the enzyme."""
        return self._ln_eta_kinetic(self.ln_conc_met, self.ln_Km)

    def _ln_eta_regulation(
        self,
        ln_conc_met: np.ndarray,
        ln_Ka: np.ndarray,
        ln_Ki: np.ndarray,
    ) -> np.ndarray:
        ln_eta_regulation = np.zeros((self.Nr, self.Ncond))
        if ln_Ka.size == 0 and ln_Ki.size == 0:
            return ln_eta_regulation

        ln_Ka_matrix = ModelBalancing._create_dense_matrix(self.A_act, ln_Ka)
        ln_Ki_matrix = ModelBalancing._create_dense_matrix(self.A_inh, ln_Ki)

        for i in range(self.Nr):
            for j in range(self.Ncond):
                if self.fluxes[i, j] == 0:
                    continue
                ln_c_over_Ka = ln_conc_met[:, j] - ln_Ka_matrix[:, i]
                ln_c_over_Ki = ln_conc_met[:, j] - ln_Ki_matrix[:, i]
                ln_act = self.A_act[:, i].T @ np.log(1.0 + np.exp(-ln_c_over_Ka))
                ln_inh = self.A_inh[:, i].T @ np.log(1.0 + np.exp(-ln_c_over_Ki))
                ln_eta_regulation[i, j] = ln_act + ln_inh
        return ln_eta_regulation

    @property
    def ln_eta_regulation(self) -> np.ndarray:
        """Calculate the regulation (allosteric) term of the enzyme."""
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
        driving_forces = self._driving_forces(ln_Keq, ln_conc_met)
        ln_kcatr = self._ln_kcatr(self.S, ln_kcatf, ln_Km, ln_Keq)
        ln_capacity = self._ln_capacity(ln_kcatf, ln_kcatr)
        ln_eta_thermodynamic = self._ln_eta_thermodynamic(driving_forces)
        ln_eta_kinetic = self._ln_eta_kinetic(ln_conc_met, ln_Km)
        ln_eta_regulation = self._ln_eta_regulation(ln_conc_met, ln_Ka, ln_Ki)
        return ln_capacity - ln_eta_thermodynamic - ln_eta_kinetic - ln_eta_regulation

    @property
    def ln_conc_enz(self) -> np.ndarray:
        """Calculate the required enzyme levels based on fluxes and rate laws."""
        return self._ln_conc_enz(
            self.ln_Keq,
            self.ln_kcatf,
            self.ln_Km,
            self.ln_Ka,
            self.ln_Ki,
            self.ln_conc_met,
        )

    def is_gmean_feasible(self) -> bool:
        """Check if the gmean  is a thermodynamically feasible solution.

        This is useful because we sometimes would like to initialize the
        optimization with the geometric means, but that can only be done if
        that point is feasible (otherwise, the dependent parameter
        conc_enz is not defined).
        """
        return (
            self._driving_forces(self.ln_gmean["Keq"], self.ln_gmean["conc_met"])
            >= (MIN_DRIVING_FORCE / RT).m_as("")
        ).all()

    @property
    def thermodynamic_constraints(self) -> scipy.optimize.LinearConstraint:
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
        for p in INDEPENDENT_VARIABLES:
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

        v = self.fluxes.T.reshape((self.Nr * self.Ncond,))
        lb = np.ones(self.Nr * self.Ncond) * -np.inf
        ub = np.ones(self.Nr * self.Ncond) * np.inf
        lb[v > 0] = (MIN_DRIVING_FORCE / RT).m_as("")
        ub[v < 0] = -(MIN_DRIVING_FORCE / RT).m_as("")

        # skip the rows in A, lb, ub which correspond to flux = 0, since these
        # reactions are unbounded above and below.
        return scipy.optimize.LinearConstraint(A[v != 0, :], lb[v != 0], ub[v != 0])

    def extend_variable_vector(self, x: np.ndarray) -> np.ndarray:
        """Add the dependent variables to a vector of the independen ones."""
        var_dict = self._get_full_variable_dictionary(x)
        var_array = []
        for p in INDEPENDENT_VARIABLES + DEPENDENT_VARIABLES:
            for x in var_dict[f"ln_{p}"].T.flat:
                var_array.append(x)
        return np.array(var_array)

    @property
    def parameter_constraints(self) -> scipy.optimize.NonlinearConstraint:
        lb = []
        ub = []
        for p in INDEPENDENT_VARIABLES + DEPENDENT_VARIABLES:
            if self.ln_lower_bound[p] is not None:
                for x in self.ln_lower_bound[p].T.flat:
                    lb.append(x)
            if self.ln_upper_bound[p] is not None:
                for x in self.ln_upper_bound[p].T.flat:
                    ub.append(x)
        lb = np.array(lb)
        ub = np.array(ub)
        return scipy.optimize.NonlinearConstraint(
            self.extend_variable_vector, lb=lb, ub=ub
        )

    def initialize_with_gmeans(self) -> None:
        """Initialize the independent parameters with their gmeans.

        Note that the dependent parameters (kcatr and ln_conc_enz) can both
        be very far from their gmeans, and that the system might not be
        thermodynamically feasible.
        """
        for p in INDEPENDENT_VARIABLES:
            self.__setattr__(f"ln_{p}", self.ln_gmean[p])

        # set the geometric means of the dependent parameters (kcatr and conc_enz)
        # to the values calculated using all the independent parameters
        for p in DEPENDENT_VARIABLES:
            self.ln_gmean[p] = self.__getattribute__(f"ln_{p}")

    def solve(self, solver: str = "SLSQP", options: Optional[dict] = None) -> None:
        """Find a local minimum of the objective function using SciPy."""
        x0 = np.array(list(self.ln_x))
        r = minimize(
            fun=self.objective_function,
            x0=x0,
            constraints=(self.thermodynamic_constraints, self.parameter_constraints),
            method=solver,
            options=options,
        )
        if not r.success:
            raise Exception(f"optimization unsuccessful because of {r.message}")

        # copy the values from the solution to the class members
        for key, val in self._variable_vector_to_dict(r.x).items():
            self.__setattr__(key, val)

    def get_z_scores(self) -> Dict[str, float]:
        """Get the z-score values for all variables."""
        res = {}
        for p in INDEPENDENT_VARIABLES + DEPENDENT_VARIABLES:
            ln_p_gmean = self.ln_gmean[p]
            ln_p_precision = self.ln_precision[p]
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

            res[p] = z
        return res

    def print_z_scores(self) -> None:
        """Print the z-score values for all variables."""
        for p, z in self.get_z_scores().items():
            print(f"{p} = {z:.2f}")

    def print_status(self) -> None:
        """Print a status report based on the current solution."""
        print("\nMetabolite concentrations (M) =\n", np.exp(self.ln_conc_met))
        print("\nEnzyme concentrations (M) =\n", np.exp(self.ln_conc_enz))
        print("\nDriving forces (RT) =\n", self.driving_forces)
        print("\nη(thr) =\n", np.exp(self.ln_eta_thermodynamic).round(2))
        print("\nη(kin) =\n", np.exp(self.ln_eta_kinetic).round(2))
        print("\nη(reg) =\n", np.exp(self.ln_eta_regulation).round(2))
        print("\n\n\n")

    def to_state_sbtab(self) -> SBtab.SBtabDocument:
        """Create a state SBtab.

        The state SBtab contains the values of the state-dependent variables,
        i.e. flux, concentrations of metabolites, concentrations of enzymes,
        and the ΔG' values.
        """
        v = self.fluxes
        c = Q_(np.exp(self.ln_conc_met), "M")
        e = Q_(np.exp(self.ln_conc_enz), "M")

        # during the optimization, reactions with 0 flux were not skipped in
        # the calculation of ln_capacity, ln_thermodynamic, ln_kinetic, and
        # ln_regulation. therefore, these reactions will have an enzyme
        # concentration of 1 (i.e. 0 in log-scale).
        # we need to replace that value with a 0 in order to return a correct
        # metabolic state.
        e[v == 0] = Q_(0.0, "M")

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
        """Create a model SBtab.

        The model SBtab contains the values of the state-independent variables,
        i.e. kcatf, kcatr, Km, Ka, and Ki.
        """
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


__all__ = ["ModelBalancing"]
