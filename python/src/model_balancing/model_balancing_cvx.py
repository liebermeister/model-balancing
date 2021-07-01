"""
A module for performing convex model balancing (using CVXPY).

"""

import itertools
import os
import warnings
from typing import Any, List, Union, Dict

import cvxpy as cp
import numpy as np
from sbtab import SBtab

from . import (
    DEFAULT_UNITS,
    DEPENDENT_VARIABLES,
    INDEPENDENT_VARIABLES,
    ALL_VARIABLES,
    STATE_VARIABLES,
    MODEL_VARIABLES,
    MIN_DRIVING_FORCE,
    MIN_FLUX,
    Q_,
    RT,
)
from .io import read_arguments_json, to_model_sbtab, to_state_sbtab


class NegativeFluxError(Exception):
    pass


class ModelBalancingConvex(object):
    """A class for performing Convex Model Balancing (using CVXPY).

    All input parameters are provided through the constructor. Then running
    .solve() will find the (global) optimum and the result will be stored in the
    cvxpy.Variable objects that represent the independent variables:
    .ln_Km .ln_Ka .ln_Ki .ln_Keq .ln_kcatf .ln_conc_met
    and also in the two dependent cvxpy.Expression objects:
    .ln_kcatr .ln_conc_enz

    Use the .value to get the solution for any of these Variables/Expressions.
    Alternatively, .to_model_sbtab and .to_state_sbtab return SBtab objects
    containing the state (met_conc, enz_conc) and the model (Km, Ki, Ka, Keq,
    kcatf, kcatr) values respectively.

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
        solver: str = "MOSEK",
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
        self.beta = beta
        self.solver = solver

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

        for p in ALL_VARIABLES:
            assert f"{p}_gmean" in kwargs
            assert f"{p}_ln_precision" in kwargs

        self._var_dict = {}
        self.ln_gmean = {}
        self.z2_scores = {}
        self.ln_lower_bound = {}
        self.ln_upper_bound = {}

        if not all(self.fluxes.flatten() > MIN_FLUX):
            raise NegativeFluxError(
                "In order to use Convex Model Balancing, all fluxes must be strictly positive"
                f"(i.e. > {MIN_FLUX:.1e})"
            )

        for p in INDEPENDENT_VARIABLES + DEPENDENT_VARIABLES:
            p_dim = kwargs[f"{p}_gmean"].size
            if p_dim > 0:
                if p in INDEPENDENT_VARIABLES:
                    self._var_dict[f"ln_{p}"] = cp.Variable(
                        shape=kwargs[f"{p}_gmean"].shape
                    )
                elif p == "conc_enz":
                    self._var_dict[f"ln_{p}"] = self.ln_conc_enz
                elif p == "kcatr":
                    self._var_dict[f"ln_{p}"] = self.ln_kcatr
                else:
                    raise Exception(f"unknown parameter: {p}")

                self.ln_gmean[p] = np.log(kwargs[f"{p}_gmean"].m_as(DEFAULT_UNITS[p]))
                assert kwargs[f"{p}_ln_precision"].shape == (p_dim, p_dim)
                displacement = self._var_dict[f"ln_{p}"] - self.ln_gmean[p]
                if len(displacement.shape) > 1:
                    displacement = displacement.flatten()
                if p == "conc_enz":
                    displacement = cp.pos(displacement)

                self.z2_scores[p] = cp.quad_form(
                    displacement, kwargs[f"{p}_ln_precision"]
                )
            else:
                self.ln_gmean[p] = None
                self._var_dict[f"ln_{p}"] = None
                self.z2_scores[p] = cp.Constant(0)

        if self.beta > 0:
            self.z2_scores["c_over_Km"] = self.beta * self.penalty_term_beta(
                self._var_dict["ln_Km"], self._var_dict["ln_conc_met"]
            )

        self.total_z2_scores = sum(self.z2_scores.values())

    @staticmethod
    def from_json(fname: str) -> "ModelBalancingConvex":
        """Initialize a ModelBalancingConvex object using a JSON file.

        See our page about the :ref:`JSON specification sheet <json>`.
        """
        args = read_arguments_json(fname)
        return ModelBalancingConvex(**args)

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
    def _create_dense_matrix(
        S: np.array,
        x: Union[np.array, cp.Variable],
    ) -> cp.Expression:
        """Converts a sparse list of affinity parameters (e.g. Km) to a matrix."""

        Nc, Nr = S.shape

        if x is None:
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

    def _driving_forces(
        self,
        ln_Keq: Union[np.array, cp.Expression],
        ln_conc_met: Union[np.array, cp.Expression],
    ) -> cp.Expression:
        return cp.vstack([ln_Keq] * self.Ncond).T - self.S.T @ ln_conc_met

    @property
    def driving_forces(self) -> cp.Expression:
        """Calculates the driving forces of all reactions."""
        return self._driving_forces(
            self._var_dict["ln_Keq"], self._var_dict["ln_conc_met"]
        )

    @staticmethod
    def _ln_kcatr(
        S: np.array,
        ln_kcatf: Union[np.array, cp.Expression],
        ln_Km: Union[np.array, cp.Expression],
        ln_Keq: Union[np.array, cp.Expression],
    ) -> cp.Expression:
        ln_Km_matrix = ModelBalancingConvex._create_dense_matrix(S, ln_Km)
        return cp.diag(S.T @ ln_Km_matrix) + ln_kcatf - ln_Keq

    @property
    def ln_kcatr(self) -> cp.Expression:
        """Calculate the kcat-reverse based on Haldane relationship constraint."""
        return ModelBalancingConvex._ln_kcatr(
            self.S,
            self._var_dict["ln_kcatf"],
            self._var_dict["ln_Km"],
            self._var_dict["ln_Keq"],
        )

    def _ln_capacity(
        self,
        ln_kcatf: Union[np.array, cp.Expression],
    ) -> cp.Expression:
        """Calculate the capacity term of the enzyme."""
        return np.log(self.fluxes.m_as("M/s")) - cp.vstack([ln_kcatf] * self.Ncond).T

    def _ln_eta_thermodynamic(
        self, driving_forces: Union[np.array, cp.Expression]
    ) -> cp.Expression:
        return cp.log(1.0 - cp.exp(-driving_forces))

    @property
    def ln_eta_thermodynamic(self) -> cp.Expression:
        """Calculate the thermodynamic term of the enzyme :math:`\eta^{thr}`."""
        return self._ln_eta_thermodynamic(self.driving_forces)

    def _ln_eta_kinetic(
        self,
        ln_conc_met: Union[np.array, cp.Expression],
        ln_Km: Union[np.array, cp.Expression],
    ) -> cp.Expression:
        S_subs = abs(self.S)
        S_prod = abs(self.S)
        S_subs[self.S > 0] = 0
        S_prod[self.S < 0] = 0

        ln_Km_matrix = ModelBalancingConvex._create_dense_matrix(self.S, ln_Km)

        ln_eta_kinetic = []
        for i in range(self.Ncond):
            ln_conc_met_matrix = cp.vstack([ln_conc_met[:, i]] * self.Nr).T
            ln_D_S = S_subs.T @ (ln_conc_met_matrix - ln_Km_matrix)
            ln_D_P = S_prod.T @ (ln_conc_met_matrix - ln_Km_matrix)

            if self.rate_law == "S":
                ln_eta_kinetic += [
                    cp.Constant(
                        np.zeros(
                            self.Nr,
                        )
                    )
                ]
            elif self.rate_law == "1S":
                ln_eta_kinetic += [-cp.diag(cp.logistic(-ln_D_S))]
            elif self.rate_law == "SP":
                ln_eta_kinetic += [-cp.diag(cp.logistic(-ln_D_S + ln_D_P))]
            elif self.rate_law == "1SP":
                ln_eta_kinetic += [-cp.diag(cp.logistic(-ln_D_S + cp.logistic(ln_D_P)))]
            elif self.rate_law == "CM":
                ln_eta_of_reaction = []
                for j in range(self.Nr):
                    B = ModelBalancingConvex._B_matrix(
                        self.Nc, S_subs[:, j], S_prod[:, j]
                    )
                    ln_eta_of_reaction += [
                        -cp.log_sum_exp(B @ (ln_conc_met[:, i] - ln_Km_matrix[:, j]))
                    ]
                ln_eta_kinetic += [
                    cp.reshape(cp.vstack(ln_eta_of_reaction), (self.Nr,))
                ]
            else:
                raise ValueError(f"unsupported rate law {self.rate_law}")

        return cp.vstack(ln_eta_kinetic).T

    @property
    def ln_eta_kinetic(self) -> cp.Expression:
        """Calculate the kinetic (saturation) term of the enzyme :math:`\eta^{kin}`."""
        return self._ln_eta_kinetic(self.ln_conc_met, self.ln_Km)

    def _ln_eta_regulation(
        self,
        ln_conc_met: Union[np.array, cp.Expression],
        ln_Ka: Union[np.array, cp.Expression],
        ln_Ki: Union[np.array, cp.Expression],
    ) -> cp.Expression:
        ln_Ka_matrix = ModelBalancingConvex._create_dense_matrix(self.A_act, ln_Ka)
        ln_Ki_matrix = ModelBalancingConvex._create_dense_matrix(self.A_inh, ln_Ki)

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
    def ln_eta_regulation(self) -> cp.Expression:
        """Calculate the regulation (allosteric) term of the enzyme :math:`\eta^{reg}`."""
        return self._ln_eta_regulation(self.ln_conc_met, self.ln_Ka, self.ln_Ki)

    def _ln_conc_enz(
        self,
        ln_Keq: Union[np.array, cp.Expression],
        ln_kcatf: Union[np.array, cp.Expression],
        ln_Km: Union[np.array, cp.Expression],
        ln_Ka: Union[np.array, cp.Expression],
        ln_Ki: Union[np.array, cp.Expression],
        ln_conc_met: Union[np.array, cp.Expression],
    ) -> cp.Expression:
        driving_forces = self._driving_forces(ln_Keq, ln_conc_met)
        ln_capacity = self._ln_capacity(ln_kcatf)
        ln_eta_thermodynamic = self._ln_eta_thermodynamic(driving_forces)
        ln_eta_kinetic = self._ln_eta_kinetic(ln_conc_met, ln_Km)
        ln_eta_regulation = self._ln_eta_regulation(ln_conc_met, ln_Ka, ln_Ki)
        return ln_capacity - ln_eta_thermodynamic - ln_eta_kinetic - ln_eta_regulation

    @property
    def ln_conc_enz(self) -> cp.Expression:
        """Calculate the required enzyme levels based on fluxes and rate laws."""
        kwargs = {f"ln_{p}": self._var_dict[f"ln_{p}"] for p in INDEPENDENT_VARIABLES}
        return self._ln_conc_enz(**kwargs)

    def penalty_term_beta(
        self,
        ln_Km: Union[np.array, cp.Expression],
        ln_conc_met: Union[np.array, cp.Expression],
    ) -> cp.Expression:
        """Calculate the penalty term for c/Km."""
        beta_z_score = 0.0
        ln_Km_matrix = ModelBalancingConvex._create_dense_matrix(self.S, ln_Km)
        for i in range(self.Nc):
            for j in range(self.Nr):
                if self.S[i, j] == 0:
                    continue
                for k in range(self.Ncond):
                    ln_c_minus_km = ln_conc_met[i, k] - ln_Km_matrix[i, j]
                    beta_z_score += ln_c_minus_km ** 2
        return beta_z_score

    def is_gmean_feasible(self) -> bool:
        """Check if the gmean  is a thermodynamically feasible solution.

        This is useful because we sometimes would like to initialize the
        optimization with the geometric means, but that can only be done if
        that point is feasible (otherwise, the dependent parameter
        conc_enz is not defined).
        """
        return (
            self._driving_forces(self.ln_gmean["Keq"], self.ln_gmean["conc_met"]).value
            >= (MIN_DRIVING_FORCE / RT).m_as("")
        ).all()

    def initialize_with_gmeans(self) -> None:
        """Initialize the independent parameters with their gmeans.

        Note that the dependent parameters (kcatr and ln_conc_enz) can both
        be very far from their gmeans, and that the system might not be
        thermodynamically feasible.
        """
        # define the independent variables and their z-scores
        for p in INDEPENDENT_VARIABLES:
            if self.ln_gmean[p] is not None:
                self._var_dict[f"ln_{p}"].value = self.ln_gmean[p]

        # self.ln_kcatr_gmean.value = self.ln_kcatr.value
        # self.ln_conc_enz_gmean.value = self.ln_conc_enz.value

    def solve(self, verbose: bool = False) -> str:
        """Use CVXPY to find the global optimum that minimizes all z-scores."""
        prob = cp.Problem(
            cp.Minimize(self.total_z2_scores),
            [self.driving_forces >= (MIN_DRIVING_FORCE / RT).m_as("")],
        )
        prob.solve(solver=self.solver, verbose=verbose)
        return prob.status

    def find_inner_point(self, verbose: bool = False) -> Any:
        """Find a point insize the feasible thermodynamic space.

        If the geometric mean is not a feasible solution, we will need this
        function in order to initialize the solver with a point inside the
        thermodynamically feasible space.
        """
        prob = cp.Problem(
            cp.Minimize(self.z2_scores["conc_met"]),
            [self.driving_forces >= (MIN_DRIVING_FORCE / RT).m_as("")],
        )
        prob.solve(solver=self.solver, verbose=verbose)
        return prob.status

    @property
    def objective_value(self) -> float:
        """Get the objective value (i.e. the sum of squares of all z-scores)."""
        return self.total_z2_scores.value

    def get_z_scores(self) -> Dict[str, np.array]:
        return {k: v.value for k, v in self.z2_scores.items()}

    def print_z_scores(self, precision: int = 2) -> None:
        """Print the z-score values for all variables."""
        for p in ALL_VARIABLES:
            z = self.z2_scores[p].value
            print(f"{p} = {z.round(precision)}")

    def print_status(self) -> None:
        """Print a status report based on the current solution."""
        print("\nMetabolite concentrations (M) =\n", np.exp(self.ln_conc_met.value))
        print("\nEnzyme concentrations (M) =\n", np.exp(self.ln_conc_enz.value))
        print("\nDriving forces (RT) =\n", self.driving_forces.value)
        print("\nη(thr) =\n", np.exp(self.ln_eta_thermodynamic.value).round(2))
        print("\nη(kin) =\n", np.exp(self.ln_eta_kinetic.value).round(2))
        print("\nη(reg) =\n", np.exp(self.ln_eta_regulation.value).round(2))
        print("\n\n\n")

    def to_state_sbtab(self) -> SBtab.SBtabDocument:
        """Create a state SBtab.

        The state SBtab contains the values of the state-dependent variables,
        i.e. flux, concentrations of metabolites, concentrations of enzymes,
        and the ΔG' values.
        """
        v = self.fluxes
        c = Q_(np.exp(self._var_dict["ln_conc_met"].value), "M")
        e = Q_(np.exp(self._var_dict["ln_conc_enz"].value), "M")
        delta_g = -self.driving_forces.value * RT
        state_sbtabdoc = to_state_sbtab(
            v,
            c,
            e,
            delta_g,
            self.metabolite_names,
            self.reaction_names,
            self.state_names,
        )
        state_sbtabdoc.set_name("CMB result")
        return state_sbtabdoc

    def to_model_sbtab(self) -> SBtab.SBtabDocument:
        """Create a model SBtab.

        The model SBtab contains the values of the state-independent variables,
        i.e. kcatf, kcatr, Km, Ka, and Ki.
        """
        kwargs = {
            "S": self.S,
            "A_act": self.A_act,
            "A_inh": self.A_inh,
            "metabolite_names": self.metabolite_names,
            "reaction_names": self.reaction_names,
            "state_names": self.state_names,
        }
        for p in MODEL_VARIABLES:
            if self.ln_gmean[p] is None:
                kwargs[p] = Q_(np.array([]), DEFAULT_UNITS[p])
                continue
            val = np.exp(self._var_dict[f"ln_{p}"].value)
            if p == "Km":
                val = self._create_dense_matrix(self.S, val)
            elif p == "Ka":
                val = self._create_dense_matrix(self.A_act, val)
            elif p == "Ki":
                val = self._create_dense_matrix(self.A_inh, val)
            kwargs[p] = Q_(val, DEFAULT_UNITS[p])

        model_sbtabdoc = to_model_sbtab(**kwargs)
        model_sbtabdoc.set_name("CMB result")
        return model_sbtabdoc


__all__ = ["ModelBalancingConvex"]
