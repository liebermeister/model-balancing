# -*- coding: utf-8 -*-
"""
This is my module brief line.

This is a more complete paragraph documenting my module.

- A list item.
- Another list item.

This section can use any reST syntax.
"""

import itertools
import os
import warnings
from typing import Any, List, Union

import cvxpy as cp
import numpy as np
from sbtab import SBtab

from . import (
    MIN_DRIVING_FORCE,
    MIN_FLUX,
    Q_,
    RT,
    DEFAULT_UNITS,
)

from .io import read_arguments_json, to_model_sbtab, to_state_sbtab


class ModelBalancingConvex(object):
    """A class for performing Convex Model Balancing (using CVXPY).
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

        for p in DEFAULT_UNITS.keys():
            assert f"{p}_gmean" in kwargs
            assert f"{p}_ln_cov" in kwargs

        self.ln_gmean = {}
        self.ln_precision = {}
        self.ln_lower_bound = {}
        self.ln_upper_bound = {}

        if not all(self.fluxes.flatten() > MIN_FLUX):
            raise ValueError(
                "In order to use Convex Model Balancing, all fluxes must be strictly positive"
                f"(i.e. > {MIN_FLUX:.1e})"
            )

        self.ln_Keq_gmean = cp.Parameter(
            shape=(self.Nr,), value=np.log(kwargs["Keq_gmean"].m_as(""))
        )
        self.ln_Keq_precision = np.linalg.pinv(kwargs["Keq_ln_cov"])
        self.ln_conc_met_gmean = cp.Parameter(
            shape=(self.Nc, self.Ncond),
            value=np.log(kwargs["conc_met_gmean"].m_as("M").reshape(self.Nc, self.Ncond)),
        )
        self.ln_conc_met_gstd = np.diag(kwargs["conc_met_ln_cov"] ** (0.5)).reshape(
            self.Nc, self.Ncond
        )
        self.ln_conc_enz_gmean = cp.Parameter(
            shape=(self.Nr, self.Ncond),
            value=np.log(kwargs["conc_enz_gmean"].m_as("M").reshape(self.Nr, self.Ncond)),
        )
        self.ln_conc_enz_gstd = np.diag(kwargs["conc_enz_ln_cov"] ** (0.5)).reshape(
            self.Nr, self.Ncond
        )
        self.ln_kcatf_gmean = cp.Parameter(
            shape=(self.Nr,), value=np.log(kwargs["kcatf_gmean"].m_as("1/s"))
        )
        self.ln_kcatf_precision = np.linalg.pinv(kwargs["kcatf_ln_cov"])
        self.ln_kcatr_gmean = cp.Parameter(
            shape=(self.Nr,), value=np.log(kwargs["kcatr_gmean"].m_as("1/s"))
        )
        self.ln_kcatr_precision = np.linalg.pinv(kwargs["kcatr_ln_cov"])
        self.ln_Km_gmean = np.log(kwargs["Km_gmean"].m_as("M"))
        if self.ln_Km_gmean.size != 0:
            self.ln_Km_precision = np.linalg.pinv(kwargs["Km_ln_cov"])
        else:
            self.ln_Km_precision = None

        self.ln_Ka_gmean = np.log(kwargs["Ka_gmean"].m_as("M"))
        if self.ln_Ka_gmean.size != 0:
            self.ln_Ka_precision = np.linalg.pinv(kwargs["Ka_ln_cov"])
        else:
            self.ln_Ka_precision = None

        self.ln_Ki_gmean = np.log(kwargs["Ki_gmean"].m_as("M"))
        if self.ln_Ki_gmean.size != 0:
            self.ln_Ki_precision = np.linalg.pinv(kwargs["Ki_ln_cov"])
        else:
            self.ln_Ki_precision = None

        assert self.ln_Keq_precision.shape == (self.Nr, self.Nr)
        assert self.ln_conc_enz_gstd.shape == (self.Nr, self.Ncond)
        assert self.ln_conc_met_gstd.shape == (self.Nc, self.Ncond)

        self.ln_conc_met = cp.Variable(shape=(self.Nc, self.Ncond))
        self.ln_Keq = cp.Variable(shape=(self.Nr,))
        self.ln_kcatf = cp.Variable(shape=(self.Nr,))

        for p in ["Km", "Ka", "Ki"]:
            ln_gmean = self.__getattribute__(f"ln_{p}_gmean")
            ln_precision = self.__getattribute__(f"ln_{p}_precision")
            if ln_gmean.size == 0:
                self.__setattr__(f"ln_{p}", np.array([]))
                self.__setattr__(f"z2_scores_{p}", cp.Constant(0))
            else:
                ln_p = cp.Variable(shape=ln_gmean.shape)
                self.__setattr__(
                    f"z2_scores_{p}",
                    ModelBalancingConvex._z_score(ln_p, ln_gmean, ln_precision),
                )
                self.__setattr__(f"ln_{p}", ln_p)

        for p in ["Keq", "kcatf", "kcatr"]:
            ln_p = self.__getattribute__(f"ln_{p}")
            ln_gmean = self.__getattribute__(f"ln_{p}_gmean")
            ln_precision = self.__getattribute__(f"ln_{p}_precision")
            self.__setattr__(
                f"z2_scores_{p}",
                ModelBalancingConvex._z_score(ln_p, ln_gmean, ln_precision),
            )

        # conc_met is given as a matrix (with conditions as columns) and therefore
        # we assume a diagonal covariance matrix (for simplicity). Instead of a
        # ln_cov matrix, we simply have the geometric means and stds arranged in
        # the same shape as the variables.
        self.z2_scores_conc_met = sum(
            cp.square(
                (self.ln_conc_met - self.ln_conc_met_gmean) / self.ln_conc_met_gstd
            ).flatten()
        )

        # ln enzyme concentrations are convex functions of the ln metabolite concentrations
        # but, since the z-scores use the square function, we have to take only the positive
        # values (otherwise the result is not convex).
        self.z2_scores_conc_enz = sum(
            cp.square(
                cp.pos(self.ln_conc_enz - self.ln_conc_enz_gmean)
                / self.ln_conc_enz_gstd
            ).flatten()
        )

        self.total_z2_scores = sum(
            [
                self.__getattribute__(f"z2_scores_{p}")
                for p in [
                    "Km",
                    "Ka",
                    "Ki",
                    "Keq",
                    "kcatf",
                    "kcatr",
                    "conc_met",
                    "conc_enz",
                ]
            ]
        )

    @staticmethod
    def from_json(fname: str) -> "ModelBalancingConvex":
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
        precision: np.array,
    ) -> cp.Expression:
        """Calculates the sum of squared Z-scores (with a covariance mat)."""
        return cp.quad_form(x - mu, precision)

    def _driving_forces(
        self,
        ln_Keq: Union[np.array, cp.Expression],
        ln_conc_met: Union[np.array, cp.Expression],
    ) -> cp.Expression:
        """Calculates the driving forces of all reactions."""
        return cp.vstack([ln_Keq] * self.Ncond).T - self.S.T @ ln_conc_met

    @property
    def driving_forces(self) -> cp.Expression:
        return self._driving_forces(self.ln_Keq, self.ln_conc_met)

    @staticmethod
    def _ln_kcatr(
        S: np.array,
        ln_kcatf: Union[np.array, cp.Expression],
        ln_Km: Union[np.array, cp.Expression],
        ln_Keq: Union[np.array, cp.Expression],
    ) -> cp.Expression:
        """Calculate the kcat-reverse based on Haldane relationship constraint."""
        ln_Km_matrix = ModelBalancingConvex._create_dense_matrix(S, ln_Km)
        return cp.diag(S.T @ ln_Km_matrix) + ln_kcatf - ln_Keq

    @property
    def ln_kcatr(self) -> cp.Expression:
        return ModelBalancingConvex._ln_kcatr(
            self.S, self.ln_kcatf, self.ln_Km, self.ln_Keq
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
        """Calculate the thermodynamic term of the enzyme."""
        return cp.log(1.0 - cp.exp(-driving_forces))

    @property
    def ln_eta_thermodynamic(self) -> cp.Expression:
        return self._ln_eta_thermodynamic(self.driving_forces)

    def _ln_eta_kinetic(
        self,
        ln_conc_met: Union[np.array, cp.Expression],
        ln_Km: Union[np.array, cp.Expression],
    ) -> cp.Expression:
        """Calculate the kinetic (saturation) term of the enzyme."""
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
        return self._ln_eta_kinetic(self.ln_conc_met, self.ln_Km)

    def _ln_eta_regulation(
        self,
        ln_conc_met: Union[np.array, cp.Expression],
        ln_Ka: Union[np.array, cp.Expression],
        ln_Ki: Union[np.array, cp.Expression],
    ) -> cp.Expression:
        """Calculate the regulation (allosteric) term of the enzyme."""
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
        """Calculate the required enzyme levels based on fluxes and rate laws."""
        driving_forces = self._driving_forces(ln_Keq, ln_conc_met)
        ln_capacity = self._ln_capacity(ln_kcatf)
        ln_eta_thermodynamic = self._ln_eta_thermodynamic(driving_forces)
        ln_eta_kinetic = self._ln_eta_kinetic(ln_conc_met, ln_Km)
        ln_eta_regulation = self._ln_eta_regulation(ln_conc_met, ln_Ka, ln_Ki)
        return ln_capacity - ln_eta_thermodynamic - ln_eta_kinetic - ln_eta_regulation

    @property
    def ln_conc_enz(self) -> cp.Expression:
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
            self._driving_forces(self.ln_Keq_gmean, self.ln_conc_met_gmean).value
            >= (MIN_DRIVING_FORCE / RT).m_as("")
        ).all()

    def initialize_with_gmeans(self) -> None:
        """Initialize the independent parameters with their gmeans.

        Note that the dependent parameters (kcatr and ln_conc_enz) can both
        be very far from their gmeans, and that the system might not be
        thermodynamically feasible.
        """
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

        # self.ln_kcatr_gmean.value = self.ln_kcatr.value
        # self.ln_conc_enz_gmean.value = self.ln_conc_enz.value

    def solve(self, verbose: bool = False) -> None:
        prob = cp.Problem(
            cp.Minimize(self.total_z2_scores),
            [self.driving_forces >= (MIN_DRIVING_FORCE / RT).m_as("")],
        )
        prob.solve(solver=self.solver, verbose=verbose)
        return prob.status

    def find_inner_point(self, verbose: bool = False) -> Any:
        prob = cp.Problem(
            cp.Minimize(self.z2_scores_conc_met),
            [self.driving_forces >= (MIN_DRIVING_FORCE / RT).m_as("")],
        )
        prob.solve(solver=self.solver, verbose=verbose)
        return prob.status

    @property
    def objective_value(self) -> float:
        return self.total_z2_scores.value

    def print_z_scores(self, precision: int = 2) -> None:
        for p in ["Km", "Ka", "Ki", "Keq", "kcatf", "conc_met", "kcatr", "conc_enz"]:
            z = self.__getattribute__(f"z2_scores_{p}").value
            print(f"{p} = {z.round(precision)}")

    def print_status(self) -> None:
        print("\nMetabolite concentrations (M) =\n", np.exp(self.ln_conc_met.value))
        print("\nEnzyme concentrations (M) =\n", np.exp(self.ln_conc_enz.value))
        print("\nDriving forces (RT) =\n", self.driving_forces.value)
        print("\nη(thr) =\n", np.exp(self.ln_eta_thermodynamic.value).round(2))
        print("\nη(kin) =\n", np.exp(self.ln_eta_kinetic.value).round(2))
        print("\nη(reg) =\n", np.exp(self.ln_eta_regulation.value).round(2))
        print("\n\n\n")

    def to_state_sbtab(self) -> SBtab.SBtabDocument:
        v = self.fluxes
        c = Q_(np.exp(self.ln_conc_met.value), "M")
        e = Q_(np.exp(self.ln_conc_enz.value), "M")
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
        kcatf = Q_(np.exp(self.ln_kcatf.value), "1/s")
        kcatr = Q_(np.exp(self.ln_kcatr.value), "1/s")
        Keq = Q_(np.exp(self.ln_Keq.value), "")
        if self.ln_Km_gmean.size != 0:
            Km = Q_(np.exp(self._create_dense_matrix(self.S, self.ln_Km).value), "M")
        else:
            Km = Q_(np.exp(self._create_dense_matrix(self.S, self.ln_Km)), "M")
        if self.ln_Ka_gmean.size != 0:
            Ka = Q_(
                np.exp(self._create_dense_matrix(self.A_act, self.ln_Ka).value), "M"
            )
        else:
            Ka = Q_(np.exp(self._create_dense_matrix(self.A_act, self.ln_Ka)), "M")
        if self.ln_Ki_gmean.size != 0:
            Ki = Q_(
                np.exp(self._create_dense_matrix(self.A_inh, self.ln_Ki).value), "M"
            )
        else:
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
        model_sbtabdoc.set_name("CMB result")
        return model_sbtabdoc

__all__ = [
    'ModelBalancingConvex'
]

