"""
This module containts a few mathematical utility functions for model balancing.
"""

import itertools
import os
import warnings
from typing import Union

import cvxpy as cp
import numpy as np
import pint

ureg = pint.UnitRegistry(system="mks")
Q_ = ureg.Quantity

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Q_([])


def B_matrix(Nc: int, col_subs: np.ndarray, col_prod: np.ndarray) -> np.ndarray:
    """Build the B matrix for the :math:`\eta^{kin}` expression.

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


def logistic(x: np.ndarray) -> np.ndarray:
    """elementwise calculation of :math:`\log(1 + e^x)`"""
    return np.log(1.0 + np.exp(x))


__all__ = [
    "logistic",
    "B_matrix",
]
