#!/usr/bin/env python3
"""Convert model balancing input files (SBtab) into the JSON problem file.

This is a stand-alone Python port of the Matlab code path

    mb_load_and_convert_model_and_data -> mb_model_and_data
      -> mb_prepare_posterior -> cvxpy_problem_data_structure

so that the Python version of model balancing can start from the same input files as
the Matlab version. It needs numpy only. Input are SBtab files describing network,
kinetic data and state data, plus (optionally) concentration constraints; the prior
table is bundled in data/mb_prior.tsv.

The JSON file has the structure documented in docs/source/json_format.rst of the Python
package, with these deliberate differences from the file written by Matlab:
  - arrays keep their shape: a (Nm, 1) matrix or a vector of length 1 stays an array
    (Matlab's jsonencode turns them into a flat vector or a number)

Command line:
    python sbtab_to_json.py network.tsv kinetic_data.tsv state_data.tsv out.json \\
        --constraints constraints.tsv --mean S1_Mean S2_Mean --std S1_Std S2_Std \\
        --fill-true-values

From Python:
    problem = build_problem("network.tsv", "kinetic_data.tsv", "state_data.tsv",
                            columns_mean=["S1_Mean", "S2_Mean"])
"""

import argparse
import json
import math
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

__version__ = "0.1.0"

DEFAULT_PRIOR_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "mb_prior.tsv")

R_GAS = 0.0083144621  # kJ/(mol K), as in Matlab function RT
TEMPERATURE = 298.15  # K
RT = R_GAS * TEMPERATURE

STANDARD_CONCENTRATION = "1 mM"

# fixed values from mb_default_options
DATA_KIN_GEOM_STD = 1.5
DATA_V_GEOM_STD = 1.2
FLUX_MAX = 100.0  # mM/s
KM_MAX = 1000.0  # mM, default upper limit for KM values (data and bounds); user option km_max
CONC_EPSILON = 1e-5  # used in mb_make_bounds if a lower and upper bound are identical
LOG_NORMAL_MIN_MEAN = 1e-14  # used in lognormal_normal2log

KINETIC_QUANTITIES = [
    "standard Gibbs free energy of reaction",
    "standard chemical potential",
    "Michaelis constant",
    "activation constant",
    "inhibitory constant",
    "equilibrium constant",
    "substrate catalytic rate constant",
    "product catalytic rate constant",
]

QUANTITY_SYNONYMS = {
    "forward catalytic rate constant": "substrate catalytic rate constant",
    "reverse catalytic rate constant": "product catalytic rate constant",
}


# --------------------------------------------------------------------------------
# SBtab reader


@dataclass
class Table:
    """One SBtab table: attributes of the declaration line and columns by name (without '!')."""

    attributes: Dict[str, str]
    columns: Dict[str, List[str]] = field(default_factory=dict)

    @property
    def nrows(self) -> int:
        return len(next(iter(self.columns.values()))) if self.columns else 0

    def has(self, name: str) -> bool:
        return name in self.columns

    def get(self, name: str) -> List[str]:
        return self.columns[name]


def _split_tables(path: str) -> List[Table]:
    tables: List[Table] = []
    header: Optional[List[str]] = None
    with open(path, "rt", encoding="utf-8-sig", newline=None) as fp:
        for raw in fp:
            cells = [c.strip() for c in raw.rstrip("\r\n").split("\t")]
            first = cells[0]
            if first.startswith("!!!"):
                continue
            if first.startswith("!!"):
                attributes = dict(re.findall(r"(\w+)\s*=\s*'([^']*)'", "\t".join(cells)))
                tables.append(Table(attributes))
                header = None
            elif first.startswith("%") or not any(cells):
                continue
            elif not tables:
                continue
            elif header is None:
                if first.startswith("!"):
                    header = [c[1:] if c.startswith("!") else c for c in cells]
                    tables[-1].columns = {h: [] for h in header if h}
            else:
                row = cells + [""] * (len(header) - len(cells))
                for name, value in zip(header, row):
                    if name:
                        tables[-1].columns[name].append(value)
    return tables


def read_sbtab_files(paths: Sequence[str]) -> Dict[str, Table]:
    """Read several SBtab files into one dictionary of tables (by TableID).

    If a TableID appears twice, the later table is used.
    """
    document: Dict[str, Table] = {}
    for path in paths:
        for table in _split_tables(path):
            table_id = table.attributes.get("TableID")
            if table_id:
                document[table_id] = table
    return document


def to_float(text: str) -> float:
    """Convert a table entry to a number (empty or non-numeric entries give NaN)."""
    try:
        return float(text)
    except (TypeError, ValueError):
        return math.nan


def to_floats(cells: Sequence[str]) -> np.ndarray:
    return np.array([to_float(c) for c in cells], dtype=float)


def index_case_insensitive(names: Sequence[str], all_names: Sequence[str]) -> np.ndarray:
    """Index of each name in all_names (first hit, ignoring case), or -1 (Matlab: label_names)."""
    lookup: Dict[str, int] = {}
    for i, name in enumerate(all_names):
        if name:
            lookup.setdefault(name.upper(), i)
    return np.array([lookup.get(n.upper(), -1) for n in names], dtype=int)


# --------------------------------------------------------------------------------
# Prior table and options (Matlab: mb_default_options)


@dataclass
class Quantity:
    """Row of the prior table."""

    symbol: str
    element: str
    scaling: str
    unit: str
    median: float
    geom_std: float
    lower: float
    upper: float
    data_std: float
    data_geom_std: float


def read_prior_table(path: str) -> Dict[str, Quantity]:
    """Read the prior table; the result is indexed by quantity type and by symbol."""
    table = _split_tables(path)[0]
    result: Dict[str, Quantity] = {}
    for i in range(table.nrows):
        q = Quantity(
            symbol=table.get("Symbol")[i],
            element=table.get("BiologicalElement")[i],
            scaling=table.get("MathematicalType")[i],
            unit=table.get("Unit")[i],
            median=to_float(table.get("PriorMedian")[i]),
            geom_std=to_float(table.get("PriorGeometricStd")[i]),
            lower=to_float(table.get("LowerBound")[i]),
            upper=to_float(table.get("UpperBound")[i]),
            data_std=to_float(table.get("DataStd")[i]),
            data_geom_std=to_float(table.get("DataGeometricStd")[i]),
        )
        result[table.get("QuantityType")[i]] = q
        result[q.symbol] = q
    return result


@dataclass
class GaussianPrior:
    mean_ln: float
    std_ln: float


def prior_of(prior_table: Dict[str, Quantity], symbol: str) -> GaussianPrior:
    q = prior_table[symbol]
    return GaussianPrior(math.log(q.median), math.log(q.geom_std))


# --------------------------------------------------------------------------------
# Network (Matlab: sbtab_to_network, network_build_from_sum_formulae)


@dataclass
class Network:
    metabolites: List[str]
    reactions: List[str]
    N: np.ndarray  # metabolites x reactions
    regulation: np.ndarray  # reactions x metabolites; +1 activation, -1 inhibition
    metabolite_kegg: List[str]
    reaction_kegg: List[str]
    metabolite_bigg: List[str]
    reaction_bigg: List[str]

    @property
    def nr(self) -> int:
        return len(self.reactions)

    @property
    def nm(self) -> int:
        return len(self.metabolites)

    # Reaction and metabolite indices of the elements of KM, KA, KI: sorted by metabolite,
    # then by reaction (like Matlab's find in a matrix reactions x metabolites)

    @property
    def km_indices(self) -> Tuple[np.ndarray, np.ndarray]:
        m, r = np.nonzero(self.N)
        return r, m

    @property
    def ka_indices(self) -> Tuple[np.ndarray, np.ndarray]:
        m, r = np.nonzero((self.regulation > 0).T)
        return r, m

    @property
    def ki_indices(self) -> Tuple[np.ndarray, np.ndarray]:
        m, r = np.nonzero((self.regulation < 0).T)
        return r, m


def _parse_side(side: str) -> List[Tuple[float, str]]:
    side = side.strip()
    if not side:
        return []
    terms = []
    for term in side.split(" + "):
        term = term.strip()
        if " " in term:
            coefficient, metabolite = term.split(" ", 1)
            terms.append((float(coefficient), metabolite.strip()))
        elif term:
            terms.append((1.0, term))
    return terms


def _parse_formula(formula: str) -> Tuple[List[Tuple[float, str]], List[Tuple[float, str]]]:
    if "<=>" not in formula:
        raise ValueError(f'Malformed reaction formula "{formula}"')
    if "," in formula:
        raise ValueError(f'Malformed reaction formula "{formula}"')
    left, right = formula.split("<=>", 1)
    return _parse_side(left), _parse_side(right)


def _parse_regulation(text: str) -> List[Tuple[float, str]]:
    """Parse an entry of column MetabolicRegulation, for example '+ A - B'."""
    return [(1.0 if m[0] == "+" else -1.0, m[1]) for m in re.findall(r"(?:^|\s)([+-])\s+(\S+)", text)]


def _column_or_empty(table: Table, name: str, n: int) -> List[str]:
    return table.get(name) if table.has(name) else [""] * n


def build_network(document: Dict[str, Table]) -> Network:
    reaction_table = document["Reaction"]
    formulas = reaction_table.get("ReactionFormula")
    sides = [_parse_formula(f) for f in formulas]
    nr = len(formulas)

    compound_table = document.get("Compound")
    if compound_table is not None:
        metabolites = compound_table.get("ID")
    else:
        metabolites = sorted({m for s, p in sides for _, m in s + p})
    nm = len(metabolites)

    if reaction_table.has("ID"):
        reactions = reaction_table.get("ID")
    elif reaction_table.has("Reaction"):
        reactions = reaction_table.get("Reaction")
    else:
        reactions = [f"R{i + 1}" for i in range(nr)]

    N = np.zeros((nm, nr))
    for j, (substrates, products) in enumerate(sides):
        for sign, terms in ((-1.0, substrates), (1.0, products)):
            for coefficient, metabolite in terms:
                i = index_case_insensitive([metabolite], metabolites)[0]
                if i < 0:
                    raise ValueError(f'Unknown metabolite "{metabolite}" in reaction {reactions[j]}')
                N[i, j] = sign * coefficient

    regulation = np.zeros((nr, nm))
    if reaction_table.has("MetabolicRegulation"):
        for j, text in enumerate(reaction_table.get("MetabolicRegulation")):
            for sign, metabolite in _parse_regulation(text):
                i = index_case_insensitive([metabolite], metabolites)[0]
                if i >= 0:
                    regulation[j, i] = sign
    regulation_table = document.get("EnzymeRegulation")
    if regulation_table is not None:
        regulation = np.zeros((nr, nm))
        ind_m = index_case_insensitive(regulation_table.get("Compound"), metabolites)
        ind_r = index_case_insensitive(regulation_table.get("Reaction"), reactions)
        signs = {"activation": 1.0, "inhibition": -1.0}
        for i, j, kind in zip(ind_m, ind_r, regulation_table.get("RegulationType")):
            if i >= 0 and j >= 0 and kind in signs:
                regulation[j, i] = signs[kind]

    metabolite_kegg = [""] * nm
    metabolite_bigg = [""] * nm
    if compound_table is not None:
        metabolite_kegg = _column_or_empty(compound_table, "Identifiers:kegg.compound", nm)
        metabolite_bigg = _column_or_empty(compound_table, "Identifiers:bigg.compound", nm)
    reaction_kegg = _column_or_empty(reaction_table, "Identifiers:kegg.reaction", nr)
    reaction_bigg = _column_or_empty(reaction_table, "Identifiers:bigg.reaction", nr)

    # water is not treated as a reactant (Matlab: option omit_water in sbtab_to_network)
    for i in range(nm):
        name = metabolites[i].upper()
        if name in ("WATER", "H2O") or metabolite_kegg[i] in ("C00001", "kegg:C00001"):
            N[i, :] = 0.0

    return Network(metabolites, reactions, N, regulation, metabolite_kegg, reaction_kegg, metabolite_bigg, reaction_bigg)


# --------------------------------------------------------------------------------
# Log-normal conversions (Matlab: lognormal_normal2log, lognormal_log2normal; convention 'arithmetic')


def lognormal_normal2log(mu: np.ndarray, sigma: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    with np.errstate(invalid="ignore", divide="ignore"):
        mu = np.where(mu < LOG_NORMAL_MIN_MEAN, LOG_NORMAL_MIN_MEAN, mu)
        mu_log = np.log(mu) - 0.5 * np.log(1.0 + (sigma / mu) ** 2)
        mu_log = np.where(np.isnan(sigma), np.log(mu), mu_log)
        sigma_log = np.sqrt(np.log(1.0 + (sigma / mu) ** 2))
    return mu_log, sigma_log


# --------------------------------------------------------------------------------
# Kinetic data (Matlab: kinetic_data_load, kinetic_data_to_kinetics)


@dataclass
class KineticEntry:
    """Data for one type of quantity, in the shape (nr, 1), (nm, 1), or (nr, nm)."""

    median: np.ndarray
    mean: np.ndarray
    std: np.ndarray
    std_ln: np.ndarray


def _nanmean(values: np.ndarray) -> float:
    finite = values[~np.isnan(values)]
    return float(finite.mean()) if finite.size else math.nan


def _load_quantity(
    table: Table, quantity_name: str, prior_table: Dict[str, Quantity], network: Network
) -> KineticEntry:
    q = prior_table[quantity_name]
    shapes = {"Species": (network.nm, 1), "Reaction": (network.nr, 1), "Reaction/Species": (network.nr, network.nm)}
    shape = shapes[q.element]

    types = [QUANTITY_SYNONYMS.get(t, t) for t in table.get("QuantityType")]
    has_unit = table.has("Unit")
    if not has_unit:
        raise ValueError('Column "Unit" missing in kinetic data table')
    rows = [i for i, t in enumerate(types) if t == quantity_name and table.get("Unit")[i] == q.unit]

    def column(name: str) -> np.ndarray:
        if table.has(name):
            return to_floats([table.get(name)[i] for i in rows])
        return np.full(len(rows), math.nan)

    if table.has("Median"):
        median = column("Median")
    elif table.has("Value"):
        median = column("Value")
    else:
        median = column("Mode")

    mean = column("Mean")
    mean_is_given = table.has("Mean") and not np.isnan(mean).all()

    # values outside the allowed range are replaced by the bounds (Matlab: kinetic_data_load,
    # option enforce_ranges); NaN values stay NaN
    mean = np.clip(mean, q.lower, q.upper)
    median = np.clip(median, q.lower, q.upper)

    if q.scaling == "Multiplicative":
        zero = (mean == 0) | (median == 0)
        mean[zero] = math.nan
        median[zero] = math.nan

    std = column("Std")
    std[std == 0] = math.nan

    # complete missing standard deviations
    missing = (np.isfinite(mean).astype(int) + np.isfinite(median).astype(int) > 0) & ~np.isfinite(std)
    if q.scaling == "Additive":
        guess = q.data_std * np.ones_like(median)
    elif mean_is_given:
        guess = mean * math.sqrt(math.exp(math.log(q.data_geom_std) ** 2) - 1.0)
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            s = math.log(q.data_geom_std)
            guess = np.sqrt((math.exp(s**2) - 1.0) * np.exp(2 * np.log(median) + s**2))
    std[missing] = guess[missing]

    std_ln = np.full(len(rows), math.nan)
    if q.scaling == "Additive":
        if mean_is_given:
            median = mean.copy()
        else:
            mean = median.copy()
    elif mean_is_given:
        mean_ln, std_ln = lognormal_normal2log(mean, std)
        median = np.where(np.isfinite(mean_ln), np.exp(mean_ln), median)
    else:
        std_ln = math.log(q.data_geom_std) * np.ones(len(rows))

    # map rows to the network
    n = len(rows)
    selected = lambda name: [table.get(name)[i] for i in rows] if table.has(name) else None
    if q.element in ("Species", "Reaction/Species"):
        ids = selected("Compound")
        c_ind = index_case_insensitive(ids, network.metabolites) if ids is not None else None
    if q.element in ("Reaction", "Reaction/Species"):
        ids = selected("Reaction")
        r_ind = index_case_insensitive(ids, network.reactions) if ids is not None else None
    if q.element == "Species":
        first = c_ind if c_ind is not None else np.full(n, -1)
        second = np.zeros(n, dtype=int)
    elif q.element == "Reaction":
        first = r_ind if r_ind is not None else np.zeros(n, dtype=int)
        second = np.zeros(n, dtype=int)
    else:
        if r_ind is None or c_ind is None:
            first, second = np.zeros(n, dtype=int), np.zeros(n, dtype=int)
        else:
            first, second = r_ind, c_ind

    entry = KineticEntry(*(np.full(shape, math.nan) for _ in range(4)))
    # several values for the same quantity are averaged
    for key in sorted({(int(a), int(b)) for a, b in zip(first, second)}):
        if key[0] < 0 or key[1] < 0:
            continue
        hit = (first == key[0]) & (second == key[1])
        entry.median[key] = _nanmean(median[hit])
        entry.mean[key] = _nanmean(mean[hit])
        entry.std[key] = _nanmean(std[hit])
        entry.std_ln[key] = _nanmean(std_ln[hit]) if q.scaling == "Multiplicative" else math.nan

    if q.scaling == "Multiplicative" and mean_is_given:
        mean_ln, entry.std_ln = lognormal_normal2log(entry.mean, entry.std)
        entry.median = np.exp(mean_ln)
    return entry


def load_kinetic_data(document: Dict[str, Table], prior_table: Dict[str, Quantity], network: Network) -> Dict[str, KineticEntry]:
    """Read kinetic data from table ParameterData (or Parameter) and map them to the network."""
    table = document.get("ParameterData", document.get("Parameter"))
    if table is None:
        raise ValueError('No kinetic data table (ID "ParameterData" or "Parameter") found')

    data = {}
    for name in KINETIC_QUANTITIES:
        data[prior_table[name].symbol] = _load_quantity(table, name, prior_table, network)

    # Keq values are converted to standard concentration mM, if necessary
    standard_concentration = table.attributes.get("StandardConcentration", "undefined")
    if standard_concentration in ("M", "1M", "1 M"):
        mismatch = network.N.sum(axis=0).reshape(-1, 1)
        keq = data["Keq"]
        factor = 1000.0**mismatch
        keq.median = keq.median * factor
        keq.mean = keq.mean * factor
        keq.std = keq.std * factor
        data["dmu0"].median = data["dmu0"].median - RT * np.log(factor)
        data["dmu0"].mean = data["dmu0"].mean - RT * np.log(factor)
    elif standard_concentration not in ("mM", "1mM", "1 mM", "undefined"):
        raise ValueError(f"Unknown standard concentration {standard_concentration}")

    # missing Keq values are computed from standard Gibbs free energies of reaction
    keq, dmu0 = data["Keq"], data["dmu0"]
    missing = np.isnan(keq.median)
    keq.median[missing] = np.exp(-dmu0.median[missing] / RT)
    keq.std_ln[missing] = dmu0.std[missing] / RT
    return data


# --------------------------------------------------------------------------------
# State data (Matlab: sbtab_load_quantity_data, state_data_to_data)


STATE_TABLES = {
    "flux": ("FluxData", "rate of reaction", "reactions"),
    "metabolite": ("MetaboliteConcentrationData", "concentration", "metabolites"),
    "enzyme": ("EnzymeConcentrationData", "concentration of enzyme", "reactions"),
}


def _id_column(kind: str, match_data_by: str) -> str:
    element = "Compound" if kind == "metabolites" else "Reaction"
    suffix = "compound" if kind == "metabolites" else "reaction"
    return {
        "ModelElementId": element,
        "KeggId": f"{element}:Identifiers:kegg.{suffix}",
        "BiggId": f"{element}:Identifiers:bigg.{suffix}",
    }[match_data_by]


def _element_ids(network: Network, kind: str, match_data_by: str) -> List[str]:
    table = {
        ("metabolites", "ModelElementId"): network.metabolites,
        ("metabolites", "KeggId"): network.metabolite_kegg,
        ("metabolites", "BiggId"): network.metabolite_bigg,
        ("reactions", "ModelElementId"): network.reactions,
        ("reactions", "KeggId"): network.reaction_kegg,
        ("reactions", "BiggId"): network.reaction_bigg,
    }
    return table[(kind, match_data_by)]


def _read_state_table(
    table: Table, quantity_type: str, id_column: str, ids: Sequence[str], columns: Sequence[str]
) -> np.ndarray:
    """Data matrix (one row for each element in ids, in this order; NaN for missing values)."""
    use_rows = [i for i, t in enumerate(table.get("QuantityType")) if t == quantity_type]
    id_list = [table.get(id_column)[i] for i in use_rows]
    M = np.full((len(ids), len(columns)), math.nan)
    for k, element in enumerate(ids):
        if not element:
            continue
        if element in id_list:
            row = use_rows[id_list.index(element)]
            for c, name in enumerate(columns):
                if table.has(name):
                    M[k, c] = to_float(table.get(name)[row])
    return M


@dataclass
class StateData:
    samples: List[str]
    X_mean: np.ndarray
    X_std: np.ndarray
    lnE_mean: np.ndarray
    lnE_std: np.ndarray
    V_mean: np.ndarray
    V_std: np.ndarray


def _ln_mean_and_std(mean: np.ndarray, std: Optional[np.ndarray], default_geom_std: float) -> Tuple[np.ndarray, np.ndarray]:
    with np.errstate(divide="ignore", invalid="ignore"):
        if std is not None and np.isfinite(std).any():
            mean_ln, std_ln = lognormal_normal2log(mean, std)
            std_ln[~np.isfinite(std_ln)] = math.log(default_geom_std)
        else:
            mean_ln = np.log(mean)
            std_ln = math.log(default_geom_std) * np.ones_like(mean)
    return mean_ln, std_ln


def load_state_data(
    document: Dict[str, Table],
    network: Network,
    prior_table: Dict[str, Quantity],
    columns_mean: Sequence[str],
    columns_std: Sequence[str],
    match_data_by: str,
) -> StateData:
    ns = len(columns_mean)
    missing = [c for c in list(columns_mean) + list(columns_std)
               if not any(document[t].has(c) for t, _, _ in STATE_TABLES.values() if t in document)]
    if missing:
        raise ValueError(f"Columns not found in the state data tables: {missing}")
    matrices = {}
    for key, (table_id, quantity_type, kind) in STATE_TABLES.items():
        if table_id not in document:
            raise ValueError(f'State data table "{table_id}" not found')
        table = document[table_id]
        id_column = _id_column(kind, match_data_by)
        ids = _element_ids(network, kind, match_data_by)
        mean = _read_state_table(table, quantity_type, id_column, ids, columns_mean)
        std = _read_state_table(table, quantity_type, id_column, ids, columns_std) if columns_std else None
        matrices[key] = (mean, std)

    c_prior = prior_table["c"]
    u_prior = prior_table["u"]

    X, X_std = matrices["metabolite"]
    X = np.clip(X, c_prior.lower, c_prior.upper)  # NaN values stay NaN
    E, E_std = matrices["enzyme"]
    if (E == 0).any():
        raise ValueError("Enzyme concentration data contain zero values")

    X_mean, X_std = _ln_mean_and_std(X, X_std, c_prior.data_geom_std)
    lnE_mean, lnE_std = _ln_mean_and_std(E, E_std, u_prior.data_geom_std)

    V_mean, V_std = matrices["flux"]
    if V_std is None:
        V_std = np.full_like(V_mean, math.nan)
    if np.isnan(V_mean).any():
        raise ValueError("Flux data contain missing values")

    return StateData([f"S{i + 1}" for i in range(ns)], X_mean, X_std, lnE_mean, lnE_std, V_mean, V_std)


# --------------------------------------------------------------------------------
# Parameterisation (Matlab: mb_define_parameterisation)


def matlab_pinv(A: np.ndarray) -> np.ndarray:
    """Pseudoinverse with the cutoff of Matlab's pinv. With the default cutoff of numpy, the zero
    singular values of a rank-deficient covariance matrix can be inverted, giving huge numbers."""
    return np.linalg.pinv(A, rcond=max(A.shape) * np.finfo(float).eps)


def rref_pivots(A: np.ndarray) -> Tuple[np.ndarray, List[int]]:
    """Reduced row echelon form and pivot columns (same algorithm and tolerance as Matlab's rref)."""
    A = np.array(A, dtype=float)
    m, n = A.shape
    tol = max(m, n) * np.finfo(float).eps * np.linalg.norm(A, np.inf)
    pivots = []
    row = 0
    for col in range(n):
        if row >= m:
            break
        p = row + int(np.argmax(np.abs(A[row:, col])))
        if abs(A[p, col]) <= tol:
            A[row:, col] = 0.0
            continue
        pivots.append(col)
        A[[row, p], col:] = A[[p, row], col:]
        A[row, col:] /= A[row, col]
        for i in range(m):
            if i != row:
                A[i, col:] -= A[i, col] * A[row, col:]
        row += 1
    return A, pivots


@dataclass
class Parameterisation:
    """Independent parameters q = [Keq_ind, KV, KM, KA, KI] and all parameters
    qall = [Keq, KV, KM, KA, KI, Kcatf, Kcatr]; qall = M_q_to_qall * q."""

    nKeqind: int
    nKM: int
    nKA: int
    nKI: int
    M_q_to_qall: np.ndarray
    names: List[str]
    index: Dict[str, np.ndarray]


def define_parameterisation(network: Network) -> Parameterisation:
    nr = network.nr
    r_km, m_km = network.km_indices
    nKM = len(r_km)
    nKA = len(network.ka_indices[0])
    nKI = len(network.ki_indices[0])

    # Wegscheider conditions: Keq = Bt * Keq_ind
    echelon, pivots = rref_pivots(network.N)
    nKeqind = len(pivots)
    Bt = echelon[:nKeqind, :].T

    nq = nKeqind + nr + nKM + nKA + nKI
    ranges = np.cumsum([0, nr, nr, nKM, nKA, nKI, nr, nr])
    labels = ["Keq", "KV", "KM", "KA", "KI", "Kcatf", "Kcatr"]
    index = {k: np.arange(ranges[i], ranges[i + 1]) for i, k in enumerate(labels)}
    nqall = int(ranges[-1])
    names = [k for k in labels for _ in index[k]]

    qi = np.cumsum([0, nKeqind, nr, nKM, nKA, nKI])
    q_keq, q_kv, q_km, q_ka, q_ki = (np.arange(qi[i], qi[i + 1]) for i in range(5))

    N_km = np.zeros((nr, nKM))
    for k, (r, m) in enumerate(zip(r_km, m_km)):
        N_km[r, k] = network.N[m, r]

    M = np.zeros((nqall, nq))
    M[np.ix_(index["Keq"], q_keq)] = Bt
    M[np.ix_(index["KV"], q_kv)] = np.eye(nr)
    M[np.ix_(index["KM"], q_km)] = np.eye(nKM)
    M[np.ix_(index["KA"], q_ka)] = np.eye(nKA)
    M[np.ix_(index["KI"], q_ki)] = np.eye(nKI)
    M[np.ix_(index["Kcatf"], q_kv)] = np.eye(nr)
    M[np.ix_(index["Kcatf"], q_keq)] = 0.5 * Bt
    M[np.ix_(index["Kcatf"], q_km)] = -0.5 * N_km
    M[np.ix_(index["Kcatr"], q_kv)] = np.eye(nr)
    M[np.ix_(index["Kcatr"], q_keq)] = -0.5 * Bt
    M[np.ix_(index["Kcatr"], q_km)] = 0.5 * N_km

    return Parameterisation(nKeqind, nKM, nKA, nKI, M, names, index)


def collect_kinetic_vectors(
    kinetic: Dict[str, KineticEntry], network: Network, par: Parameterisation
) -> Tuple[np.ndarray, np.ndarray]:
    """Vectors of log data values and their standard deviations (Matlab: data.qall)."""
    n = len(par.names)
    mean, std = np.full(n, math.nan), np.full(n, math.nan)

    def fill(symbol: str, key: str, selector) -> None:
        with np.errstate(divide="ignore", invalid="ignore"):
            mean[par.index[key]] = np.log(selector(kinetic[symbol].median))
        std[par.index[key]] = selector(kinetic[symbol].std_ln)

    fill("Keq", "Keq", lambda x: x[:, 0])
    fill("KM", "KM", lambda x: x[network.km_indices])
    fill("KA", "KA", lambda x: x[network.ka_indices])
    fill("KI", "KI", lambda x: x[network.ki_indices])
    fill("Kcatf", "Kcatf", lambda x: x[:, 0])
    fill("Kcatr", "Kcatr", lambda x: x[:, 0])

    std[~np.isfinite(std)] = math.log(DATA_KIN_GEOM_STD)
    std[np.isnan(mean)] = math.nan
    return mean, std


# --------------------------------------------------------------------------------
# Prior + data (Matlab: mb_make_prior, mb_prepare_posterior)


def combine(prior_mean: np.ndarray, prior_std: np.ndarray, data_mean: np.ndarray, data_std: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Elementwise combination of prior and data; elements without valid data keep the prior."""
    mean, std = np.array(prior_mean, dtype=float), np.array(prior_std, dtype=float)
    ok = np.isfinite(data_mean) & np.isfinite(data_std)
    variance = 1.0 / (1.0 / data_std[ok] ** 2 + 1.0 / prior_std[ok] ** 2)
    mean[ok] = variance * (data_mean[ok] / data_std[ok] ** 2 + prior_mean[ok] / prior_std[ok] ** 2)
    std[ok] = np.sqrt(variance)
    return mean, std


def kinetic_preposterior(
    par: Parameterisation,
    prior_table: Dict[str, Quantity],
    qall_mean: np.ndarray,
    qall_std: np.ndarray,
    use_kinetic_data: str,
) -> Dict[str, np.ndarray]:
    nr = len(par.index["Keq"])
    parts = [
        ("Keq", par.nKeqind),
        ("KV", nr),
        ("KM", par.nKM),
        ("KA", par.nKA),
        ("KI", par.nKI),
    ]
    priors = [prior_of(prior_table, symbol) for symbol, _ in parts]
    prior_mean = np.concatenate([p.mean_ln * np.ones(n) for p, (_, n) in zip(priors, parts)])
    prior_std = np.concatenate([p.std_ln * np.ones(n) for p, (_, n) in zip(priors, parts)])
    prior_prec = np.diag(1.0 / prior_std**2)

    if use_kinetic_data == "all":
        select = np.arange(len(qall_mean))
    elif use_kinetic_data == "only_Keq":
        select = par.index["Keq"]
    elif use_kinetic_data == "none":
        select = np.array([], dtype=int)
    else:
        raise ValueError(f"Unknown option use_kinetic_data = {use_kinetic_data}")

    ok = select[np.isfinite(qall_mean[select]) & np.isfinite(qall_std[select])]
    q_prec = prior_prec.copy()
    right = prior_prec @ prior_mean
    if len(ok):
        A = par.M_q_to_qall[ok, :]
        data_prec = np.diag(1.0 / qall_std[ok] ** 2)
        q_prec = q_prec + A.T @ data_prec @ A
        right = right + A.T @ data_prec @ qall_mean[ok]
    q_mean = np.linalg.solve(q_prec, right)
    q_cov = np.linalg.inv(q_prec)

    qall_post_mean = par.M_q_to_qall @ q_mean
    qall_post_cov = par.M_q_to_qall @ q_cov @ par.M_q_to_qall.T

    result = {"mean": qall_post_mean, "cov": qall_post_cov, "prec": matlab_pinv(qall_post_cov)}
    for key in ["KM", "KA", "KI", "Kcatf", "Kcatr"]:
        i = par.index[key]
        result[f"prec_{key}"] = np.linalg.inv(qall_post_cov[np.ix_(i, i)]) if len(i) else np.zeros((0, 0))
    i = par.index["Keq"]
    result["prec_Keq"] = np.linalg.pinv(qall_post_cov[np.ix_(i, i)])  # singular if Wegscheider conditions exist
    return result


# --------------------------------------------------------------------------------
# Bounds (Matlab: mb_make_bounds)


def kinetic_bounds(par: Parameterisation, prior_table: Dict[str, Quantity]) -> Tuple[np.ndarray, np.ndarray]:
    """Lower and upper bounds for the vector of log kinetic constants qall."""
    lower, upper = np.zeros(len(par.names)), np.zeros(len(par.names))
    for key in ["Keq", "KV", "KM", "KA", "KI", "Kcatf", "Kcatr"]:
        lower[par.index[key]] = math.log(prior_table[key].lower)
        upper[par.index[key]] = math.log(prior_table[key].upper)
    return lower, upper


def metabolite_bounds(
    document: Dict[str, Table], network: Network, prior_table: Dict[str, Quantity]
) -> Tuple[np.ndarray, np.ndarray]:
    table = document.get("ConcentrationConstraint")
    if table is None:
        c = prior_table["c"]
        return math.log(c.lower) * np.ones(network.nm), math.log(c.upper) * np.ones(network.nm)

    # Unlike the Matlab code, values are matched by metabolite (Matlab uses the order of the rows)
    ids = _column_or_empty(table, "Compound", table.nrows)
    index = index_case_insensitive(network.metabolites, ids)
    if (index < 0).any():
        missing = [m for m, i in zip(network.metabolites, index) if i < 0]
        raise ValueError(f"No concentration constraint for metabolites: {missing}")
    c_min = to_floats(table.get("Min"))[index]
    c_max = to_floats(table.get("Max"))[index]
    same = c_min == c_max
    c_max[same] = c_max[same] * math.exp(CONC_EPSILON)
    return np.log(c_min), np.log(c_max)


# --------------------------------------------------------------------------------
# JSON problem structure (Matlab: cvxpy_problem_data_structure)


def clean(x):
    """Numpy arrays to nested lists, with None for non-finite numbers."""
    if isinstance(x, dict):
        return {k: clean(v) for k, v in x.items()}
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (float, np.floating)):
        return float(x) if math.isfinite(x) else None
    return x


def _kinetic_constant(unit, index, prec, post_mean, q_min, q_max, data_mean, data_std) -> dict:
    mean_ln = post_mean[index]
    with np.errstate(over="ignore"):
        return {
            "unit": unit,
            "true": np.full(len(index), math.nan),
            "bounds_ln": {"min": q_min[index], "max": q_max[index]},
            "data_ln": {"mean": data_mean[index], "std": data_std[index]},
            "bounds": {"min": np.exp(q_min[index]), "max": np.exp(q_max[index])},
            "combined": {"mean_ln": mean_ln, "prec_ln": prec, "geom_mean": np.exp(mean_ln)},
        }


def _state_quantity(prior_mean, prior_std, data_mean, data_std, ln_min, ln_max) -> dict:
    mean, std = combine(prior_mean, prior_std, data_mean, data_std)
    with np.errstate(divide="ignore"):
        return {
            "unit": "mM",
            "true": [],
            "bounds_ln": {"min": ln_min, "max": ln_max},
            "prior_ln": {"mean": prior_mean, "std": prior_std},
            "data_ln": {"mean": data_mean, "std": data_std},
            "bounds": {"min": np.exp(ln_min), "max": np.exp(ln_max)},
            "combined": {"geom_mean": np.exp(mean), "geom_std": np.exp(std)},
        }


def build_problem(
    network_file: str,
    kinetic_data_file: str,
    state_data_file: str,
    constraint_file: Optional[str] = None,
    columns_mean: Sequence[str] = ("Mean",),
    columns_std: Sequence[str] = (),
    match_data_by: str = "ModelElementId",
    use_kinetic_data: str = "all",
    prior_file: str = DEFAULT_PRIOR_FILE,
    km_max: float = KM_MAX,
) -> dict:
    """Read model balancing input files and return the problem as a dictionary (for JSON).

    km_max is the upper limit (in mM) for KM values: larger data values are replaced by it,
    and it is the upper bound for KM in the problem. It replaces the value from the prior table.
    """
    prior_table = read_prior_table(prior_file)
    prior_table["KM"].upper = km_max
    document = read_sbtab_files([network_file, kinetic_data_file, state_data_file])
    constraints = read_sbtab_files([constraint_file]) if constraint_file else {}

    network = build_network(document)
    par = define_parameterisation(network)
    kinetic = load_kinetic_data(document, prior_table, network)
    state = load_state_data(document, network, prior_table, columns_mean, columns_std, match_data_by)
    ns = len(state.samples)

    qall_mean, qall_std = collect_kinetic_vectors(kinetic, network, par)
    post = kinetic_preposterior(par, prior_table, qall_mean, qall_std, use_kinetic_data)
    q_min, q_max = kinetic_bounds(par, prior_table)

    # hide data that are not supposed to be used
    hide = np.full(len(qall_mean), use_kinetic_data != "all")
    if use_kinetic_data == "only_Keq":
        hide[par.index["Keq"]] = False
    shown_mean, shown_std = qall_mean.copy(), qall_std.copy()
    shown_mean[hide] = math.nan
    shown_std[hide] = math.nan

    constants = {}
    units = {"Keq": "depends on reaction stoichiometry", "Kcatf": "1/s", "Kcatr": "1/s", "KM": "mM", "KA": "mM", "KI": "mM"}
    for key, unit in units.items():
        constants[key] = _kinetic_constant(unit, par.index[key], post[f"prec_{key}"], post["mean"], q_min, q_max, shown_mean, shown_std)

    c = prior_of(prior_table, "c")
    u = prior_of(prior_table, "u")
    x_min, x_max = metabolite_bounds({**document, **constraints}, network, prior_table)

    X_shape = state.X_mean.shape
    metabolites = _state_quantity(
        c.mean_ln * np.ones(X_shape), c.std_ln * np.ones(X_shape), state.X_mean, state.X_std,
        np.tile(x_min[:, None], (1, ns)), np.tile(x_max[:, None], (1, ns)),
    )
    E_shape = state.lnE_mean.shape
    ln_e_min = np.full(E_shape, -math.inf)
    enzymes = _state_quantity(
        u.mean_ln * np.ones(E_shape), u.std_ln * np.ones(E_shape), state.lnE_mean, state.lnE_std,
        ln_e_min, math.log(prior_table["u"].upper) * np.ones(E_shape),
    )

    problem = {
        "standard_concentration": STANDARD_CONCENTRATION,
        "state_names": state.samples,
        "use_pseudo_values": False,
        "network": {
            "metabolite_names": network.metabolites,
            "reaction_names": network.reactions,
            "stoichiometric_matrix": network.N,
            "activation_matrix": (network.regulation > 0).astype(float).T,
            "inhibition_matrix": (network.regulation < 0).astype(float).T,
        },
        "kinetic_constants": {
            "all": {
                "names": par.names,
                "mean_ln": post["mean"],
                "prec_ln": post["prec"],
                "geom_mean": np.exp(post["mean"]),
            },
            **constants,
        },
        "metabolite_concentrations": metabolites,
        "enzyme_concentrations": enzymes,
        "reaction_fluxes": {
            "unit": "mM/s",
            "true": [],
            "data": {"mean": state.V_mean, "std": state.V_std},
        },
    }
    return clean(problem)


def fill_true_values(problem: dict) -> dict:
    """Replace the unknown true values by the geometric means of the preposterior.

    Real data have no true values, but io.read_arguments_from_json (Python package, October 2026)
    fails if they are missing. Remove this function once the reader accepts missing values.
    The true values are not used in the optimisation.
    """
    for constant in problem["kinetic_constants"].values():
        if "true" in constant and "combined" in constant:
            constant["true"] = constant["combined"]["geom_mean"]
    for key in ("metabolite_concentrations", "enzyme_concentrations"):
        problem[key]["true"] = problem[key]["combined"]["geom_mean"]
    return problem


def write_json(problem: dict, filename: str) -> None:
    with open(filename, "wt") as fp:
        json.dump(problem, fp, allow_nan=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert model balancing SBtab input files into a JSON problem file")
    parser.add_argument("network_file")
    parser.add_argument("kinetic_data_file")
    parser.add_argument("state_data_file")
    parser.add_argument("output_file")
    parser.add_argument("--constraints", default=None, help="SBtab file with table ConcentrationConstraint")
    parser.add_argument("--mean", nargs="+", default=["Mean"], help="columns with mean values (one for each state)")
    parser.add_argument("--std", nargs="*", default=[], help="columns with standard deviations")
    parser.add_argument("--match-data-by", default="ModelElementId", choices=["ModelElementId", "KeggId", "BiggId"])
    parser.add_argument("--use-kinetic-data", default="all", choices=["all", "only_Keq", "none"])
    parser.add_argument("--prior-file", default=DEFAULT_PRIOR_FILE)
    parser.add_argument("--km-max", type=float, default=KM_MAX,
                        help="upper limit for KM values in mM, for data and bounds (default: %(default)s)")
    parser.add_argument("--fill-true-values", action="store_true",
                        help="insert preposterior means as true values (needed by the current Python reader)")
    args = parser.parse_args()
    problem = build_problem(
        args.network_file, args.kinetic_data_file, args.state_data_file, args.constraints,
        args.mean, args.std, args.match_data_by, args.use_kinetic_data, args.prior_file, args.km_max,
    )
    if args.fill_true_values:
        fill_true_values(problem)
    write_json(problem, args.output_file)


if __name__ == "__main__":
    main()
