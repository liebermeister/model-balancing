#!/usr/bin/env python3
"""Generate the SBtab files of the MATLAB/Python comparison test model (fixed random seed).

The model contains what the published example models lack:
  activator (R1), inhibitor (R3), stoichiometry 2 (R2), several products (R3),
  a zero flux (R6 in S1), a flux reversal (R5 in S3), an always negative flux (R4),
  missing metabolite data (B in all states, A in S2) and missing enzyme data (R3 in S2, R6 in all states).
Run once; the generated files are kept in test_model/. Needs numpy only.
"""
import os
import numpy as np

rng = np.random.default_rng(1)
here = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_model")

mets = ["Xa", "Xb", "Xc", "Xd", "Xe", "Xf", "Xg", "A", "B", "C"]
external = {"Xa", "Xb", "Xc", "Xd", "Xe", "Xf", "Xg"}
rxns = ["R1", "R2", "R3", "R4", "R5", "R6"]
formulas = ["Xa <=> A", "2 A <=> B", "B <=> C + Xd", "Xg <=> C", "A <=> Xe", "B <=> Xf"]
regulation = {"R1": ("Xb", "activation"), "R3": ("Xc", "inhibition")}
states = ["S1", "S2", "S3"]

N = np.zeros((len(mets), len(rxns)))
for j, f in enumerate(formulas):
    left, right = f.split(" <=> ")
    for side, sign in ((left, -1), (right, 1)):
        for term in side.split(" + "):
            parts = term.split(" ")
            n, m = (1, parts[0]) if len(parts) == 1 else (int(parts[0]), parts[1])
            N[mets.index(m), j] = sign * n

# fluxes (mM/s); stationary in A, B, C
V = np.array([[1.0, 1.0, 0.5],   # R1
              [0.4, 0.3, 0.3],   # R2
              [0.4, 0.2, 0.3],   # R3
              [-0.4, -0.2, -0.3],  # R4 (always negative)
              [0.2, 0.4, -0.1],  # R5 (reverses in S3)
              [0.0, 0.1, 0.0]])  # R6 (zero in S1 and S3)
assert np.abs(N[[mets.index(m) for m in "ABC"]] @ V).max() < 1e-12

# true log concentrations (mM)
X = np.log(np.array([[2.0, 2.2, 1.8], [0.5, 0.6, 0.4], [0.3, 0.2, 0.5], [0.2, 0.25, 0.3],
                     [0.4, 0.3, 2.0], [0.05, 0.05, 0.08], [3.0, 3.5, 2.5],
                     [0.9, 0.8, 0.7], [0.3, 0.25, 0.28], [0.6, 0.7, 0.65]]))

# Keq: place ln Keq so that all flux directions are thermodynamically feasible (margin 0.5)
Q = N.T @ X
lnKeq = np.zeros(len(rxns))
for j in range(len(rxns)):
    lo = max([Q[j, s] for s in range(3) if V[j, s] > 0], default=-np.inf)
    hi = min([Q[j, s] for s in range(3) if V[j, s] < 0], default=np.inf)
    assert lo < hi, rxns[j]
    lnKeq[j] = (lo + hi) / 2 if np.isfinite(lo + hi) else (lo + 1.0 if np.isfinite(lo) else hi - 1.0)

km_pairs = [(i, j) for j in range(len(rxns)) for i in range(len(mets)) if N[i, j] != 0]
lnKM = np.log(0.5) + 0.8 * rng.standard_normal(len(km_pairs))
KM = np.full(N.shape, np.nan)
for (i, j), v in zip(km_pairs, lnKM):
    KM[i, j] = v
lnKV = np.log(100.0) + 0.7 * rng.standard_normal(len(rxns))
sumNlnKM = np.array([sum(N[i, j] * KM[i, j] for i in range(len(mets)) if N[i, j] != 0) for j in range(len(rxns))])
lnKcatf = lnKV + 0.5 * (lnKeq - sumNlnKM)
lnKcatr = lnKV - 0.5 * (lnKeq - sumNlnKM)
lnKA, lnKI = np.log(0.4), np.log(0.5)
actr, inhr = mets.index("Xb"), mets.index("Xc")

# enzyme levels from the common modular rate law
E = np.zeros(V.shape)
for j in range(len(rxns)):
    for s in range(3):
        if V[j, s] == 0:
            continue
        ls = [i for i in range(len(mets)) if N[i, j] < 0]
        lp = [i for i in range(len(mets)) if N[i, j] > 0]
        u = lambda l: np.exp(sum(abs(N[i, j]) * (X[i, s] - KM[i, j]) for i in l))
        p = lambda l: np.prod([(1 + np.exp(X[i, s] - KM[i, j])) ** abs(N[i, j]) for i in l])
        D = p(ls) + p(lp) - 1
        w = (np.exp(lnKcatf[j]) * u(ls) - np.exp(lnKcatr[j]) * u(lp)) / D
        if j == 0:
            c = np.exp(X[actr, s]); w *= c / (c + np.exp(lnKA))
        if j == 2:
            c = np.exp(X[inhr, s]); w *= 1 / (1 + c / np.exp(lnKI))
        E[j, s] = V[j, s] / w
assert (E[V != 0] > 0).all()

# noisy data, with missing values
Cdata = np.exp(X + 0.1 * rng.standard_normal(X.shape))
Edata = np.exp(np.log(np.where(E > 0, E, 1.0)) + 0.1 * rng.standard_normal(E.shape))
Cdata[mets.index("B"), :] = np.nan
Cdata[mets.index("A"), 1] = np.nan
Edata[2, 1] = np.nan
Edata[5, :] = np.nan
Cstd, Estd = 0.2 * Cdata, 0.2 * Edata
Vstd = 0.1 * np.abs(V) + 0.01


def num(x):
    return "NaN" if not np.isfinite(x) else f"{x:.10g}"


def write(name, text):
    with open(os.path.join(here, name), "wt") as fp:
        fp.write(text)


lines = ["!!!SBtab SBtabVersion='1.1' Document='test_model'", "",
         "!!SBtab TableID='Reaction' TableType='Reaction' TableName='Reaction'",
         "!ID\t!ReactionFormula\t!IsReversible\t!MetabolicRegulation"]
for r, f in zip(rxns, formulas):
    sign = {"activation": "+", "inhibition": "-"}
    reg = f"{sign[regulation[r][1]]} {regulation[r][0]}" if r in regulation else ""
    lines.append(f"{r}\t{f}\tTrue\t{reg}")
lines += ["", "!!SBtab TableID='Compound' TableType='Compound' TableName='Compound'", "!ID\t!Name\t!IsConstant"]
lines += [f"{m}\t{m}\t{'True' if m in external else 'False'}" for m in mets]
write("test_model_network.tsv", "\n".join(lines) + "\n")

lines = ["!!SBtab TableID='ParameterData' TableType='Quantity'", "!QuantityType\t!Compound\t!Reaction\t!Unit\t!Median"]
# only part of the constants are measured
for k, (i, j) in enumerate(km_pairs):
    if k % 2 == 0:
        lines.append(f"Michaelis constant\t{mets[i]}\t{rxns[j]}\tmM\t{num(np.exp(KM[i, j]))}")
lines.append(f"activation constant\tXb\tR1\tmM\t{num(np.exp(lnKA))}")
lines.append(f"inhibitory constant\tXc\tR3\tmM\t{num(np.exp(lnKI))}")
for j, r in enumerate(rxns):
    lines.append(f"equilibrium constant\t\t{r}\tdimensionless\t{num(np.exp(lnKeq[j]))}")
    if j != 3:
        lines.append(f"substrate catalytic rate constant\t\t{r}\t1/s\t{num(np.exp(lnKcatf[j]))}")
    if j % 2 == 0:
        lines.append(f"product catalytic rate constant\t\t{r}\t1/s\t{num(np.exp(lnKcatr[j]))}")
write("test_model_kinetic_data.tsv", "\n".join(lines) + "\n")

cols = "\t".join(f"{s}_Mean\t{s}_Std" for s in states)


def block(name, tid, unit, idcol, qt, ids, M, S):
    out = [f"!!SBtab TableName='{name}' TableID='{tid}' TableType='QuantityMatrix' Unit='{unit}'",
           f"!QuantityType\t!{idcol}\t{cols}"]
    for k, i in enumerate(ids):
        out.append(f"{qt}\t{i}\t" + "\t".join(f"{num(M[k, s])}\t{num(S[k, s])}" for s in range(3)))
    return out + [""]


lines = ["!!!SBtab SBtabVersion='1.1'", ""]
lines += block("Metabolic fluxes", "FluxData", "mM/s", "Reaction", "rate of reaction", rxns, V, Vstd)
lines += block("Metabolite concentrations", "MetaboliteConcentrationData", "mM", "Compound", "concentration", mets, Cdata, Cstd)
lines += block("Enzyme concentrations", "EnzymeConcentrationData", "mM", "Reaction", "concentration of enzyme", rxns, Edata, Estd)
write("test_model_state_data.tsv", "\n".join(lines))

lines = ["!!!SBtab SBtabVersion='1.1' Document='test_model'", "",
         "!!SBtab TableID='ConcentrationConstraint' TableType='Quantity' Unit='mM'", "!QuantityType\t!Compound\t!Min\t!Max"]
lines += [f"concentration\t{m}\t0.000001\t1000" for m in mets]
write("test_model_constraints.tsv", "\n".join(lines) + "\n")
