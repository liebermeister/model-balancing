## Solve Model balancing problem with kinetic parameters as variables

#%%

import numpy as np
from model_balancing import ModelBalancing
from util import Q_

#%%

use_regulation = True
args = {}

# choose solver

args["solver"] = "SLSQP"

#%%

# Network topology

args["S"] = np.matrix(
    [
        [-1, 0, 0],
        [1, -1, -1],
        [0, 1, 1]
    ], dtype=float
)

if use_regulation:
    args["A_act"] = np.matrix(
        [
            [0, 0, 1],
            [0, 0, 0],
            [0, 0, 0]
        ], dtype=float
    )
    args["A_inh"] = np.matrix(
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ], dtype=float
    )
else:
    args["A_act"] = np.matrix(
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ], dtype=float
    )
    args["A_inh"] = np.matrix(
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ], dtype=float
    )

#%%

# Priors for Kinetic and Thermodynamic parameters

args["Keq_gmean"] = Q_([1, 1, 1], "")  # geometric mean (assuming a standard concentration of 1M)
args["Keq_ln_cov"] = np.array(
    [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]
)  # log-scale covariance

args["kcatf_gmean"] = Q_([20.0, 8.0, 10.0], "1/s")  # geometric mean
args["kcatf_ln_cov"] = np.array(
    [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]
) # log-scale covariance

args["kcatr_gmean"] = Q_([3.5e-3, 6.2e4, 3.1e-1], "1/s")  # geometric mean
args["kcatr_ln_cov"] = np.array(
    [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]
) # log-scale covariance

# the K parameters are only for existing connections in the S, A_act and A_inh matrices.
# the order of values is "metabolite-first".
args["Km_gmean"] = Q_([1e-2, 1e-4, 1e-4, 1e-3, 1e-1, 1e-1], "M")
args["Km_ln_cov"] = np.array(
    [
        [1, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1],
    ]
) # log-scale covariance

args["Ka_gmean"] = Q_([1e-3, 5e-3], "M")
args["Ka_ln_cov"] = np.array(
    [
        [1, 0],
        [0, 1],
    ]
)

args["Ki_gmean"] = Q_([1e-3], "M")
args["Ki_ln_cov"] = np.array(
    [
        [1],
    ]
)

#%%

# condition dependent data (columns represent conditions)

args["fluxes"] = Q_(
    [
        [2.0, 1.5, 2.5, 2.5],
        [1.0, 1.0, 2.0, 0.5],
        [1.0, 0.5, 0.5, 2.0],
    ], "mM/s"
)

args["conc_enz_gmean"] = Q_(
    [
        [1e-3, 1e-3, 1e-3, 3e-3],
        [2e-3, 3e-3, 1e-3, 2e-3],
        [1e-3, 2e-3, 2e-3, 1e-3],
    ], "M")  # geometric mean

args["conc_enz_gstd"] = np.array(
    [
        [1.1, 1.1, 1.1, 1.1],
        [1.1, 1.1, 1.1, 1.1],
        [1.1, 1.1, 1.1, 1.1],
    ]
)  # geometric standard deviation

args["conc_met_gmean"] = Q_(
    [
        [3e-3, 4e-3, 9e-4, 3e-4],
        [2e-3, 3e-3, 6e-4, 2e-4],
        [1e-3, 2e-3, 3e-4, 1e-4],
    ], "M"
)  # geometric mean

args["conc_met_gstd"] = np.array(
    [
        [1.1, 1.1, 1.1, 1.1],
        [1.1, 1.1, 1.1, 1.1],
        [1.1, 1.1, 1.1, 1.1],
    ]
)  # geometric standard deviation

#%%

args["rate_law"] = "CM"

mb = ModelBalancing(**args)

# check that the driving forces at the mode are all positive
assert mb.is_gmean_feasible()

mb.initialize_with_gmeans()
print(f"initial total squared Z-scores = {mb.objective_value}")
mb.print_z_scores(5)
mb.print_status()

constr = mb._thermodynamic_constraints()
print(constr.A)

#%%
mb.solve()
print(f"optimized total squared Z-scores = {mb.objective_value}")

print("-"*25 + " All Z-scores " + "-"*25)
mb.print_z_scores(5)

mb.print_status()
