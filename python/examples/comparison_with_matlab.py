import os
import warnings
from typing import Union, List
import json
import numpy as np
from model_balancing import ModelBalancing, ModelBalancingConvex
import cvxpy as cp
import pandas as pd


#%%
os.chdir("/home/eladn/git/model-balancing")

#json_fnames = ["three_chain_model_artificial_noisy_state_noisy_kinetic.json"]
#json_fnames = ["e_coli_artificial_noisefree_state_noisefree_Keq.json"]
#json_fnames = ["double_branch_model_artificial_noisefree_state_noisy_Keq.json"]
json_fnames = ["branch_point_model_artificial_noisefree_state_noisefree_Keq.json"]
# json_fnames = os.listdir("cvxpy/examples/JSON")

for json_fname in json_fnames:
    if json_fname.endswith(".json"):
        example_name = json_fname.replace(".json", "")
    else:
        continue

    print(f"Analyzing example: {example_name}")

    mbc = ModelBalancingConvex.from_json(json_fname)
    # if not mbc.is_gmean_feasible():
    #    print("geometric mean is not a feasible solution")
    #    continue
    mbc.solve(verbose=False)
    print(
        f"Convex optimization (equivalent to α = 0) ... optimized total squared Z-scores = {mbc.objective_value:.3f}"
    )
    mbc.to_state_sbtab().write(f"res/{example_name}_state_convex.tsv")
    mbc.to_model_sbtab().write(f"res/{example_name}_model_convex.tsv")

    mb = ModelBalancing.from_json(json_fname)
    for a in [0.0, 0.001, 0.01, 0.1, 0.5, 1.0]:
        # warnings.filterwarnings("ignore", category=RuntimeWarning)

        # initialize solver with the Convex optimization solution
        for p in ModelBalancing.INDEPENDENT_VARIABLES:
            if mbc.__getattribute__(f"ln_{p}").size != 0:
                mb.__setattr__(f"ln_{p}", mbc.__getattribute__(f"ln_{p}").value)

        print(f"Solving using non-convex solver, α = {a:5.1g} ... ", end="")
        mb.alpha = a
        mb.solve()
        print(f"optimized total squared Z-scores = {mb.objective_value:.3f}")
        mb.to_state_sbtab().write(f"res/{example_name}_state_alpha_{a:.1g}.tsv")
        mb.to_model_sbtab().write(f"res/{example_name}_model_alpha_{a:.1g}.tsv")
