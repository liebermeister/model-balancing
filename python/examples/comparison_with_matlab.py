import os
from model_balancing import (
    ModelBalancingConvex,
)
from model_balancing.model_balancing import ModelBalancing
import cvxpy as cp
from model_balancing.io import read_arguments_json
import numpy as np

os.chdir("/home/eladn/git/model-balancing")

json_fnames = []
# json_fnames += ["three_chain_model_artificial_noisy_state_noisy_kinetic.json"]
# json_fnames += ["e_coli_artificial_noisefree_state_noisefree_Keq.json"]
# json_fnames += ["e_coli_noor_2016_no_kinetic_data_balanced.json"]
json_fnames += ["branch_point_model_artificial_noisefree_state_noisefree_kinetic.json"]
#json_fnames += ["branch_point_model_artificial_noisy_state_noisy_kinetic.json"]

# json_fnames = os.listdir("cvxpy/examples/JSON")

json_fnames = [f"examples/JSON/{n}" for n in json_fnames]

for json_fname in json_fnames:
    if json_fname.endswith(".json"):
        example_name = os.path.basename(json_fname.replace(".json", ""))
    else:
        continue

    print(f"Analyzing example: {example_name}")

    args = read_arguments_json(json_fname)
    initial_point = {f"ln_{p}": None for p in ModelBalancing.INDEPENDENT_VARIABLES}

    if all(args["fluxes"].flatten() > 0):
        mbc = ModelBalancingConvex(**args)
        # if not mbc.is_gmean_feasible():
        #    print("geometric mean is not a feasible solution")
        #    continue

        mbc.initialize_with_gmeans()
        if mbc.find_inner_point(verbose=False) != cp.OPTIMAL:
            print("Cannot find an inner point given the constraints")
            continue

        mbc.solve(verbose=False)
        print(
            f"Convex optimization (equivalent to α = 0) ... optimized total "
            f"squared Z-scores = {mbc.objective_value:.3f}"
        )
        with open(f"python/res/{example_name}_state_convex.tsv", "wt") as fp:
            fp.write(mbc.to_state_sbtab().to_str())

        with open(f"python/res/{example_name}_model_convex.tsv", "wt") as fp:
            fp.write(mbc.to_model_sbtab().to_str())

        for p in ModelBalancing.INDEPENDENT_VARIABLES:
            if mbc.__getattribute__(f"ln_{p}").size != 0:
                initial_point[f"ln_{p}"] = mbc.__getattribute__(f"ln_{p}").value

    mb = ModelBalancing(**args)

    for a in [0.0, 0.001, 0.01, 0.1, 0.5, 1.0]:
        # warnings.filterwarnings("ignore", category=RuntimeWarning)

        # initialize solver with the Convex optimization solution
        for p in ModelBalancing.INDEPENDENT_VARIABLES:
            if initial_point[f"ln_{p}"] is not None:
                mb.__setattr__(f"ln_{p}", initial_point[f"ln_{p}"])

        print(f"Solving using non-convex solver, α = {a:5.1g} ... ", end="")
        mb.alpha = a
        mb.solve(solver="SLSQP", options={"maxiter": 1000, "disp": False})
        mb.solve()
        print(f"optimized total squared Z-scores = {mb.objective_value:.3f}")
        with open(f"python/res/{example_name}_state_alpha_{a:.1g}.tsv", "wt") as fp:
            fp.write(mb.to_state_sbtab().to_str())
        with open(f"python/res/{example_name}_model_alpha_{a:.1g}.tsv", "wt") as fp:
            fp.write(mb.to_model_sbtab().to_str())
