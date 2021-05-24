import os
from path import Path
from model_balancing import (
    ModelBalancing,
    ModelBalancingConvex,
    INDEPENDENT_VARIABLES,
)
import cvxpy as cp
from model_balancing.io import read_arguments_json
import numpy as np

os.chdir("/home/eladn/git/model-balancing")

# json_fnames = []
# json_fnames += ["double_branch_model_artificial_noisy_state_no_kinetic.json"]
# json_fnames = [f"examples/JSON/{n}" for n in json_fnames]

json_fnames = Path("examples/JSON").listdir("*.json")


for json_fname in json_fnames:
    if json_fname.endswith(".json"):
        example_name = os.path.basename(json_fname.replace(".json", ""))
    else:
        continue

    print(f"Analyzing example: {example_name}")

    args = read_arguments_json(json_fname)

    try:
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
        with open(f"python/res/{example_name}_convex_state.tsv", "wt") as fp:
            fp.write(mbc.to_state_sbtab().to_str())

        with open(f"python/res/{example_name}_convex_model.tsv", "wt") as fp:
            fp.write(mbc.to_model_sbtab().to_str())

        initial_point = {
            f"ln_{p}": mbc.__getattribute__(f"ln_{p}").value
            for p in INDEPENDENT_VARIABLES
            if mbc.__getattribute__(f"ln_{p}").size != 0
        }
    except ValueError as e:
        print(str(e))
        initial_point = {}
    except cp.error.SolverError as e:
        print(str(e))
        initial_point = {}

    mb = ModelBalancing(**args)

    for a in [0.0, 0.001, 0.01, 0.1, 0.5, 1.0]:
        # warnings.filterwarnings("ignore", category=RuntimeWarning)

        # initialize solver with the Convex optimization solution
        for k, v in initial_point.items():
            mb.__setattr__(k, v)

        print(f"Solving using non-convex solver, α = {a:5.1g} ... ", end="")
        mb.alpha = a
        mb.solve(solver="SLSQP", options={"maxiter": 1000, "disp": False})
        mb.solve()
        print(f"optimized total squared Z-scores = {mb.objective_value:.3f}")
        with open(f"python/res/{example_name}_alpha_{a:.1g}_state.tsv", "wt") as fp:
            fp.write(mb.to_state_sbtab().to_str())
        with open(f"python/res/{example_name}_alpha_{a:.1g}_model.tsv", "wt") as fp:
            fp.write(mb.to_model_sbtab().to_str())
