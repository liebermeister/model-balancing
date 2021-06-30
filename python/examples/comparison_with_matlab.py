import os
import sys
import warnings
import pandas as pd

import cvxpy as cp
from model_balancing import ALL_VARIABLES
from model_balancing.io import read_arguments_json
from model_balancing.model_balancing_cvx import ModelBalancingConvex, NegativeFluxError
from model_balancing.model_balancing_noncvx import ModelBalancing
from path import Path

os.chdir("/home/eladn/git/model-balancing")

json_fnames = [f"examples/JSON/branch_point_model_elad.json"]
# json_fnames = Path("examples/JSON").listdir("*.json")

INITIALIZE_WITH_CONVEX = True

z_scores_data = []

for json_fname in json_fnames:
    if json_fname.endswith(".json"):
        example_name = os.path.basename(json_fname.replace(".json", ""))
    else:
        continue

    print(f"Analyzing example: {example_name}")

    args = read_arguments_json(json_fname)
    
    initial_point = {}
    if INITIALIZE_WITH_CONVEX:
        try:
            mbc = ModelBalancingConvex(**args)
            mbc.initialize_with_gmeans()
            if not mbc.is_gmean_feasible():
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
                f"ln_{p}": mbc._var_dict[f"ln_{p}"].value
                for p in ALL_VARIABLES
                if mbc._var_dict[f"ln_{p}"] is not None
            }
        except NegativeFluxError as e:
            print(str(e))
        except cp.error.SolverError as e:
            print(str(e))

        result_dict = mbc.get_z_scores()
        result_dict["JSON"] = example_name
        result_dict["CVXPY"] = True
        result_dict["alpha"] = 0.0
        result_dict["beta"] = 0.0
        result_dict["objective"] = mbc.objective_value
        z_scores_data.append(result_dict)

    mb = ModelBalancing(**args)

    # for a in [0.0, 0.001, 0.01, 0.1, 0.5, 1.0]:
    for a in [0.0, 1.0]:
        for b in [0.0, 1.0]:
            warnings.filterwarnings("ignore", category=RuntimeWarning)

            mb.alpha = a
            mb.beta = b

            # initialize solver with the Convex optimization solution
            mb._var_dict.update(initial_point)
            print(
                f"Solving using non-convex solver, α = {a:5.1g}, β = {b:5.1g} ... ",
                end="",
            )
            mb.solve(solver="SLSQP", options={"maxiter": 1000, "disp": False})
            print(f"optimized total squared Z-scores = {mb.objective_value:.3f}")

            result_dict = mb.get_z_scores()
            result_dict["JSON"] = example_name
            result_dict["CVXPY"] = False
            result_dict["alpha"] = a
            result_dict["beta"] = b
            result_dict["objective"] = mb.objective_value
            z_scores_data.append(result_dict)

            with open(f"python/res/{example_name}_alpha_{a:.1g}_state.tsv", "wt") as fp:
                fp.write(mb.to_state_sbtab().to_str())
            with open(f"python/res/{example_name}_alpha_{a:.1g}_model.tsv", "wt") as fp:
                fp.write(mb.to_model_sbtab().to_str())

df = pd.DataFrame.from_dict(z_scores_data).set_index(["JSON", "alpha", "beta"])
df.round(5).to_csv("python/res/z_score_report.csv")
