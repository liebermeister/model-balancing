import os
import warnings
import pandas as pd
import time
import datetime
import tarfile
import sys

import cvxpy as cp
from model_balancing import ALL_VARIABLES
from model_balancing.io import read_arguments_json
from model_balancing.model_balancing_cvx import ModelBalancingConvex, NegativeFluxError
from model_balancing.model_balancing_noncvx import ModelBalancing
from path import Path

# change to the root path of the model-balancing package
os.chdir(Path(__file__).parent.parent.parent)

timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

json_fnames = Path("examples/JSON").listdir("*.json")


def get_convex_solution(args: dict) -> dict:
    mbc = ModelBalancingConvex(**args)
    mbc.initialize_with_gmeans()
    if not mbc.is_gmean_feasible():
        if mbc.find_inner_point(verbose=False) not in \
                cp.settings.SOLUTION_PRESENT:
            print("Cannot find an inner point given the constraints")
            return None

    if mbc.solve(verbose=False) in cp.settings.SOLUTION_PRESENT:
        return {
            f"ln_{p}": mbc._var_dict[f"ln_{p}"].value
            for p in ALL_VARIABLES
            if mbc._var_dict[f"ln_{p}"] is not None
        }

z_scores_data = []
with tarfile.open(
    f"python/res/model_balancing_{timestamp}.tar.gz", "w:gz"
) as result_tarfile:
    for json_fname in json_fnames:
        if json_fname.endswith(".json"):
            example_name = os.path.basename(json_fname.replace(".json", ""))
        else:
            continue

        print(f"Analyzing example: {example_name}")

        args = read_arguments_json(json_fname)
        args["beta"] = 0.0

        for init_with_mbc in [True, False]:
            for a in [0.0, 1.0]:
                result_dict = {
                    "JSON": example_name,
                    "CVXPY": False,
                    "alpha": a,
                }
                args["alpha"] = a
                mb = ModelBalancing(**args)
                if init_with_mbc:
                    initial_point = get_convex_solution(args)
                    if initial_point is None:
                        continue
                    mb._var_dict.update(initial_point)
                    result_dict["initialization"] = "convex solution"
                    print(
                        f"Initializing non-convex solver with CVXPY solution, "
                        f"α = {a:5.1g} ... ",
                        end="",
                    )
                else:
                    mb.initialize_with_gmeans()
                    result_dict["initialization"] = "geometric means"
                    print(
                        f"Initializing non-convex solver with geometric means, "
                        f"α = {a:5.1g} ... ",
                        end="",
                    )
                tic = time.perf_counter()
                status = mb.solve(
                    solver="SLSQP", options={"maxiter": 1000, "disp": False}
                )
                toc = time.perf_counter()
                result_dict["runtime"] = toc - tic
                result_dict["status"] = status
                if status != "optimal":
                    print(f"optimization failed: status = {status}")
                    z_scores_data.append(result_dict)
                else:
                    print(
                        f"optimized total squared Z-scores = {mb.objective_value:.3f}"
                    )
                    result_dict.update(mb.get_z_scores())
                    result_dict["objective"] = mb.objective_value
                    z_scores_data.append(result_dict)

                    state_fname = f"/tmp/{example_name}_a{a:.1g}_state.tsv"
                    model_fname = f"/tmp/{example_name}_a{a:.1g}_model.tsv"
                    with open(state_fname, "wt") as fp:
                        fp.write(mb.to_state_sbtab().to_str())
                    result_tarfile.add(state_fname)
                    os.remove(state_fname)
                    with open(model_fname, "wt") as fp:
                        fp.write(mb.to_model_sbtab().to_str())
                    result_tarfile.add(model_fname)
                    os.remove(model_fname)

    df = pd.DataFrame.from_dict(z_scores_data).set_index(["JSON", "alpha"])
    summary_fname = f"/tmp/summary.csv"
    df.round(5).to_csv(summary_fname)
    result_tarfile.add(summary_fname)
