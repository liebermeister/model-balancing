import os
import warnings
import pandas as pd
import time
import datetime
import tarfile
import tempfile

import cvxpy as cp
from model_balancing import ALL_VARIABLES
from model_balancing.io import read_arguments_json
from model_balancing.model_balancing_cvx import ModelBalancingConvex, NegativeFluxError
from model_balancing.model_balancing_noncvx import ModelBalancing
from path import Path

os.chdir("/home/eladn/git/model-balancing")

timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

# json_fnames = [f"examples/JSON/branch_point_model_elad.json"]
json_fnames = Path("examples/JSON").listdir("*.json")

INITIALIZE_WITH_CONVEX = True

z_scores_data = []
with tarfile.open(f"python/res/z_score_{timestamp}.tar.gz", "w:gz") as result_tarfile:
    for json_fname in json_fnames:
        if json_fname.endswith(".json"):
            example_name = os.path.basename(json_fname.replace(".json", ""))
        else:
            continue

        print(f"Analyzing example: {example_name}")

        args = read_arguments_json(json_fname)

        for b in [0.0, 1e-2]:
            initial_point = None
            if INITIALIZE_WITH_CONVEX:
                try:
                    args["beta"] = b
                    mbc = ModelBalancingConvex(**args)
                    mbc.initialize_with_gmeans()
                    if not mbc.is_gmean_feasible():
                        if mbc.find_inner_point(verbose=False) != cp.OPTIMAL:
                            print("Cannot find an inner point given the constraints")
                            continue

                    print(
                        f"Solving using convex solver (similar to α = 0), "
                        f"β = {b:5.1g} ... ",
                        end="",
                    )
                    tic = time.perf_counter()
                    status = mbc.solve(verbose=False)
                    toc = time.perf_counter()

                    result_dict = {}
                    result_dict["JSON"] = example_name
                    result_dict["CVXPY"] = True
                    result_dict["status"] = status
                    result_dict["runtime"] = toc - tic
                    result_dict["alpha"] = None
                    result_dict["beta"] = b
                    if status not in cp.settings.SOLUTION_PRESENT:
                        print(f"CVXPY failed: {status}")
                    else:
                        print(
                            f"optimized total squared Z-scores = {mbc.objective_value:.3f}"
                        )

                        result_dict.update(mbc.get_z_scores())
                        result_dict["objective"] = mbc.objective_value
                        z_scores_data.append(result_dict)

                        state_fname = f"/tmp/{example_name}_convex_b{b:.1g}_state.tsv"
                        model_fname = f"/tmp/{example_name}_convex_b{b:.1g}_model.tsv"
                        with open(state_fname, "wt") as fp:
                            fp.write(mbc.to_state_sbtab().to_str())
                        result_tarfile.add(state_fname)
                        os.remove(state_fname)
                        with open(model_fname, "wt") as fp:
                            fp.write(mbc.to_model_sbtab().to_str())
                        result_tarfile.add(model_fname)
                        os.remove(model_fname)

                        initial_point = {
                            f"ln_{p}": mbc._var_dict[f"ln_{p}"].value
                            for p in ALL_VARIABLES
                            if mbc._var_dict[f"ln_{p}"] is not None
                        }
                except NegativeFluxError as e:
                    print(str(e))
                except cp.error.SolverError as e:
                    print(str(e))

            mb = ModelBalancing(**args)

            for a in [0.0, 0.001, 0.01, 0.1, 0.5, 1.0]:
                warnings.filterwarnings("ignore", category=RuntimeWarning)

                mb.alpha = a

                # initialize solver with the Convex optimization solution
                if initial_point is not None:
                    mb._var_dict.update(initial_point)
                else:
                    mb.initialize_with_gmeans()

                print(
                    f"Solving using non-convex solver, α = {a:5.1g}, β = {b:5.1g} ... ",
                    end="",
                )
                tic = time.perf_counter()
                status = mb.solve(
                    solver="SLSQP", options={"maxiter": 1000, "disp": False}
                )
                toc = time.perf_counter()
                print(f"optimized total squared Z-scores = {mb.objective_value:.3f}")

                result_dict = mb.get_z_scores()
                result_dict["JSON"] = example_name
                result_dict["CVXPY"] = False
                result_dict["status"] = status
                result_dict["runtime"] = toc - tic
                result_dict["alpha"] = a
                result_dict["beta"] = b
                result_dict["objective"] = mb.objective_value
                z_scores_data.append(result_dict)

                state_fname = f"/tmp/{example_name}_a{a:.1g}_b{b:.1g}_state.tsv"
                model_fname = f"/tmp/{example_name}_a{a:.1g}_b{b:.1g}_model.tsv"
                with open(state_fname, "wt") as fp:
                    fp.write(mb.to_state_sbtab().to_str())
                result_tarfile.add(state_fname)
                os.remove(state_fname)
                with open(model_fname, "wt") as fp:
                    fp.write(mb.to_model_sbtab().to_str())
                result_tarfile.add(model_fname)
                os.remove(model_fname)

    df = pd.DataFrame.from_dict(z_scores_data).set_index(["JSON", "alpha", "beta"])
    summary_fname = f"/tmp/summary.csv"
    df.round(5).to_csv(summary_fname)
    result_tarfile.add(summary_fname)
