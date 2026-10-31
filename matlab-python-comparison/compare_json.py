#!/usr/bin/env python3
"""Compare JSON files from the Matlab exporter with those from python-sbtab-to-json.

usage: compare_json.py RESULT_DIR
Differences in array shape (Matlab flattens matrices with one column and vectors with one
element) are ignored; numbers must agree within tolerance. Exit code 1 if anything differs.
"""
import json
import os
import sys

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "python-sbtab-to-json"))
import sbtab_to_json  # noqa: E402

TOLERANCE = 1e-8


def arr(x):
    return np.array(x, dtype=float).flatten()


def is_numeric(x):
    a = np.array(x, dtype=object).flatten()
    return all(v is None or isinstance(v, (int, float, bool)) for v in a)


def walk(a, b, path, out):
    if isinstance(a, dict) and isinstance(b, dict):
        for key in sorted(set(a) | set(b)):
            if key not in a or key not in b:
                out.append(f"{path}/{key}: only in {'Matlab' if key in a else 'Python'}")
            else:
                walk(a[key], b[key], f"{path}/{key}", out)
        return
    if is_numeric(a) and is_numeric(b):
        x, y = arr(a), arr(b)
        if x.size != y.size:
            out.append(f"{path}: size {x.size} vs {y.size}")
        elif (np.isnan(x) != np.isnan(y)).any():
            out.append(f"{path}: {int((np.isnan(x) != np.isnan(y)).sum())} entries are undefined in only one file")
        else:
            ok = ~np.isnan(x)
            if ok.any():
                d = np.abs(x[ok] - y[ok]) / (1 + np.abs(x[ok]))
                if d.max() > TOLERANCE:
                    out.append(f"{path}: max relative difference {d.max():.3g}")
        return
    if a != b and [a] != b and a != [b]:
        out.append(f"{path}: {str(a)[:50]} vs {str(b)[:50]}")


def main(result_dir):
    manifest = json.load(open(os.path.join(result_dir, "manifest.json")))
    failed = 0
    for case in manifest:
        as_list = lambda x: [x] if isinstance(x, str) else list(x)
        p = sbtab_to_json.build_problem(
            case["network"], case["kinetic"], case["state"], case["constraints"] or None,
            as_list(case["mean"]), as_list(case["std"]), use_kinetic_data=case["use_kinetic_data"],
        )
        sbtab_to_json.write_json(p, os.path.join(result_dir, case["name"] + "_python.json"))
        m = json.load(open(os.path.join(result_dir, case["name"] + "_problem.json")))
        out = []
        walk(m, json.loads(json.dumps(p)), "", out)
        print(f"{case['name']:32s} {'identical' if not out else str(len(out)) + ' differences'}")
        for line in out[:12]:
            print("    " + line)
        failed += bool(out)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
