"""Tables and fits from the cadobench raw data (V5, V6, V7 and the phase-timing table).

Floating point is used for summaries and least-squares fits only.

Usage: python -m cadobench.analyze [--v6]   (prints markdown tables; writes experiments/raw/cadobench_analysis.json)
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys

import numpy as np

from .relstats import filtering_numbers, incidence_stats
from .validate import RAW

GROUPS = ("polyselect", "sieve", "filter", "linalg", "sqrt", "other")


def load_runs():
    runs = {}
    for path in sorted(glob.glob(os.path.join(RAW, "cadobench_c*_run*.json"))):
        r = json.load(open(path))
        runs.setdefault(r["digits"], []).append(r)
    return dict(sorted(runs.items()))


def L13(N: int) -> float:
    """x^(1/3) (log x)^(2/3) with x = ln N (so that L[1/3, c] = exp(c * L13))."""
    x = math.log(N)
    return x ** (1 / 3) * math.log(x) ** (2 / 3)


def phase_table(runs) -> str:
    lines = ["| digits | total wall (s) | polyselect | sieve | filter | linalg | sqrt | other | CPU (s) | peak RSS (MB) |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for d, rs in runs.items():
        for r in rs:
            g = r["wall_by_group_s"]
            tot = r["wall_total_s"]
            cpu = None
            for k, v in r["cado_reported"].items():
                if "Total cpu/elapsed time for entire" in v:
                    cpu = float(v.split()[-1].split("/")[0])
            cells = [f"{g.get(k, 0):.1f} ({100 * g.get(k, 0) / tot:.0f}%)" for k in GROUPS]
            lines.append(f"| c{d} | {tot:.1f} | " + " | ".join(cells) +
                         f" | {cpu:.0f} | {r['memwatch']['peak_tree_rss_bytes'] >> 20} |")
    return "\n".join(lines)


def v5(runs) -> dict:
    out = {}
    for d, rs in runs.items():
        r = rs[0]
        f = filtering_numbers(r["workdir"])
        g = r["wall_by_group_s"]
        out[d] = {**f, "linalg_share": g.get("linalg", 0) / r["wall_total_s"],
                  "sieve_share": g.get("sieve", 0) / r["wall_total_s"],
                  "filter_share": g.get("filter", 0) / r["wall_total_s"],
                  "sieve_wall": g.get("sieve", 0), "linalg_wall": g.get("linalg", 0),
                  "matrix_over_unique": f.get("matrix_N", 0) / f["unique_relations"] if f.get("unique_relations") else None}
    return out


def v7(runs, dmin: int = 70) -> dict:
    ds = [d for d in runs if d >= dmin]
    xs = np.array([L13(int(runs[d][0]["N"])) for d in ds])
    res = {"digits": ds}
    for key, get in (("total", lambda r: r["wall_total_s"]), ("sieve", lambda r: r["wall_by_group_s"]["sieve"]),
                     ("linalg", lambda r: r["wall_by_group_s"]["linalg"])):
        ys = np.log([np.mean([get(r) for r in runs[d]]) for d in ds])
        c, logA = np.polyfit(xs, ys, 1)
        # local constants between consecutive sizes
        loc = [float((ys[k + 1] - ys[k]) / (xs[k + 1] - xs[k])) for k in range(len(ds) - 1)]
        res[key] = {"c_fit": float(c), "logA": float(logA), "local_c": loc}
    return res


def v6(runs, dmin: int = 70) -> dict:
    out = {}
    for d, rs in runs.items():
        if d < dmin:
            continue
        r = rs[0]
        p = r["params"]
        out[d] = incidence_stats(r["workdir"], int(p["tasks.lim0"]), int(p["tasks.lim1"]))
        print(f"[v6] c{d} done", file=sys.stderr, flush=True)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--v6", action="store_true")
    a = ap.parse_args(argv)
    runs = load_runs()
    res = {"phase_table_md": phase_table(runs), "v5": v5(runs), "v7": v7(runs)}
    if a.v6:
        res["v6"] = v6(runs)
        json.dump(res["v6"], open(os.path.join(RAW, "cadobench_v6_incidence.json"), "w"), indent=1)
    print(res["phase_table_md"])
    print(json.dumps({"v5": res["v5"], "v7": res["v7"]}, indent=1))
    json.dump(res, open(os.path.join(RAW, "cadobench_analysis.json"), "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
