"""Pre-registered search (notes/cycloharvey_prereg.md, sections 2-4).

  .venv/bin/python -m src.cycloharvey.search analytic
  .venv/bin/python -m src.cycloharvey.search constructive [--procs 28]
  .venv/bin/python -m src.cycloharvey.search guessed [--procs 28]
  .venv/bin/python -m src.cycloharvey.search summarize

MODEL CODE: N is a size parameter (N = 2^n + 2^(n/2) + 1); no factoring.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

from src.common.metrics import append_jsonl, read_jsonl
from .covering import balanced_range, class_basis, harvey_cost, monomial_basis, sweep
from .lattice import (all_supports, class_lattice, exponent_with_guesses,
                      support_stats)

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW_C = os.path.join(ROOT, "experiments", "raw", "cycloharvey_constructive.jsonl")
RAW_G = os.path.join(ROOT, "experiments", "raw", "cycloharvey_guessed.jsonl")
ANA = os.path.join(ROOT, "experiments", "analysis", "cycloharvey_search.json")
CFG = os.path.join(ROOT, "experiments", "configs", "cycloharvey_search.json")

CONFIG = {
    "supports": [[1], [-1, 1], [1, 2], [-1, 1, 2], [-2, -1, 1], [-2, -1, 1, 2],
                 [-3, -1, 1, 3], [-1, 1, 2, 3], [-3, -2, -1, 1, 2, 3]],
    "sizes_full": [40, 48, 56, 64],
    "sizes_sampled": [72, 80],
    "sample_positions": [[1, 8], [3, 8], [5, 8], [7, 8]],
    "sample_rel_len_log2": -6,
    "W_log2_range": "ceil(0.12 n) .. ceil(0.32 n)",
    "fit_sizes": [56, 64, 72, 80],
    "guess_M": [5, 12, 60],
    "guess_sizes": [48, 56],
    "guess_class": "smallest unit r = 1",
    "balanced_ratio": 2,
}


def N_of(n: int) -> int:
    return (1 << n) + (1 << (n // 2)) + 1


def W_grid(n: int):
    return [1 << j for j in range(math.ceil(0.12 * n), math.ceil(0.32 * n) + 1)]


def _task(args):
    I, n, W, M, r = args
    N = N_of(n)
    plo, phi = balanced_range(N)
    if M == 1:
        basis = monomial_basis(I, N)
    else:
        L = class_lattice(I, N, r, M)
        basis = class_basis(I, N, L["rows"], L["scale"])
    t = time.time()
    if n in CONFIG["sizes_sampled"] and M == 1:
        span = phi - plo
        seg = span >> (-CONFIG["sample_rel_len_log2"])
        tot, cov, checks = 0, 0, 0
        for a, b in CONFIG["sample_positions"]:
            x0 = plo + span * a // b - seg // 2
            res = sweep(basis, x0, x0 + seg, W)
            tot += res["n_targets"]
            cov += seg + 1
            checks += res["checks"]
        n_t = tot * (span + 1) / cov          # extrapolated (float: fit input only)
        mode = "sampled"
    else:
        res = sweep(basis, plo, phi, W, M=M, r=r)
        n_t, checks, mode = res["n_targets"], res["checks"], "full"
    return {"I": list(I), "n": n, "W": W, "M": M, "r": r, "n_targets": n_t,
            "checks": checks, "mode": mode, "wall_s": round(time.time() - t, 2)}


def run_constructive(procs: int):
    tasks = []
    for n in CONFIG["sizes_full"] + CONFIG["sizes_sampled"]:
        for I in CONFIG["supports"]:
            for W in W_grid(n):
                tasks.append((tuple(I), n, W, 1, 0))
    _run(tasks, RAW_C, procs)


def run_guessed(procs: int):
    tasks = []
    for n in CONFIG["guess_sizes"]:
        for M in CONFIG["guess_M"]:
            for I in CONFIG["supports"]:
                for W in W_grid(n):
                    tasks.append((tuple(I), n, W, M, 1))
    _run(tasks, RAW_G, procs)


def _run(tasks, path, procs):
    done = set()
    if os.path.exists(path):
        for r in read_jsonl(path):
            done.add((tuple(r["I"]), r["n"], r["W"], r["M"], r["r"]))
    tasks = [t for t in tasks if t not in done]
    # big tasks first
    tasks.sort(key=lambda t: (-t[1], t[2]))
    print(f"{len(tasks)} tasks", flush=True)
    with Pool(procs) as pool:
        for rec in pool.imap_unordered(_task, tasks):
            append_jsonl(path, rec)
            print(rec["I"], rec["n"], rec["W"], rec["M"], rec["n_targets"], rec["wall_s"], flush=True)


def analytic():
    rows = []
    for I in all_supports():
        st = support_stats(I)
        rows.append({"I": list(I), "K": st["K"], "S": st["S"],
                     "exponent": str(st["exponent"]), "exponent_f": float(st["exponent"]),
                     "delta_eff": str(st["delta_eff"]), "a": str(st["a"]),
                     "delta": str(st["delta"]),
                     # with guesses at kappa = K(K+1)/2: slope in mu is exactly 0
                     "mu_slope_F_p": str(exponent_with_guesses(I, Fraction(1)) - st["exponent"]),
                     "mu_slope_cyclotomic": str(exponent_with_guesses(I, Fraction(1), D_exponent=Fraction(1)) - st["exponent"])})
    rows.sort(key=lambda r: r["exponent_f"])
    return rows


def fit(xs, ys):
    nx = len(xs)
    mx, my = sum(xs) / nx, sum(ys) / nx
    sxx = sum((x - mx) ** 2 for x in xs)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx


def summarize():
    out = {"analytic_top10": analytic()[:10]}
    recs = read_jsonl(RAW_C) if os.path.exists(RAW_C) else []
    by = {}
    for r in recs:
        by.setdefault((tuple(r["I"]), r["n"]), []).append(r)
    cons = {}
    for (I, n), rs in by.items():
        best = None
        for r in rs:
            nt = max(1, math.ceil(r["n_targets"]))
            c = harvey_cost(nt, r["W"])
            if best is None or c["cost"] < best["cost"]:
                best = {"cost": c["cost"], "m": c["m"], "W": r["W"], "n_targets": nt}
        cons.setdefault(str(list(I)), {})[n] = best
    slopes = {}
    for I, d in cons.items():
        ns = [n for n in CONFIG["fit_sizes"] if n in d]
        if len(ns) >= 3:
            slopes[I] = fit(ns, [math.log2(d[n]["cost"]) for n in ns])
        allns = sorted(d)
        if len(allns) >= 3:
            slopes[I + "_all"] = fit(allns, [math.log2(d[n]["cost"]) for n in allns])
    out["constructive"] = {I: {str(n): v for n, v in sorted(d.items())} for I, d in cons.items()}
    out["constructive_slopes"] = slopes
    g = read_jsonl(RAW_G) if os.path.exists(RAW_G) else []
    gb = {}
    from .calibration import phi
    for r in g:
        key = (str(r["I"]), r["n"], r["M"])
        c = harvey_cost(max(1, r["n_targets"]), r["W"])
        tot = phi(r["M"]) * c["cost"]                   # phi(M) independent runs
        if key not in gb or tot < gb[key]["cost_total"]:
            gb[key] = {"cost_total": tot, "per_class": c["cost"], "W": r["W"],
                       "n_targets": r["n_targets"]}
    out["guessed"] = {f"{k[0]}|n={k[1]}|M={k[2]}": v for k, v in sorted(gb.items())}
    with open(ANA, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"slopes": slopes}, indent=1))


if __name__ == "__main__":
    cmd = sys.argv[1]
    procs = int(sys.argv[sys.argv.index("--procs") + 1]) if "--procs" in sys.argv else 28
    os.makedirs(os.path.dirname(CFG), exist_ok=True)
    with open(CFG, "w") as f:
        json.dump(CONFIG, f, indent=1)
    if cmd == "analytic":
        rows = analytic()
        print(json.dumps(rows[:12], indent=1))
    elif cmd == "constructive":
        run_constructive(procs)
    elif cmd == "guessed":
        run_guessed(procs)
    elif cmd == "summarize":
        summarize()
