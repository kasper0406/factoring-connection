"""Experiment E2 (pre-registered, notes/rankdegen_prereg.md): degenerate rank-3 designs vs rank 2, toy sizes.

Measured only, NOT asymptotic. Reuses rank3.e2 by import (same N, base-m pair, lattice pair, boxes, smoothness test).
Harness: N = pq generated with a primality test; designs receive N and N-derived polynomials only.
Run: PYTHONPATH=. .venv/bin/python -m src.rankdegen.e2   (writes experiments/raw/rankdegen_E2.jsonl)
"""
from __future__ import annotations

import json
import random
import time

from flint import fmpz_mat

from src.rank3 import e2 as R3

OUT = "experiments/raw/rankdegen_E2.jsonl"


def short_quadratic_kernel(N: int, m: int):
    """Shortest (LLL) h = h0 + h1 x + h2 x^2 with h(m) = 0 mod N and h2 != 0 (a quadratic element of ker phi')."""
    L = fmpz_mat([[N, 0, 0], [-(m % N), 1, 0], [-(m * m % N), 0, 1]]).lll()
    rows = [[int(L[i, j]) for j in range(3)] for i in range(3)]
    rows.sort(key=lambda h: sum(c * c for c in h))
    for h in rows:
        if h[2] != 0:
            return h
    raise RuntimeError("no quadratic kernel element")


def admissible_quadratic_kernel(N: int, m: int):
    """CORRECTED design (rankdegen_2 #11): the LLL-selected quadratic h with h(m) = 0 mod N that is irreducible over Q
    (non-square discriminant), so not a multiple of a linear kernel element. Searches small combinations of the
    LLL-reduced basis and returns the shortest admissible one found (an LLL-selected candidate, not a proven minimum)."""
    import itertools as it, math as mt
    L = fmpz_mat([[N, 0, 0], [-(m % N), 1, 0], [-(m * m % N), 0, 1]]).lll()
    rows = [[int(L[i, j]) for j in range(3)] for i in range(3)]
    best = None
    for c in it.product(range(-3, 4), repeat=3):
        h = [sum(ci * row[j] for ci, row in zip(c, rows)) for j in range(3)]
        if h[2] == 0:
            continue
        disc = h[1] * h[1] - 4 * h[0] * h[2]
        if disc >= 0 and mt.isqrt(disc) ** 2 == disc:
            continue
        if best is None or sum(x * x for x in h) < sum(x * x for x in best):
            best = h
    return best


def main(corrected: bool = False) -> None:
    rows = []
    out = OUT.replace(".jsonl", "_corrected.jsonl") if corrected else OUT
    with open(out, "w") as fhd:
        for n in (40, 50, 60):
            for s in range(3):
                rnd = random.Random(f"{R3.SEED}-E2-{n}-{s}")
                N, p, q = R3.make_N(rnd, n)
                f_rat, g_rat, m = R3.base_m(N, 3)
                h1, h2 = R3.lattice_pair(N, m, 3)
                mm = m % N
                nu = admissible_quadratic_kernel(N, mm) if corrected else short_quadratic_kernel(N, mm)
                res = {"n": n, "seed": s, "N": N, "m": m, "h1": h1, "h2": h2, "nu": nu}
                t0 = time.time()
                res["R2-rat"] = R3.evaluate(f_rat, g_rat, [[1], [0, -1]], R3.box2())
                res["R2-nn"] = R3.evaluate(h1, h2, [[1], [0, -1]], R3.box2())
                res["D-C1"] = R3.evaluate(h1, h2, [[1], [0, -1], [0]], R3.box3(), gens2=[[1], [0, -1], nu])
                res["D-FP"] = R3.evaluate(h1, h2, [[1], [0, -1], [0]], R3.box3(), gens2=[[1], [0, -1], [-mm, 1]])
                res["D-RK"] = R3.evaluate(f_rat, g_rat, [[1], [0, -1], [0]], R3.box3(), gens2=[[1], [0, -1], [N]])
                res["secs"] = round(time.time() - t0, 1)
                rows.append(res)
                fhd.write(json.dumps(res) + "\n")
                fhd.flush()
                print(n, s, {k: (v["yield"], v["median_log2F"], v["n_cand"]) for k, v in res.items()
                             if isinstance(v, dict)}, res["secs"], "s", flush=True)
    v6 = v7 = 0
    for res in rows:
        best2 = "R2-rat" if res["R2-rat"]["yield"] >= res["R2-nn"]["yield"] else "R2-nn"
        y2 = res[best2]["yield"]
        m2 = res[best2]["median_log2F"]
        yd = max(res[k]["yield"] for k in ("D-C1", "D-FP", "D-RK"))
        v6 += yd > 2 * y2
        v7 += any(res[k]["median_log2F"] is not None and res[k]["median_log2F"] < m2 - 2
                  for k in ("D-C1", "D-FP", "D-RK"))
        print("n", res["n"], "seed", res["seed"], "best2", best2, y2, "best degenerate", yd,
              "ratio", round(yd / max(y2, 1), 4))
    print(f"P6 violations {v6}/{len(rows)}; P7 violations {v7}/{len(rows)}")


if __name__ == "__main__":
    import sys
    main(corrected="--corrected" in sys.argv)
