"""Experiment E1 (pre-registered, notes/rank3_prereg.md): exact checks of T1 (divisibility of Delta_F),
T3 (box count) and T4 (rank-r Liouville) on linked rank-r modules. Harness code (factors only toy resultants
of the SNFS-like family to choose N; never factors N in an algorithm).

Run: .venv/bin/python -m src.rank3.e1     (writes experiments/raw/rank3_E1.jsonl)
"""
from __future__ import annotations

import json
import math
import random
import statistics
import sys
import time
from typing import Dict, List

import flint
import sympy

from src.rank3.forms import (box_abs_values, count_bound_terms, delta_F, discriminant, liouville_check,
                             predicted_valuation, resultant, vp)

SEED = 20261020
OUT = "experiments/raw/rank3_E1.jsonl"
CELLS = [(3, 3, 3), (4, 3, 3), (4, 4, 3), (4, 4, 4)]
CELLS_C = [(3, 3, 3), (4, 4, 3)]


def irreducible(c) -> bool:
    fac = flint.fmpz_poly(c).factor()[1]
    return len(fac) == 1 and fac[0][1] == 1


def rand_prime(rnd: random.Random, bits: int) -> int:
    return int(sympy.nextprime(rnd.randrange(1 << (bits - 1), 1 << bits)))  # primality test used (harness)


def poly_with_root(rnd: random.Random, m: int, d: int, N: int) -> List[int]:
    c = [0] * d + [1]
    for k in range(1, d):
        c[k] = rnd.randrange(-5, 6)
    val = sum(c[k] * pow(m, k, N) for k in range(1, d + 1)) % N
    c0 = (-val) % N
    c[0] = c0 - N if c0 > N // 2 else c0
    return c


def inst_A(rnd: random.Random, d: int, d2: int):
    while True:
        p, q = rand_prime(rnd, 10), rand_prime(rnd, 11)
        if p == q:
            continue
        N = p * q
        m = rnd.randrange(N)
        f, g = poly_with_root(rnd, m, d, N), poly_with_root(rnd, m, d2, N)
        if f == g or not irreducible(f) or not irreducible(g):
            continue
        dd = discriminant(f) * discriminant(g)
        if dd % p == 0 or dd % q == 0 or resultant(f, g) == 0:
            continue
        return {"N": N, "primes": [p, q], "m": m, "f": f, "g": g}


def inst_B(rnd: random.Random, d: int, d2: int):
    D = d + d2
    while True:
        f = [rnd.randrange(-3, 4) for _ in range(d)] + [1]
        g = [rnd.randrange(-3, 4) for _ in range(d2)] + [1]
        if f == g or not irreducible(f) or not irreducible(g):
            continue
        R = abs(resultant(f, g))
        if R == 0:
            continue
        dd = discriminant(f) * discriminant(g)
        primes = [p for p in sympy.factorint(R) if p > D and dd % p != 0]  # harness: toy resultant
        N = math.prod(primes) if primes else 1
        if N < 1000:
            continue
        return {"N": N, "primes": sorted(primes), "m": None, "f": f, "g": g}


def module_power(r: int):
    P = [[0] * j + [1] for j in range(r)]
    return P, P


def module_C(rnd: random.Random, r: int, d: int, d2: int, m: int, N: int):
    mm = m - N if m > N // 2 else m
    P, Q = [], []
    for _ in range(r):
        Pj = [rnd.randrange(-3, 4) for _ in range(min(d, d2))]
        s = rnd.randrange(-2, 3)
        Qj = Pj[:] + [0] * max(0, 2 - len(Pj))
        Qj[0] -= s * mm
        Qj[1] += s
        P.append(Pj)
        Q.append(Qj)
    return P, Q


def run_one(kind: str, rnd: random.Random, d: int, d2: int, r: int, seed: int) -> Dict:
    t0 = time.time()
    inst = inst_B(rnd, d, d2) if kind == "B" else inst_A(rnd, d, d2)
    f, g, N = inst["f"], inst["g"], inst["N"]
    P, Q = module_C(rnd, r, d, d2, inst["m"], N) if kind == "C" else module_power(r)
    # linking check (A, C): P_j(m) == Q_j(m) mod N
    if inst["m"] is not None:
        m = inst["m"]
        ev = lambda h: sum(c * pow(m, k, N) for k, c in enumerate(h)) % N
        assert all(ev(Pj) == ev(Qj) for Pj, Qj in zip(P, Q))
    D = d + d2
    info = delta_F(f, P, g, Q, r)
    pred = predicted_valuation(r, D)
    vals_p = {str(p): vp(info["Q2"], p) * info["exp"] for p in inst["primes"]}
    t1_ok = info["Q2"] != 0 and all(v >= pred for v in vals_p.values())
    liou = liouville_check(f, P, g, Q, r, info)
    Es = [3] if r == 4 else [3, 6]
    counts = []
    for E in Es:
        vals = box_abs_values(f, P, g, Q, r, E)
        lo = max(0, vals[0].bit_length() - 1)
        hi = vals[-1].bit_length()
        import bisect
        for k in range(lo, hi + 1):
            X = 1 << k
            c = bisect.bisect_right(vals, X)
            main, err = count_bound_terms(r, D, E, X, info)
            ok = bool(c <= main.upper() + err)
            eligible = bool(main.lower() > err)
            sat = float(c / main.mid()) if eligible else None
            counts.append({"E": E, "log2X": k, "count": c, "main_mid": float(main.mid()), "err": err,
                           "ok": ok, "eligible": eligible, "sat": sat})
    sats = [x["sat"] for x in counts if x["sat"] is not None]
    return {"kind": kind, "d": d, "d2": d2, "r": r, "seed": seed, "N": N, "primes": inst["primes"], "f": f,
            "g": g, "P": P, "Q": Q, "Delta_Q2_bits": abs(info["Q2"]).bit_length(), "Delta_exp": info["exp"],
            "vp": vals_p, "pred": pred, "T1_ok": t1_ok, "T4": liou,
            "T3_ok": all(x["ok"] for x in counts), "max_sat": max(sats) if sats else None,
            "n_eligible": len(sats), "counts": counts, "secs": round(time.time() - t0, 2)}


def main() -> None:
    rows = []
    with open(OUT, "w") as fh:
        for kind, cells in (("A", CELLS), ("B", CELLS), ("C", CELLS_C)):
            for (d, d2, r) in cells:
                for s in range(3):
                    rnd = random.Random(f"{SEED}-{kind}-{d}-{d2}-{r}-{s}")
                    row = run_one(kind, rnd, d, d2, r, s)
                    rows.append(row)
                    fh.write(json.dumps(row) + "\n")
                    fh.flush()
                    print(kind, d, d2, r, s, "N", row["N"], "vp", row["vp"], ">=", row["pred"], "T1", row["T1_ok"],
                          "T3", row["T3_ok"], "T4", row["T4"]["holds"], "sat", row["max_sat"], row["secs"], "s",
                          flush=True)
    # summary of predictions
    p1 = sum(not r_["T1_ok"] for r_ in rows)
    p2 = sum(not r_["T3_ok"] for r_ in rows)
    p3 = sum(not r_["T4"]["holds"] for r_ in rows)
    B = [r_ for r_ in rows if r_["kind"] == "B" and r_["max_sat"] is not None]
    p4 = sum(r_["max_sat"] >= 1e-3 for r_ in B)
    print(f"P1 violations {p1}; P2 violations {p2}; P3 violations {p3}; "
          f"P4: {p4}/{len(B)} eligible B instances with sat >= 1e-3")
    for (d, d2, r) in CELLS:
        sa = [x["max_sat"] for x in rows if x["kind"] == "A" and (x["d"], x["d2"], x["r"]) == (d, d2, r)
              and x["max_sat"] is not None]
        sb = [x["max_sat"] for x in rows if x["kind"] == "B" and (x["d"], x["d2"], x["r"]) == (d, d2, r)
              and x["max_sat"] is not None]
        print("P5 cell", (d, d2, r), "median A", statistics.median(sa) if sa else None,
              "median B", statistics.median(sb) if sb else None)


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
