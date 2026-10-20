"""E1 (pre-registered in notes/resbarrier_prereg.md): exact check of Lemma C2 on GNFS and SNFS pairs.

Harness/diagnostic. Run:  .venv/bin/python -m src.resbarrier.e1
Writes experiments/raw/resbarrier_E1.jsonl and prints a summary. Floats appear only in the summary.
"""
from __future__ import annotations

import json
import os
import random
import statistics
from typing import Dict, List

from src.common.instances import random_prime
from src.resbarrier.forms import (binary_disc, box_values, base_m_pair, c2_bound_holds, c2_saturation,
                                  count_le, polymul, res_check, snfs_pair)

SEED = 20261020
E = 256
OUT = "experiments/raw/resbarrier_E1.jsonl"


def iroot(n: int, k: int) -> int:
    lo, hi = 0, 1 << (n.bit_length() // k + 1)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** k <= n:
            lo = mid
        else:
            hi = mid - 1
    return lo


def gnfs_instance(N: int, d: int):
    m = iroot(N, d + 1)
    for shift in range(0, 64):
        f, g = base_m_pair(N, d, m + shift)
        if f[-1] != 0 and binary_disc(f) != 0:
            return f, g, m + shift
    raise RuntimeError("no base-m pair")


def run_instance(kind: str, N: int, f, g, meta: Dict) -> Dict:
    D = (len(f) - 1) + (len(g) - 1)
    F = polymul(f, g)
    disc = abs(binary_disc(F))
    rc = res_check(f, g, N)
    vals = box_values(F, E)
    lo, hi = vals[0].bit_length() - 1, vals[-1].bit_length()
    rows = []
    viol_N2 = viol_disc = 0
    smax_N2 = smax_disc = 0.0
    for k in range(max(lo, 0), hi + 1, 4):
        X = 1 << k
        c = count_le(vals, X)
        okN = c2_bound_holds(c, D, E, X, N * N)
        okd = c2_bound_holds(c, D, E, X, disc)
        viol_N2 += (not okN)
        viol_disc += (not okd)
        sN = c2_saturation(c, D, E, X, N * N)
        sd = c2_saturation(c, D, E, X, disc)
        smax_N2 = max(smax_N2, sN)
        smax_disc = max(smax_disc, sd)
        rows.append({"log2X": k, "count": c, "ok_N2": int(okN), "ok_disc": int(okd), "sat_N2": sN, "sat_disc": sd})
    return {"kind": kind, **meta, "N": str(N), "f": f, "g": g, "D": D, "E": E, "res_divisible": rc["divisible"],
            "hadamard_ok": rc["hadamard_ok"], "log2_disc": disc.bit_length(), "log2_N2": (N * N).bit_length(),
            "violations_N2": viol_N2, "violations_disc": viol_disc, "smax_N2": smax_N2, "smax_disc": smax_disc,
            "grid": rows}


def main() -> None:
    rng = random.Random(SEED)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    recs: List[Dict] = []
    for n in (40, 56, 72):
        for d in (2, 3, 4):
            for rep in range(3):
                p = random_prime(rng, n // 2)
                q = random_prime(rng, n - n // 2)
                N = p * q
                f, g, m = gnfs_instance(N, d)
                recs.append(run_instance("GNFS", N, f, g, {"n": n, "d": d, "rep": rep, "m": m}))
    for e in (40, 56, 72):
        for d in (2, 3, 4):
            N, f, g = snfs_pair(e, d)
            recs.append(run_instance("SNFS", N, f, g, {"n": e, "d": d, "rep": 0}))
    with open(OUT, "w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
    tot = sum(r["violations_N2"] for r in recs)
    print(f"P1: total violations of C2 (N^2 floor) = {tot}; (exact disc floor) = {sum(r['violations_disc'] for r in recs)}")
    print(f"res divisible all: {all(r['res_divisible'] for r in recs)}, hadamard ok all: {all(r['hadamard_ok'] for r in recs)}")
    for n in (40, 56, 72):
        for d in (2, 3, 4):
            gs = [r["smax_N2"] for r in recs if r["kind"] == "GNFS" and r["n"] == n and r["d"] == d]
            ss = [r["smax_N2"] for r in recs if r["kind"] == "SNFS" and r["n"] == n and r["d"] == d]
            gd = [r["smax_disc"] for r in recs if r["kind"] == "GNFS" and r["n"] == n and r["d"] == d]
            print(f"n={n} d={d}: smax GNFS median {statistics.median(gs):.3g}  SNFS {ss[0]:.3g}  "
                  f"(GNFS with exact disc floor {statistics.median(gd):.3g})")
    g72 = statistics.median([r["smax_N2"] for r in recs if r["kind"] == "GNFS" and r["n"] == 72 and r["d"] == 3])
    s72 = [r["smax_N2"] for r in recs if r["kind"] == "SNFS" and r["n"] == 72 and r["d"] == 3][0]
    print(f"P2: SNFS/GNFS saturation ratio at n=72,d=3: {s72 / g72 if g72 else float('inf'):.3g} (pass if >= 8)")


if __name__ == "__main__":
    main()


def supplementary_irreducible_snfs() -> None:
    """POST-HOC, disclosed (critic A #9): 6/9 pre-registered SNFS instances had reducible f = x^d - 1. Re-run the
    SNFS cells with N = 2^e + 3, f = x^d + 3*2^(dk-e), g = x - 2^k (f(2^k) = 2^(dk-e) N), keeping only
    irreducible f (checked by flint). Writes experiments/raw/resbarrier_E1_supp.jsonl."""
    import flint
    recs = []
    for e in (40, 56, 72):
        for d in (2, 3, 4):
            k = -(-e // d)
            N = (1 << e) + 3
            f = [3 << (d * k - e)] + [0] * (d - 1) + [1]
            g = [-(1 << k), 1]
            fac = flint.fmpz_poly(f).factor()[1]
            irreducible = len(fac) == 1 and fac[0][1] == 1
            if not irreducible:
                continue
            recs.append(run_instance("SNFS_irr", N, f, g, {"n": e, "d": d, "rep": 0}))
    with open("experiments/raw/resbarrier_E1_supp.jsonl", "w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
    for r in recs:
        print(f"SNFS_irr n={r['n']} d={r['d']}: violations={r['violations_N2']} smax_N2={r['smax_N2']:.3g}")
