"""HARNESS, POST HOC (not pre-registered): smoothness-observable peaks versus B and n (Prop SM).

For balanced semiprimes of n bits and factor-base bounds B, measure the max non-exceptional relative peak of the
centred B-smooth indicator on the orbit of a, its ratio to the orbit density rho_H, and to the random baseline ln r / r.
Usage: .venv/bin/python -m src.periodfind.posthoc_smooth   -> experiments/raw/periodfind_posthoc_smooth.jsonl
"""
from __future__ import annotations

import json
import random

import gmpy2
import numpy as np
import sympy

from src.common.instances import balanced_semiprime
from src.periodfind import spectra as sp
from src.periodfind.run_experiments import pick_bases


def smooth_indicator(y: np.ndarray, B: int) -> np.ndarray:
    z = y.copy()
    for ell in sympy.primerange(2, B + 1):
        while True:
            m = (z % ell == 0) & (z > 1)
            if not m.any():
                break
            z[m] //= ell
    return (z == 1).astype(np.float64)


def main():
    out = open("experiments/raw/periodfind_posthoc_smooth.jsonl", "w")
    for n in (18, 20, 22, 24):
        for seed in range(4):
            inst = balanced_semiprime(n, seed)
            N = inst.N
            p, q = inst.factors
            rng = random.Random(f"periodfind-posthoc|{n}|{seed}")
            bases = pick_bases(N, p, q, n, rng, count=1, max_r=2 ** 23)
            if not bases:
                continue
            a, info = bases[0]
            r = info["r"]
            orb = sp.orbit(a, N, r)
            ja = gmpy2.jacobi(a, N)
            for B in (10, 20, 50, 200, 1000, 5000):
                ind = smooth_indicator(orb, B)
                rho = float(ind.mean())
                if rho in (0.0, 1.0):
                    continue
                F = ind - rho
                met = sp.metrics(sp.spectrum(F), r, ja)
                rec = {"n": n, "seed": seed, "N": N, "a": a, "r": r, "B": B, "piB": int(sympy.primepi(B)),
                       "rho_H": rho, "M_peak": met["M_peak"], "peak_over_rho": met["M_peak"] / rho,
                       "peak_over_baseline": met["M_peak"] / met["rand_baseline"],
                       "argmax_order": met["argmax_order"]}
                out.write(json.dumps(rec) + "\n")
                out.flush()
                print(n, seed, B, "%.3g %.3g %.1f" % (rho, rec["peak_over_rho"], rec["peak_over_baseline"]),
                      flush=True)


if __name__ == "__main__":
    main()
