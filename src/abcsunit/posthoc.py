"""POST-HOC diagnostics (disclosed in notes/abcsunit_results.md §6; not pre-registered).

PH1: dispersion of the E1 pair count.  Every primitive pair (u', v') generates the Ψ(X/v') pairs (g u', g v')
     with g smooth, so pairs = Σ_prim m(v'), m(v') = Ψ(X/v').  Under a Poisson model for primitive pairs the
     standard deviation of the pair count is sqrt(Σ m²); we report z = (pairs − M1)/sqrt(Σ m²).
PH2: Z-rank of the E2 relation set (primitive spanning-forest relations, first R), against Lemma FR.
"""
from __future__ import annotations

import json
import math
import os
import sys
from math import gcd

import flint
import numpy as np

from src.abcsunit import core
from src.abcsunit.run_experiment import X_of, y_of, RAW
from src.common.instances import balanced_semiprime


def main():
    rows = [json.loads(l) for l in open(os.path.join(RAW, "abcsunit_E123.jsonl"))]
    out = open(os.path.join(RAW, "abcsunit_posthoc.jsonl"), "w")
    for r in rows:
        n, kap, seed = r["n"], r["kappa"], r["seed"]
        Kn, Kd = map(int, r["K"].split("/"))
        inst = balanced_semiprime(n, seed)
        N = inst.N
        P = core.primes_upto(y_of(n, kap))
        X = X_of(N, Kn, Kd)
        vals = core.smooth_numbers(X, P)
        bk = core.congruent_buckets(vals, N)
        msq = 0
        nprim = 0
        forest = []
        for b in bk:
            par = list(range(len(b)))

            def f(t):
                while par[t] != t:
                    t = par[t]
                return t

            for j in range(1, len(b)):
                for i in range(j):
                    if gcd(b[i], b[j]) == 1:
                        nprim += 1
                        m = int(np.searchsorted(vals, X // b[j], side="right"))
                        msq += m * m
                        if f(i) != f(j):
                            par[f(i)] = f(j)
                            forest.append((b[i], b[j]))
        forest.sort(key=lambda t: (t[1], t[0]))
        forest = forest[: 3 * len(P) + 30]
        E = [[a - c for a, c in zip(core.exponent_vector(v, P), core.exponent_vector(u, P))] for u, v in forest]
        zrank = flint.fmpz_mat(E).rank() if E else 0
        rec = {"n": n, "kappa": kap, "K": r["K"], "seed": seed, "pairs": r["pairs"], "M1": r["M1"],
               "primitive": nprim, "sd": math.sqrt(msq), "z": (r["pairs"] - r["M1"]) / math.sqrt(msq) if msq else None,
               "E2_relations": len(forest), "E2_Zrank": zrank, "s": len(P), "rank": r["rank"],
               "z_image_dim": r["z_image_dim"], "fruitful": r["fruitful"]}
        out.write(json.dumps(rec) + "\n")
        print(rec, flush=True)


if __name__ == "__main__":
    main()
