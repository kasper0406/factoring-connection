"""Exact index of the integer-valued (guessed-class) lattice L_{r,M} over the
monomial lattice, for several supports and moduli (HNF; MODEL CODE).

  .venv/bin/python -m src.cycloharvey.kappa_table

Writes experiments/analysis/cycloharvey_kappa.json.  N is a fixed public
size parameter (1000003*1000033); only N mod M^e enters.
"""
import json
import math
import os

from .lattice import class_lattice_index

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "experiments", "analysis", "cycloharvey_kappa.json")


def main():
    N = 1000003 * 1000033
    out = []
    for I in [(-1, 1), (1, 2), (-1, 1, 2), (-2, -1, 1, 2), (1,), (-1,), (-3, -1, 1, 3),
              (-1, 1, 2, 3), (-3, -2, -1, 1, 2, 3)]:
        for M in [5, 7, 11, 13, 4, 8, 9, 12, 30]:
            idx = [class_lattice_index(I, N, r, M)["index"] for r in range(1, M) if math.gcd(r, M) == 1]
            ext = None
            if len(I) <= 3 and M <= 8:
                ext = sorted(set(class_lattice_index(I, N, r, M, emax=len(I) + 1)["index"]
                                 for r in range(1, M) if math.gcd(r, M) == 1))
            lg = sorted(set(round(math.log(i) / math.log(M), 3) for i in idx))
            rec = {"I": I, "M": M, "KK1_2": len(I) * (len(I) + 1) // 2, "log_M_index": lg,
                   "emax_K_plus_1_indices": ext, "indices": sorted(set(idx))}
            print(rec, flush=True)
            out.append(rec)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1, default=str)


if __name__ == "__main__":
    main()
