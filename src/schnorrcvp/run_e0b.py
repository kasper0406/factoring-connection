"""E0b (MODEL; added after the pre-registration and disclosed as such): concrete cost of the prime-lattice method.

For each (bits, n): log2 V(n) = number of lattice vectors an IDEAL enumerator (increasing |w|, unit cost) must
output to see n + 1 fac-relations (psi.ideal_enumeration_cost, H_S naive -- an over-estimate of primitive
relations, so V is optimistic for the method).  Lattice cost model (heuristic, optimistic for the method):
one sieve-CVP call in dimension n costs 2^{0.292 n} and yields up to 2^{0.2075 n} distinct close vectors
(BDGL sieving list size), so
    log2 T(n) = 0.292 n + max(0, log2 V(n) - 0.2075 n).
Reference points (o(1) dropped, crude): QS L[1/2, 1], GNFS L[1/3, (64/9)^{1/3}].

Usage: .venv/bin/python -m src.schnorrcvp.run_e0b -> experiments/raw/schnorrcvp_E0b.jsonl
"""
from __future__ import annotations

import json
import math

import mpmath as mp

from .core import first_primes
from .psi import ideal_enumeration_cost

BITS = [64, 128, 200, 256, 400, 512, 800, 1024]
NS = [10, 20, 30, 50, 75, 100, 150, 200, 300, 400, 600, 800, 1200]


def Lcrude(bits, a, c):
    x = bits * math.log(2)
    return c * x ** a * math.log(x) ** (1 - a) / math.log(2)


def main():
    out = open("experiments/raw/schnorrcvp_E0b.jsonl", "w")
    for bits in BITS:
        best = None
        for n in NS:
            r = ideal_enumeration_cost(first_primes(n), bits * mp.log(2))
            lv = r.get("log2V")
            T = None if lv is None else 0.292 * n + max(0.0, lv - 0.2075 * n)
            rec = {"bits": bits, "n": n, "log2V": lv, "log2T_sieve": T, "lnX": r.get("lnX"),
                   "log2_QS": Lcrude(bits, 0.5, 1.0), "log2_GNFS": Lcrude(bits, 1 / 3, (64 / 9) ** (1 / 3))}
            out.write(json.dumps(rec) + "\n"); out.flush()
            print(rec, flush=True)
            if T is not None and (best is None or T < best[1]):
                best = (n, T)
            if lv is not None and 0.292 * n > (best[1] if best else 1e9) + 30:
                break
        out.write(json.dumps({"bits": bits, "best_n": best[0] if best else None,
                              "best_log2T": best[1] if best else None}) + "\n"); out.flush()


if __name__ == "__main__":
    main()
