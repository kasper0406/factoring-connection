"""POST HOC (not pre-registered) diagnostics for E2's two failed predictions.

(i) P2b: is the Gauss-split deficit explained by coprimality of (u, v)? Null = independent same-bit-length integers
    redrawn until coprime.
(ii) P2d: are multiplier splits duplicates? The split (u', v') of t*s is a duplicate relation when
    u'/v' = s*u/v as rationals, i.e. u'*v == s*u*v' (then r^2 s v' = u' is the same relation as r^2 v = u).
Run: .venv/bin/python -m src.smoothbuild.posthoc_e2   -> experiments/raw/smoothbuild_E2_posthoc.jsonl
"""
from __future__ import annotations

import json
import math
import random
from fractions import Fraction

from src.common.instances import balanced_semiprime

from .core import gauss_split, primes_upto
from .e2_class import is_smooth, same_bitlength

OUT = "experiments/raw/smoothbuild_E2_posthoc.jsonl"


def main() -> None:
    mults = primes_upto(60)[:16]
    with open(OUT, "w") as fh:
        for bits in (40, 48, 56, 64):
            y = 1 << (bits // 4)
            agg = {"hits": 0, "null_indep": 0, "null_coprime": 0, "trials": 0,
                   "m2_targets": 0, "m2_hit_targets": 0, "m2_split_hits": 0, "m2_distinct_hits": 0,
                   "m2_dup_of_base": 0, "m2_splits": 0}
            for seed in range(6):
                N = balanced_semiprime(bits, seed).N
                rng = random.Random(f"smoothbuild-E2post|{bits}|{seed}")
                for _ in range(20000):
                    r = rng.randrange(2, N - 1)
                    t = r * r % N
                    u, v = gauss_split(t, N)
                    agg["trials"] += 1
                    agg["hits"] += is_smooth(u, y) and is_smooth(v, y)
                    a, b = same_bitlength(rng, u), same_bitlength(rng, v)
                    agg["null_indep"] += is_smooth(a, y) and is_smooth(b, y)
                    while True:
                        a, b = same_bitlength(rng, u), same_bitlength(rng, v)
                        if math.gcd(a, b) == 1:
                            break
                    agg["null_coprime"] += is_smooth(a, y) and is_smooth(b, y)
                for _ in range(2000):
                    r = rng.randrange(2, N - 1)
                    t = r * r % N
                    u0, v0 = gauss_split(t, N)
                    base = Fraction(u0, v0)
                    rels = set()
                    agg["m2_targets"] += 1
                    for s in mults:
                        u, v = gauss_split(t * s % N, N)
                        agg["m2_splits"] += 1
                        if v and Fraction(u, v) == s * base:
                            agg["m2_dup_of_base"] += 1
                        if v and is_smooth(u, y) and is_smooth(v, y):
                            agg["m2_split_hits"] += 1
                            rels.add(Fraction(u, v) / s)  # relation r^2 = (u/v)/s
                    agg["m2_hit_targets"] += bool(rels)
                    agg["m2_distinct_hits"] += len(rels)
            agg.update({"bits": bits, "y": y,
                        "O/E_indep": agg["hits"] / agg["null_indep"],
                        "O/E_coprime": agg["hits"] / agg["null_coprime"],
                        "dup_of_base_frac": agg["m2_dup_of_base"] / agg["m2_splits"],
                        "distinct_per_hit": agg["m2_distinct_hits"] / max(1, agg["m2_split_hits"])})
            fh.write(json.dumps(agg) + "\n")
            print(json.dumps(agg), flush=True)


if __name__ == "__main__":
    main()
