"""HARNESS: Track A (c) false-zero catalogue and section-structure census.

For small semiprimes N = p*q (and a few Carmichael / 3-prime N), for every
M = 2^k <= 2N and every residue r with S_{r,M}(F_N) != 0 (explicit binomials),
test whether the section is detected (nonzero) by cheap quotient evaluations:
  scalar_all   : S(c) != 0 for SOME c in R_N  (exhaustive, decided per prime
                 component: a function F_l -> F_l vanishes iff S(0)=0 and the
                 sums of coefficients over exponent classes j = e mod (l-1),
                 j >= 1, all vanish -- harness criterion, cross-checked by brute
                 force on small N in tests)
  aks(s,a)     : S(Z) mod (N, Z^s - a) != 0, for s in S_LIST, a in A_LIST
  monic2       : S mod (N, h) for T random monic quadratics h (Agrawal-Biswas /
                 Arvind-Mukhopadhyay-Srinivasan style random monic quotient)
A "false zero" for a test family = nonzero section that every test in the
family declares zero.

Also records the census used for Props. C1-C3 (notes/trackA_analysis.md):
number of nonzero classes, max number of terms per nonzero section, number of
pure classes (hit by exactly one prime's support).

Usage:
  .venv/bin/python -m src.trackA.zero_tests --maxN 4000 --out experiments/raw/trackA_zero_tests.jsonl
"""
from __future__ import annotations

import argparse
import random
from math import comb, gcd
from typing import Dict, List

from ..common.instances import CARMICHAEL_SMALL, is_prime
from ..common.metrics import append_jsonl
from ..seed_checks import factorize_reference

S_LIST = (2, 3, 4)
A_LIST = (1, 2, 3)


def sections_of(N: int, coeffs: List[int], M: int) -> Dict[int, Dict[int, int]]:
    out: Dict[int, Dict[int, int]] = {}
    for i in range(1, N):
        a = coeffs[i]
        if a:
            out.setdefault(i % M, {})[i // M] = a
    return out


def vanishes_on_field(sec: Dict[int, int], l: int) -> bool:
    """Does sum_j s_j c^j vanish for all c in F_l?"""
    s0 = sec.get(0, 0) % l
    if s0:
        return False
    acc = [0] * (l - 1)
    for j, a in sec.items():
        if j >= 1:
            acc[j % (l - 1)] = (acc[j % (l - 1)] + a) % l
    return not any(acc)


def scalar_all_zero(sec, primes) -> bool:
    return all(vanishes_on_field(sec, l) for l in primes)


def aks_zero(sec, N, s, a) -> bool:
    w = [0] * s
    for j, x in sec.items():
        w[j % s] = (w[j % s] + x * pow(a, j // s, N)) % N
    return not any(w)


def monic2_zero(sec, N, h) -> bool:
    u, v = h  # Z^2 = -u Z - v
    jmax = max(sec)
    z = (1, 0)  # Z^0 = 1 + 0 Z
    res0 = res1 = 0
    for j in range(jmax + 1):
        if j in sec:
            res0 = (res0 + sec[j] * z[0]) % N
            res1 = (res1 + sec[j] * z[1]) % N
        # multiply z by Z: (c0 + c1 Z) Z = c0 Z + c1 Z^2 = -c1 v + (c0 - c1 u) Z
        z = ((-z[1] * v) % N, (z[0] - z[1] * u) % N)
    return res0 == 0 and res1 == 0


def census(N: int, primes, T: int = 4, seed: int = 0) -> Dict:
    rng = random.Random(f"zt|{N}|{seed}")
    coeffs = [comb(N, i) % N for i in range(N + 1)]
    rows = []
    M = 1
    while M <= 2 * N:
        secs = sections_of(N, coeffs, M)
        nz = len(secs)
        maxterms = max((len(s) for s in secs.values()), default=0)
        pure = 0
        fz_scalar = fz_aks = fz_monic = fz_all = 0
        examples = []
        hs = [(rng.randrange(N), rng.randrange(N)) for _ in range(T)]
        for r, sec in secs.items():
            owners = {l for l in primes if any(x % l for x in sec.values())}
            if len(owners) == 1:
                pure += 1
            sz = scalar_all_zero(sec, primes)
            az = all(aks_zero(sec, N, s, a) for s in S_LIST for a in A_LIST)
            mz = all(monic2_zero(sec, N, h) for h in hs)
            fz_scalar += sz
            fz_aks += az
            fz_monic += mz
            fz_all += sz and az and mz
            if sz and len(examples) < 3:
                examples.append({"r": r, "terms": len(sec), "exps": sorted(sec)[:6]})
        rows.append({"M": M, "nonzero_classes": nz, "all_classes_hit": nz == M,
                     "max_terms": maxterms, "pure_classes": pure,
                     "false_zero_scalar": fz_scalar, "false_zero_aks": fz_aks,
                     "false_zero_monic2": fz_monic, "false_zero_all": fz_all,
                     "scalar_examples": examples})
        M *= 2
    return {"N": N, "primes": list(primes), "levels": rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxN", type=int, default=4000)
    ap.add_argument("--out", default="experiments/raw/trackA_zero_tests.jsonl")
    a = ap.parse_args()
    Ns = []
    for N in range(15, a.maxN + 1, 2):
        f = factorize_reference(N)
        if len(f) == 2 and all(e == 1 for e in f.values()):
            p, q = sorted(f)
            if q < 2 * p:  # same-size (balanced) family
                Ns.append((N, "balanced", [p, q]))
    for N in CARMICHAEL_SMALL[:5]:
        Ns.append((N, "carmichael", sorted(factorize_reference(N))))
    for N in (3 * 5 * 7, 7 * 11 * 13, 11 * 13 * 17):
        Ns.append((N, "three_prime", sorted(factorize_reference(N))))
    for N, fam, primes in Ns:
        rec = census(N, primes)
        rec["family"] = fam
        rec["experiment"] = "trackA_zero_tests"
        append_jsonl(a.out, rec)
        tot = {k: sum(r[k] for r in rec["levels"]) for k in
               ("false_zero_scalar", "false_zero_aks", "false_zero_monic2", "false_zero_all")}
        print(N, fam, primes, tot, flush=True)


if __name__ == "__main__":
    main()
