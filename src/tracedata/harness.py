"""HARNESS / DIAGNOSTIC (factor-aware): local data via PARI, instances, curve family.  Never imported by
src/tracedata/reductions.py or src/tracedata/decode.py."""
from __future__ import annotations

import random
from functools import lru_cache
from typing import Dict, List, Tuple

import cypari2
import gmpy2

from src.common.instances import balanced_semiprime, random_prime

pari = cypari2.Pari()
pari.allocatemem(64 * 10 ** 6, 10 ** 9, silent=True)
pari.default("nbthreads", 1)

# Weierstrass a-invariants [a1, a2, a3, a4, a6]
CURVES: Dict[str, List[int]] = {
    "C1_11a1": [0, -1, 1, -10, -20],      # torsion Z/5
    "C2_37a1": [0, 0, 1, -1, 0],
    "C3_2_3": [0, 0, 0, 2, 3],            # torsion Z/2
    "CM1_32a2": [0, 0, 0, -1, 0],         # y^2 = x^3 - x, CM by Z[i]
    "T7_26b1": [1, -1, 1, -3, 3],         # torsion Z/7
    "C5_92": [0, 0, 0, -1, 1],            # y^2 = x^3 - x + 1, conductor 92
}
TORSION_M = {"C1_11a1": 5, "T7_26b1": 7, "C3_2_3": 2}


@lru_cache(maxsize=None)
def ell(name: str):
    return pari.ellinit(CURVES[name])


def disc(name: str) -> int:
    return int(ell(name)[11])


def good(N: int, name: str) -> bool:
    from math import gcd
    return gcd(N, 6 * disc(name)) == 1


def ap(name: str, p: int) -> int:
    return int(pari.ellap(ell(name), p))


def a_pk(a: int, p: int, k: int) -> int:
    """Hecke recursion a_{p^k} = a_p a_{p^{k-1}} - p a_{p^{k-2}} (good p, weight 2)."""
    u0, u1 = 1, a
    if k == 0:
        return 1
    for _ in range(k - 1):
        u0, u1 = u1, a * u1 - p * u0
    return u1


def cm_weight3_ap(p: int) -> int:
    """a_p of eta(4z)^6 in S_3(Gamma_0(16), chi_-4): 0 if p = 3 mod 4, else 2(a^2 - b^2) with p = a^2 + b^2, a odd."""
    if p % 4 == 3:
        return 0
    for a in range(1, int(gmpy2.isqrt(p)) + 1, 2):
        b2 = p - a * a
        r = int(gmpy2.isqrt(b2))
        if r * r == b2:
            return 2 * (a * a - b2)
    raise ValueError("no representation")


def cm_weight3_ap_fast(p: int) -> int:
    """Same via Cornacchia (PARI qfbsolve), for large p."""
    if p % 4 == 3:
        return 0
    sol = pari.qfbsolve(pari.Qfb(1, 0, 1), p)
    x, y = int(sol[0]), int(sol[1])
    a, b = (x, y) if x % 2 else (y, x)
    return 2 * (a * a - b * b)


def eta4_6_coeffs(n: int) -> List[int]:
    """q-expansion coefficients of eta(4z)^6 = q prod (1 - q^{4m})^6 up to q^n (for tests)."""
    c = [0] * (n + 1)
    c[1] = 1
    for m in range(1, n // 4 + 1):
        for _ in range(6):
            for k in range(n, 4 * m - 1, -1):
                c[k] -= c[k - 4 * m]
    return c


def roots_psi3(name: str, p: int) -> int:
    E = ell(name)
    f = pari.elldivpol(E, 3)
    return int(len(pari.polrootsmod(f, p)))


# ------------------------------------------------------------------ the fixed curve family for E3

@lru_cache(maxsize=1)
def family(maxc: int = 200) -> List[Tuple[int, int]]:
    L = [(A, B) for A in range(1, maxc + 1) for B in range(1, maxc + 1) if 4 * A ** 3 + 27 * B ** 2 != 0]
    L.sort(key=lambda ab: (max(ab), ab[0], ab[1]))
    return L


_FAM: Dict[Tuple[int, int], object] = {}


def fam_ell(A: int, B: int):
    E = _FAM.get((A, B))
    if E is None:
        E = pari.ellinit([0, 0, 0, A, B])
        _FAM[(A, B)] = E
    return E


def fam_traces(p: int, K: int) -> List[int]:
    # good reduction: |4A^3 + 27B^2| < 2^25 for A, B <= 200, so every odd prime p > 2^25 is good
    return [int(pari.ellap(fam_ell(A, B), p)) for (A, B) in family()[:K]]


def max_family_disc(K: int) -> int:
    return max(abs(16 * (4 * A ** 3 + 27 * B ** 2)) for (A, B) in family()[:K])


def census_primes(bb: int = 16) -> List[int]:
    lo, hi = 3 << (bb - 2), 1 << bb
    out = []
    x = int(gmpy2.next_prime(lo - 1))
    while x < hi:
        out.append(x)
        x = int(gmpy2.next_prime(x))
    return out


def prime_list(rng: random.Random, bits: int, count: int) -> List[int]:
    out = []
    while len(out) < count:
        p = random_prime(rng, bits)
        if p not in out:
            out.append(p)
    return out
