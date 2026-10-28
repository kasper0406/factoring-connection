"""HARNESS / DIAGNOSTIC (factor-aware): instances, traces via PARI, class numbers.  Never imported by
src/cuspcoef/decode.py or src/cuspcoef/family.py."""
from __future__ import annotations

import random
from typing import List, Tuple

import cypari2
import gmpy2

from src.common.instances import random_prime

pari = cypari2.Pari()
pari.allocatemem(256 * 10 ** 6, 10 ** 9, silent=True)
pari.default("nbthreads", 1)

_ELL = {}


def ell_AB(A: int, B: int):
    E = _ELL.get((A, B))
    if E is None:
        E = pari.ellinit([0, 0, 0, A, B])
        _ELL[(A, B)] = E
    return E


def ap_AB(A: int, B: int, p: int) -> int:
    return int(pari.ellap(ell_AB(A, B), p))


def ap_curve(ainvs: List[int], p: int) -> int:
    key = tuple(ainvs)
    E = _ELL.get(key)
    if E is None:
        E = pari.ellinit(list(ainvs))
        _ELL[key] = E
    return int(pari.ellap(E, p))


def instance(seed: int, nbits: int) -> Tuple[int, int, int]:
    """p, q independent uniform primes in [2^{n/2-1}, 2^{n/2}) (random_prime of n/2 bits), p != q."""
    rng = random.Random(seed)
    while True:
        p = random_prime(rng, nbits // 2)
        q = random_prime(rng, nbits // 2)
        if p != q:
            p, q = min(p, q), max(p, q)
            return p * q, p, q


def hurwitz(D: int) -> int:
    """12 H(D) as an integer (PARI qfbhclassno returns H(D) as a rational with denominator dividing 6... use 12H)."""
    return int(12 * pari.qfbhclassno(D))


def is_prime(m: int) -> bool:
    return bool(gmpy2.is_prime(m))
