"""The mod-M^2 sign-splitting identity (OPERATIONAL-SAFE number theory).

Notation follows Harvey--Hittmeir 2022 (arXiv 2105.11105), Prop. 3.3: weights
(a, b) act on a*q + b*p, p == sigma (mod M), and the lock is

    b * sigma^2 == a * N   (mod M)                                  (lock)

(the user's note writes the same thing as  a r^2 == b N  with the roles of
(a, b) swapped).  Under the lock,

    tau(sigma) := a*N*sigma^{-1} + b*sigma   (mod M^2)

is well defined for sigma mod M and equals a*y + b*x (mod M^2) for EVERY pair of
integers with x*y == N (mod M^2), x == sigma (mod M)  (Lemma 2 in the notes).

Everything here uses only N, M and public small primes; no factor of N is ever
used.  All arithmetic is exact (Python ints).
"""
from __future__ import annotations

from itertools import product
from math import gcd
from typing import Dict, List, Sequence, Tuple


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------
def small_primes(limit: int) -> List[int]:
    out = []
    for n in range(2, limit + 1):
        if all(n % d for d in range(2, int(n ** 0.5) + 1)):
            out.append(n)
    return out


def factor_squarefree(M: int, primes: Sequence[int] = None) -> List[int]:
    """Prime factors of a squarefree public modulus M (trial division by small
    primes; M is a product of small public primes, never related to N)."""
    out, x, d = [], M, 2
    while d * d <= x:
        if x % d == 0:
            out.append(d)
            x //= d
            if x % d == 0:
                raise ValueError("M not squarefree")
        d += 1
    if x > 1:
        out.append(x)
    return out


def sqrt_mod_prime(c: int, m: int) -> List[int]:
    """All square roots of c mod a small prime m (brute force; m is tiny)."""
    c %= m
    return [x for x in range(m) if (x * x - c) % m == 0]


def crt_idempotents(moduli: Sequence[int]) -> List[int]:
    """E_i mod Q = prod(moduli): E_i == 1 mod moduli[i], == 0 mod the others."""
    Q = 1
    for q in moduli:
        Q *= q
    out = []
    for q in moduli:
        R = Q // q
        out.append((R * pow(R, -1, q)) % Q)
    return out


def crt(residues: Sequence[int], moduli: Sequence[int]) -> int:
    Q = 1
    for q in moduli:
        Q *= q
    return sum(r * e for r, e in zip(residues, crt_idempotents(moduli))) % Q


# --------------------------------------------------------------------------
# the offset tau and its CRT / sign structure
# --------------------------------------------------------------------------
def lock_holds(a: int, b: int, sigma: int, N: int, M: int) -> bool:
    return (b * sigma * sigma - a * N) % M == 0


def tau(a: int, b: int, sigma: int, N: int, M: int) -> int:
    """tau(sigma) = a N sigma^{-1} + b sigma mod M^2 (any integer rep of sigma)."""
    M2 = M * M
    return (a * N * pow(sigma, -1, M2) + b * sigma) % M2


def gamma_of(sigma: int, N: int, M: int) -> int:
    """gamma = N / sigma^2 mod M: the HH lattice is b == gamma a (mod M), and it
    depends on sigma only through sigma^2."""
    return (N * pow(sigma * sigma, -1, M)) % M


def square_class_roots(gamma: int, N: int, primes: Sequence[int]) -> List[List[int]]:
    """For each prime m_i | M, the roots sigma_i mod m_i of sigma^2 == N/gamma.
    Returns [] for a prime where N/gamma is a non-residue (no units in class)."""
    out = []
    for m in primes:
        c = (N * pow(gamma, -1, m)) % m
        out.append(sqrt_mod_prime(c, m))
    return out


def sign_decomposition(a: int, b: int, roots: Sequence[int], primes: Sequence[int],
                       N: int) -> Dict:
    """Given one root r_i mod each prime m_i (all satisfying the lock mod m_i),
    return Delta_0 (part of primes with r == -r, i.e. m = 2) and Delta_i for odd
    primes such that

        tau(sigma(eps)) == Delta_0 + sum_i eps_i Delta_i   (mod M^2),

    sigma(eps) = CRT(eps_i r_i).  Delta's are integers in [0, M^2)."""
    sq = [m * m for m in primes]
    E = crt_idempotents(sq)
    M = 1
    for m in primes:
        M *= m
    M2 = M * M
    d0, deltas, odd = 0, [], []
    for m, r, e in zip(primes, roots, E):
        ti = tau(a, b, r, N, m)               # tau_i(r_i) mod m_i^2
        if m == 2:
            d0 = (d0 + e * ti) % M2
        else:
            deltas.append((e * ti) % M2)
            odd.append(m)
    return {"Delta0": d0, "Delta": deltas, "odd_primes": odd, "M": M, "M2": M2}


def sigma_of_eps(eps: Sequence[int], roots: Sequence[int], primes: Sequence[int]) -> int:
    """sigma(eps) = CRT(eps_i r_i); entries of eps for m = 2 are ignored."""
    res, k = [], 0
    for m, r in zip(primes, roots):
        if m == 2:
            res.append(r % 2)
        else:
            res.append((eps[k] * r) % m)
            k += 1
    return crt(res, primes)


def carry(eps: Sequence[int], dec: Dict) -> int:
    """c(eps) with tau(sigma(eps)) = Delta0 + sum eps_i Delta_i - M^2 c(eps),
    tau in [0, M^2)."""
    s = dec["Delta0"] + sum(e * d for e, d in zip(eps, dec["Delta"]))
    return s // dec["M2"]


def carry_range(dec: Dict) -> Tuple[int, int]:
    s_lo = dec["Delta0"] - sum(dec["Delta"])
    s_hi = dec["Delta0"] + sum(dec["Delta"])
    return s_lo // dec["M2"], s_hi // dec["M2"]


def all_eps(t: int):
    return product((1, -1), repeat=t)
