"""Quadratic character columns (notes/lv6_results.md §2.3): first-degree primes (r, s) of K, sampled
UNIFORMLY among first-degree primes with norm r in (x/2, x], and chi_(r,s)(a - b alpha) = ((a - b s)/r).

The theorem uses log x ~ C_d (log N)^2 (Lagarias-Odlyzko range); the toy uses a fixed bit size (disclosed
in the pre-registration).  Primality: sympy.isprime (BPSW; the theorem uses AKS).  Operational.
"""
from __future__ import annotations

import random
from typing import List, Sequence, Tuple

import flint
import sympy


def sample_first_degree_prime(f: Sequence[int], xbits: int, rng: random.Random) -> Tuple[int, int, int]:
    """Uniform first-degree prime (r, s) of K = Q(alpha) with r in (2^(xbits-1), 2^xbits].
    Rejection sampler: r uniform; if r prime, accept with prob n_r/d and return a uniform root s of f mod r.
    Returns (r, s, trials)."""
    d = len(f) - 1
    lo, hi = 2 ** (xbits - 1) + 1, 2 ** xbits
    trials = 0
    while True:
        trials += 1
        r = rng.randrange(lo, hi + 1)
        if not sympy.isprime(r):
            continue
        ctx = flint.fmpz_mod_poly_ctx(r)
        P = ctx([c % r for c in f])
        roots = [int(z) for z, mult in P.roots() if mult == 1]
        if rng.randrange(d) < len(roots):
            return r, roots[rng.randrange(len(roots))], trials


def chi(r: int, s: int, a: int, b: int) -> int:
    """Value in {0, 1} (additive notation: 1 means Legendre symbol -1) of chi_(r,s)(a - b alpha)."""
    v = sympy.jacobi_symbol((a - b * s) % r, r)
    if v == 0:
        raise ValueError("prime divides a - b s")
    return 0 if v == 1 else 1


def chi_elem_value(r: int, s: int, fd: int, elem: Sequence[int]) -> int:
    """Legendre symbol (in {0,1}) of an element of Z[beta] at the prime (r, beta -> fd*s); diagnostic."""
    x = 0
    bs = fd * s % r
    for c in reversed(elem):
        x = (x * bs + c) % r
    v = sympy.jacobi_symbol(x, r)
    if v == 0:
        raise ValueError("prime divides element")
    return 0 if v == 1 else 1
