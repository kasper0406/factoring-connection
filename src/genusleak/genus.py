"""genusleak genus (OPERATIONAL reductions; oracles are supplied by the caller).

- `ambiguous_reduced_forms(D)`: all primitive reduced ambiguous forms of D < 0, by divisor enumeration of 4kN.
  EXPONENTIAL (needs divisors of N): used only in tests/harness as ground truth for Lemma 0.
- `qr_via_principal_genus(N, k, a, oracle, rng)`: Prop GEN, QR_N <= GEN-SQ.  Given a with (a/N) = 1, build a
  prime m == a (mod N), m == 1 (mod 8k), and the prime form f_m of D = -4kN; then
  [f_m is a square in Cl(D)] == [(a/p) = 1] == [a is a square mod N].
  Uses a public primality test (gmpy2.is_prime, BPSW) on m.
"""
from __future__ import annotations

import random
from math import gcd, isqrt
from typing import Callable, List

import gmpy2

from .forms import Form, is_reduced, reduce_form, sqrt_mod_prime

__all__ = ["ambiguous_reduced_forms", "qr_via_principal_genus"]


def ambiguous_reduced_forms(D: int) -> List[Form]:
    """All primitive reduced forms (a, b, c) of discriminant D < 0 with b = 0, a = b or a = c (exhaustive)."""
    M = -D  # = 4ac - b^2
    out = set()
    for a in range(1, isqrt(M) + 1):
        # b = 0: 4ac = M
        if M % (4 * a) == 0:
            out.add((a, 0, M // (4 * a)))
        # b = a: 4ac - a^2 = M  -> c = (M + a^2)/(4a)
        if (M + a * a) % (4 * a) == 0:
            out.add((a, a, (M + a * a) // (4 * a)))
        # c = a: 4a^2 - b^2 = M, 0 < b < a
        r = 4 * a * a - M
        if r > 0:
            b = isqrt(r)
            if b * b == r and 0 < b < a:
                out.add((a, b, a))
    return sorted(f for f in out if is_reduced(f) and gcd(gcd(f[0], f[1]), f[2]) == 1)


def qr_via_principal_genus(N: int, k: int, a: int, oracle: Callable[[Form], bool], rng: random.Random,
                           max_tries: int = 100000):
    """Return (answer, m, tries): answer = True iff a is a QR mod N (given (a/N) = 1), using one call to the
    principal-genus oracle.  Returns a factor instead if a gcd with N is nontrivial."""
    g = gcd(a, N)
    if 1 < g < N:
        return ("factor", g, 0)
    assert gmpy2.jacobi(a, N) == 1
    D = -4 * k * N
    mod = 8 * k
    # m0 == a (mod N), m0 == 1 (mod 8k)
    m0 = (a + N * ((1 - a) * pow(N, -1, mod) % mod)) % (N * mod)
    for tries in range(1, max_tries + 1):
        m = m0 + N * mod * rng.randrange(1, N)
        if gmpy2.is_prime(m):
            # (D/m) = 1 automatically (see docstring); build the prime form of norm m
            r = sqrt_mod_prime(D, m)
            if r % 2:
                r = m - r
            assert (r * r - D) % (4 * m) == 0
            f = reduce_form((m, r, (r * r - D) // (4 * m)))
            return (oracle(f), m, tries)
    raise RuntimeError("no prime found")
