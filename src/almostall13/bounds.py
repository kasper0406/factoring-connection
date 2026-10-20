"""Exact helper computations for notes/almostall13_results.md.  DIAGNOSTIC / analysis only.

- Lemma M1-AA: divisors of R in a dyadic window and Rankin's sigma = 1/2 bound (exact lower bound of the bound).
- Lemma MP: moment information with relative error eta cannot certify Pr[X = 0] < eta (mixture with a point mass).
- Theorem NR-AA: the cost constant (floats for display only, as in src/nested/barrier.nr_constant).
"""
from __future__ import annotations

from fractions import Fraction
from math import isqrt
from typing import Dict, List, Sequence

import sympy


def window_divisors(R: int, X: int) -> int:
    """HARNESS: number of divisors of |R| in [X, 2X) (factors R with sympy; never used by an algorithm)."""
    return sum(1 for n in sympy.divisors(abs(R)) if X <= n < 2 * X)


def rankin_half_lower(primes: Sequence[int], X: int) -> Fraction:
    """An exact LOWER bound for Rankin's bound (2X)^{1/2} * prod_{p | R} (1 - p^{-1/2})^{-1}.

    Uses (2X)^{1/2} >= isqrt(2X) and (1 - p^{-1/2})^{-1} >= 1 + p^{-1/2} >= 1 + 1/(isqrt(p) + 1).
    Since #divisors of R in [X, 2X) <= sum_{n | R} (2X/n)^{1/2} = Rankin bound, any count exceeding this lower bound
    would refute the code, not the lemma; a count below it is consistent with the lemma.
    """
    v = Fraction(isqrt(2 * X))
    for p in set(primes):
        v *= 1 + Fraction(1, isqrt(p) + 1)
    return v


def rankin_half_exact_holds(R: int, X: int) -> bool:
    """Checks #div_[X,2X)(R) <= sum_{n | R} sqrt(2X/n) exactly (compare squares termwise is not needed: each
    divisor n in the window contributes sqrt(2X/n) >= 1, and the sum over all divisors only adds nonnegative terms)."""
    cnt = window_divisors(R, X)
    lower = 0
    for n in sympy.divisors(abs(R)):
        if n <= 2 * X:
            lower += isqrt((2 * X) // n)          # floor(sqrt(2X/n)) <= sqrt(2X/n); >= 1 for n <= 2X
    return cnt <= lower


def mp_mixture(law: Dict[int, Fraction], eta: Fraction, kmax: int) -> Dict:
    """Lemma MP: X' = 0 with probability eta, else X ~ law.  Returns Pr[X' = 0] and the relative moment changes
    E[X'^k]/E[X^k] - 1 = -eta for every k >= 1 (exact)."""
    assert sum(law.values()) == 1
    mom = [sum(p * Fraction(x) ** k for x, p in law.items()) for k in range(1, kmax + 1)]
    mom2 = [(1 - eta) * m for m in mom]
    rel = [m2 / m - 1 for m, m2 in zip(mom, mom2)]
    p0 = eta + (1 - eta) * law.get(0, Fraction(0))
    return {"p0": p0, "rel_changes": rel}


def nraa_constant(b: float, b1: float, s: float, delta: float, c_chi: float) -> float:
    """c_NR-AA = max(2b, b + c_chi*delta, (b + 1/b1)/(1 - 2^{-1/3}), (b + 2s)/(1 - (delta b1)^{1/3})).

    c_chi is the (absolute, uncomputed) constant of the unconditional factorial-bit character primes
    (lv5unif section 9 at d = delta D'); floats for display only."""
    tau = delta * b1
    assert 0 < tau < 1
    return max(2 * b, b + c_chi * delta, (b + 1 / b1) / (1 - 2 ** (-1 / 3)), (b + 2 * s) / (1 - tau ** (1 / 3)))


def m1aa_pairs_log2(d: int, Hf: int, Hg: int) -> Fraction:
    """log2 of an upper bound for the number of pairs (f, g) with ||f|| <= Hf (degree d) and ||g|| <= Hg, as an
    exact rational upper bound via bit lengths: (d+1)*bitlen(2Hf+1) + 2*bitlen(2Hg+1)."""
    return Fraction((d + 1) * (2 * Hf + 1).bit_length() + 2 * (2 * Hg + 1).bit_length())


def primes_of(R: int) -> List[int]:
    """HARNESS: prime divisors of |R| (sympy)."""
    return list(sympy.factorint(abs(R)).keys())
