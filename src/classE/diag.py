"""DIAGNOSTIC (knows the factors).  Never imported by operational code (src/classE/{structure,fd,battery}.py).

* `vis_predicate(a, p, q)`: Prop VIS necessary condition -- a = +-r'/s' (mod q) with r's' < q/p.
* `planted_unit(p, q, a_p, a_q)`: the unit with given CRT components.
* `gar_phase(x, p, q)`: Lemma GAR (leakmap section 4.2): rep(x)/N = (x_p - x_q) theta + x_q/N (mod 1).
* `sqrt_uv(N, p, q)`: Prop SQ-CF: N mod r = -uv (mod r), u = r - p, v = q - r.
* `dedekind_freq_local(a, p, q)`, `dedekind_space_local(a, p, q)`: the two "local parts" of s(a, N).
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd, isqrt
from typing import Tuple

from src.newfamily.recip import dedekind6


def planted_unit(p: int, q: int, a_p: int, a_q: int) -> int:
    N = p * q
    return (a_p * q * pow(q, -1, p) + a_q * p * pow(p, -1, q)) % N


def vis_predicate(a: int, p: int, q: int) -> bool:
    """True iff a = +-r/s (mod q) for some positive integers r, s with r s < q/p (Prop VIS necessity side)."""
    bound = q / p  # harness only; exact comparison below
    for s in range(1, int(bound) + 2):
        for r in range(1, int(bound) + 2):
            if r * s * p < q:
                for sg in (1, -1):
                    if (a * s - sg * r) % q == 0:
                        return True
    return False


def gar_phase(x: int, p: int, q: int) -> Tuple[Fraction, Fraction]:
    """(lhs, rhs) of Lemma GAR mod 1, as exact rationals in [0, 1)."""
    N = p * q
    xp, xq = x % p, x % q
    theta = Fraction(pow(q, -1, p), p)
    lhs = Fraction(x % N, N)
    rhs = (xp - xq) * theta + Fraction(xq, N)
    return lhs, rhs - (rhs.numerator // rhs.denominator)


def sqrt_uv(N: int, p: int, q: int) -> Tuple[int, int, int]:
    r = isqrt(N)
    return r, r - p, q - r


def dedekind_space_local(a: int, p: int, q: int) -> Fraction:
    """s(a,p) + s(a,q): the spatial (sawtooth-block) local part (newfamily Prop DS1/DSF)."""
    return Fraction(dedekind6(a % p, p), 6 * p) + Fraction(dedekind6(a % q, q), 6 * q)


def dedekind_freq_local(a: int, p: int, q: int) -> Fraction:
    """s(a,p)/q + s(a,q)/p: the frequency (cotangent-block) local part."""
    return Fraction(dedekind6(a % p, p), 6 * p * q) + Fraction(dedekind6(a % q, q), 6 * q * p)
