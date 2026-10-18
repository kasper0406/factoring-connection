"""OPERATIONAL ECM prefilter over Z/N (sees only N and public parameters).

Families (fixed, explicit, n-independent):
  'suyama'  sigma = 6, 7, 8, ...:  u = sigma^2 - 5, v = 4 sigma, x0 = u^3, z0 = v^3,
            a24 = (A+2)/4 = (v-u)^3 (3u+v) / (16 u^3 v)      (Montgomery, x-only)
  'mont'    control: A_j = j + 3, x0 = 2, z0 = 1, a24 = (A+2)/4

Stage 1 with STAGED gcds: the scalar M_B = prod_{l <= B} l^{e_l} (e_l = max e with
l^e <= X) is applied one prime power at a time, in increasing order of l, and gcd(Z, N) is
taken after each step (step 0 = the curve-parameter / degeneracy gcd).  A curve separates
p and q iff the first step at which Z vanishes mod p differs from the one mod q
(notes/ecmavg_results.md, Lemma S).

Parameters (n-only): b_N = ceil(n/2); B = floor(2^{b_N * beta}) for beta in {1/2, 1/3, 1/4};
X = 2^{ceil((n+1)/2)} + 2^{ceil((n+1)/4)+1} + 1  (>= q + 1 + 2 sqrt q for every q <= sqrt(2N)).
"""
from __future__ import annotations

from fractions import Fraction as F
from math import gcd, isqrt
from typing import Dict, Iterator, List, Optional, Tuple

from gmpy2 import iroot
from sympy import primerange  # prime list up to B (public, factor-free)


def stage_steps(B: int, X: int) -> List[Tuple[int, int]]:
    """[(l, l^e)] for primes l <= B, e = max e >= 1 with l^e <= X, increasing l."""
    out = []
    for l in primerange(2, B + 1):
        pe = l
        while pe * l <= X:
            pe *= l
        out.append((l, pe))
    return out


def bound_B(b: int, beta: F) -> int:
    """B = floor(2^{b * beta}) exactly (beta = num/den)."""
    beta = F(beta)
    return max(2, int(iroot(2 ** (b * beta.numerator), beta.denominator)[0]))


def params_for_n(n: int, beta: F) -> Tuple[int, int]:
    """(B, X) for an n-bit input: B = floor(2^{ceil(n/2) beta})."""
    bN = (n + 1) // 2
    X = 2 ** (-(-(n + 1) // 2)) + 2 ** (-(-(n + 1) // 4) + 1) + 1
    return bound_B(bN, beta), X


def curve_params_int(family: str, j: int) -> Tuple[Tuple[int, int], Tuple[int, int], int]:
    """Integer data of curve j: ((x0, z0), (a24num, a24den), degeneracy integer D).
    The curve is degenerate mod p iff p | D (singular, or undefined parameters)."""
    if family == "suyama":
        s = 6 + j
        u, v = s * s - 5, 4 * s
        num, den = (v - u) ** 3 * (3 * u + v), 16 * u ** 3 * v
        # A + 2 = 4 num/den; A - 2 = 4 (num - 4 den/4 ...): A-2 = 4(num - den)/den
        D = den * (v - u) * (3 * u + v) * (num - den)
        return (u ** 3, v ** 3), (num, den), D
    if family == "mont":
        A = j + 3
        num, den = A + 2, 4
        D = 2 * (A - 2) * (A + 2)
        return (2, 1), (num, den), D
    raise ValueError(family)


def _dbl(X, Z, a24, N):
    t1 = (X + Z) * (X + Z) % N
    t2 = (X - Z) * (X - Z) % N
    t3 = (t1 - t2) % N
    return t1 * t2 % N, t3 * (t2 + a24 * t3) % N


def _add(X1, Z1, X2, Z2, Xd, Zd, N):
    u = (X1 - Z1) * (X2 + Z2) % N
    v = (X1 + Z1) * (X2 - Z2) % N
    return Zd * (u + v) * (u + v) % N, Xd * (u - v) * (u - v) % N


def ladder(s: int, X: int, Z: int, a24: int, N: int, ops: Optional[Dict[str, int]] = None):
    """x-only Montgomery ladder [s](X:Z) mod N (s >= 1)."""
    R0, R1 = (X, Z), _dbl(X, Z, a24, N)
    for bit in bin(s)[3:]:
        if bit == "1":
            R0, R1 = _add(*R0, *R1, X, Z, N), _dbl(*R1, a24, N)
        else:
            R0, R1 = _dbl(*R0, a24, N), _add(*R0, *R1, X, Z, N)
        if ops is not None:
            ops["mul"] = ops.get("mul", 0) + 11
    return R0


def ecm_prefilter(N: int, J: int, beta: F = F(1, 3), family: str = "suyama",
                  ops: Optional[Dict[str, int]] = None) -> Tuple[Optional[int], int]:
    """Run curves 0..J-1 sequentially; stop at the first proper divisor.
    Returns (divisor or None, number of curves run).  Every gcd outcome is checked."""
    n = N.bit_length()
    B, Xb = params_for_n(n, F(beta))
    steps = stage_steps(B, Xb)
    if ops is None:
        ops = {}
    for j in range(J):
        (x0, z0), (num, den), D = curve_params_int(family, j)
        g = gcd(D % N, N)
        ops["gcd"] = ops.get("gcd", 0) + 1
        if 1 < g < N:
            return g, j + 1
        if g == N:
            continue
        a24 = num * pow(den, -1, N) % N        # den invertible: den | D and gcd(D, N) = 1
        X, Z = x0 % N, z0 % N
        for (_l, pe) in steps:
            X, Z = ladder(pe, X, Z, a24, N, ops)
            g = gcd(Z, N)
            ops["gcd"] = ops.get("gcd", 0) + 1
            if 1 < g < N:
                return g, j + 1
            if g == N:
                break
    return None, J
