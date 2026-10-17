"""Track A (b): closed-form section evaluation in quotient algebras.  OPERATIONAL.

All functions take only N and public parameters.  Sections are of
F_N = (1+X)^N - 1 - X^N:  S_{r,M}(Y) = sum_j a_{r+jM} Y^j.

* ``section_values(N, M, c)``: [S_{r,M}(c) for r in range(M)] from one
  computation of (1+X)^N mod (N, X^M - c)  (rank-M algebra, ~1.5 log2 N
  products of degree-<M polynomials).
* ``section_values_ext(N, M, s, a)``: S_{r,M}(Z) mod (N, Z^s - a) for all r, as
  an M x s table, from (1+X)^N mod (X^{Ms} - a)  (rank Ms; AKS-style).
* ``descent_split(N, M, c, s)``: the descent identity
  X^{2M} - s^2 = (X^M - s)(X^M + s): children values at c = s^2 from parent
  values at +s and -s (requires 2 to be a unit; 2s checked by gcd).
"""
from __future__ import annotations

from typing import List, Optional

from ..common.metrics import Metrics
from ..common.ring import FlintQuotient, inv_mod, power_of_linear


def _endpoint_fix(N: int, M: int, c: int, vals: List[int]) -> List[int]:
    """Subtract the endpoint monomials 1 and X^N of (1+X)^N (reduced mod X^M - c)."""
    vals = list(vals)
    vals[0] = (vals[0] - 1) % N
    vals[N % M] = (vals[N % M] - pow(c, N // M, N)) % N
    return vals


def section_values(N: int, M: int, c: int, m: Optional[Metrics] = None) -> List[int]:
    """S_{r,M}(c) for r = 0..M-1 (sections of F_N), via R_N[X]/(X^M - c)."""
    c %= N
    if M == 1:
        # R_N[X]/(X - c): value F_N(c)
        v = (pow(1 + c, N, N) - 1 - pow(c, N, N)) % N
        if m is not None:
            m.inc("q_mul", 2 * N.bit_length())
            m.inc("q_mul_dimsum", 2 * N.bit_length())
        return [v]
    Q = FlintQuotient(N, [(-c) % N] + [0] * (M - 1) + [1], m)
    vals = Q.coeffs(power_of_linear(Q, N))
    return _endpoint_fix(N, M, c, vals)


def section_values_ext(N: int, M: int, s: int, a: int, m: Optional[Metrics] = None) -> List[List[int]]:
    """Table T[r][t] = [Z^t] (S_{r,M}(Z) mod (N, Z^s - a)), r < M, t < s.

    Uses X^M -> Z, so (1+X)^N is reduced modulo X^{Ms} - a (rank M*s)."""
    L = M * s
    a %= N
    Q = FlintQuotient(N, [(-a) % N] + [0] * (L - 1) + [1], m)
    v = Q.coeffs(power_of_linear(Q, N))
    v = _endpoint_fix(N, L, a, v)
    return [[v[r + t * M] for t in range(s)] for r in range(M)]


def descent_split(N: int, M: int, s: int, m: Optional[Metrics] = None):
    """Children of every parent class r mod M, evaluated at c = s^2, computed from
    parent evaluations at +s and -s:
        S_{r,2M}(s^2)   = (P_r(s) + P_r(-s)) / 2
        S_{r+M,2M}(s^2) = (P_r(s) - P_r(-s)) / (2 s)
    where P_r = S_{r,M}.  Raises NonUnit if 2s is not a unit (then gcd is a
    factor or s = 0).  Total algebra rank used: M + M = 2M.
    """
    plus = section_values(N, M, s, m)
    minus = section_values(N, M, (-s) % N, m)
    i2 = inv_mod(2, N, m)
    i2s = inv_mod(2 * s, N, m)
    even = [(x + y) * i2 % N for x, y in zip(plus, minus)]
    odd = [(x - y) * i2s % N for x, y in zip(plus, minus)]
    return even, odd
