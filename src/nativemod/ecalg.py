"""nativemod OPERATIONAL: the Edixhoven-Couveignes / Schoof algebra at a composite index, as a PRIM program.

A = (Z/N)[x]/(g), g a FIXED monic integer polynomial (here the monic division polynomial psi_l / l of a fixed
curve; its coefficients are public configuration).  At p the components of A are F_{p^f}; the literal map
z -> z^N equals Frob_p o Pow_q there (not a ring map), so every read-out below is a PRIM program with an
exponent gate in extension components followed by a coordinate (trace) functional or a Jacobi symbol: the
fragment that leakmap §3.5 / prim §7 leave OPEN.

Inputs: N and g only.  Every Z/N value used in a Jacobi symbol is gcd-checked first (rule 3).
"""
from __future__ import annotations

from math import gcd
from typing import Dict, List, Optional, Sequence, Tuple

import gmpy2


class Split(Exception):
    def __init__(self, d: int, where: str):
        super().__init__(f"split {d} at {where}")
        self.d = d
        self.where = where


def _mulmod_poly(a: Sequence[int], b: Sequence[int], g: Sequence[int], N: int) -> List[int]:
    """(a*b) mod (g, N); g monic of degree D given low->high with g[D] = 1."""
    D = len(g) - 1
    prod = [0] * (2 * D - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                prod[i + j] += ai * bj
    for k in range(2 * D - 2, D - 1, -1):
        c = prod[k] % N
        if c:
            for i in range(D):
                prod[k - D + i] -= c * g[i]
        prod[k] = 0
    return [x % N for x in prod[:D]]


def powmod_x(e: int, g: Sequence[int], N: int) -> List[int]:
    D = len(g) - 1
    result = [1] + [0] * (D - 1)
    base = [0, 1] + [0] * (D - 2) if D >= 2 else [0]
    if D == 1:
        base = [(-g[0]) % N]
    for bit in bin(e)[2:]:
        result = _mulmod_poly(result, result, g, N)
        if bit == "1":
            # multiply by x: shift
            result = _mulmod_poly(result, base, g, N)
    return result


def power_sums(g: Sequence[int], count: int, N: Optional[int] = None) -> List[int]:
    """Power sums s_0..s_{count-1} of the roots of the monic polynomial g (Newton's identities, which are
    polynomial identities over Z, so they hold over Z/N: s_k = Tr_{A/(Z/N)}(x^k)).  Reduced mod N if given."""
    D = len(g) - 1
    # e_k with g = x^D + c_{D-1} x^{D-1} + ... : c_{D-k} = (-1)^k e_k
    c = list(g)
    s = [D]
    for k in range(1, count):
        # s_k + c_{D-1} s_{k-1} + ... + c_{D-k+1} s_1 + k c_{D-k} = 0  (k <= D);  for k > D no k c term
        tot = 0
        for i in range(1, min(k, D) + 1):
            if i < k:
                tot += c[D - i] * s[k - i]
        if k <= D:
            tot += k * c[D - k]
        s.append(-tot if N is None else (-tot) % N)
    return s


def trace(z: Sequence[int], s: Sequence[int], N: int) -> int:
    return sum(zi * si for zi, si in zip(z, s)) % N


def _checked_jacobi(v: int, N: int, where: str) -> int:
    d = gcd(v, N)
    if d != 1:
        if d != N:
            raise Split(d, where)
        return 0
    return int(gmpy2.jacobi(v, N))


def statistics(N: int, g: Sequence[int], s: Optional[Sequence[int]] = None) -> Dict[str, int]:
    """Read-outs of y = x^N in A = (Z/N)[x]/(g): S1 = (Tr y | N), S2 = (Tr xy | N), S3 = (Tr y^2 | N),
    S4 = [g(y) = 0 in A].  Raises Split on any proper gcd."""
    D = len(g) - 1
    gN = [c % N for c in g]
    if s is None:
        s = power_sums(gN, D, N)
    y = powmod_x(N, gN, N)
    x = [0, 1] + [0] * (D - 2)
    out = {}
    out["S1"] = _checked_jacobi(trace(y, s, N), N, "S1")
    out["S2"] = _checked_jacobi(trace(_mulmod_poly(x, y, gN, N), s, N), N, "S2")
    out["S3"] = _checked_jacobi(trace(_mulmod_poly(y, y, gN, N), s, N), N, "S3")
    # S4: g(y) via Horner in A
    acc = [0] * D
    for c in reversed(gN):
        acc = _mulmod_poly(acc, y, gN, N)
        acc[0] = (acc[0] + c) % N
    nz = [a for a in acc if a]
    if not nz:
        out["S4"] = 1
    else:
        gg = 0
        for a in nz:
            gg = gcd(gg, a)
        dd = gcd(gg, N)
        if 1 < dd < N:
            raise Split(dd, "S4")
        out["S4"] = 0
    out["TrY"] = trace(y, s, N)
    return out


def monic_division_poly(psi: Sequence[int], l: int, N: int) -> List[int]:
    """psi_l (low->high, leading coefficient l) made monic mod N; a proper gcd(l, N) is a split."""
    d = gcd(l, N)
    if 1 < d < N:
        raise Split(d, "monic")
    inv = pow(l, -1, N)
    return [(c * inv) % N for c in psi]
