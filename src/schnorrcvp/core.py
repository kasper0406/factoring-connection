"""Core exact objects for Schnorr/Adleman prime-number lattices.  Agent `schnorrcvp`, phase 19.

Exact integer arithmetic throughout.  Nothing here factors N or uses its factors: the only
factorisations performed are trial divisions of residuals by the public factor base
(the smoothness test that defines a relation), which is part of the operational method.

Conventions
-----------
* Factor base P = (p_1, ..., p_n), the first n primes.  An exponent vector e in Z^n gives
      u = prod_{e_i > 0} p_i^{e_i},   v = prod_{e_i < 0} p_i^{-e_i},   w = u - v*N.
  u and v are coprime by construction.
* A *fac-relation* (Schnorr's term) is such a triple with w != 0 and |w| P-smooth.  It gives the
  congruence u = w (mod N) between two P-smooth integers (sign of w recorded separately).
* The *relation lattice* Lambda_N = { (a, s) in Z^n x Z : (-1)^s prod p_i^{a_i} = 1 (mod N) }.
  A fac-relation (u, v, w) corresponds to the vector exps(u) - exps(w) (with sign bit) in Lambda_N;
  its weighted l1 norm sum |a_i| ln p_i equals ln(u*|w|) / gcd effects (exactly ln(u|w|) when
  gcd(u, w) = 1).
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Dict, List, Optional, Sequence, Tuple


def first_primes(n: int) -> List[int]:
    """The first n primes (exact, by a simple sieve that grows until n primes are found)."""
    if n <= 0:
        return []
    limit = max(16, int(n * 1.3 * (max(2, n).bit_length())) + 10)
    while True:
        sieve = bytearray([1]) * (limit + 1)
        sieve[0:2] = b"\x00\x00"
        for i in range(2, isqrt(limit) + 1):
            if sieve[i]:
                sieve[i * i :: i] = bytearray(len(sieve[i * i :: i]))
        ps = [i for i in range(limit + 1) if sieve[i]]
        if len(ps) >= n:
            return ps[:n]
        limit *= 2


def primes_upto(B: int) -> List[int]:
    if B < 2:
        return []
    sieve = bytearray([1]) * (B + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, isqrt(B) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(sieve[i * i :: i]))
    return [i for i in range(B + 1) if sieve[i]]


def factor_over(x: int, P: Sequence[int]) -> Optional[Dict[int, int]]:
    """Exponents of |x| over P if |x| is P-smooth (x != 0), else None.  Trial division only."""
    if x == 0:
        return None
    x = abs(x)
    out: Dict[int, int] = {}
    for p in P:
        if x == 1:
            break
        if x % p == 0:
            k = 0
            while x % p == 0:
                x //= p
                k += 1
            out[p] = k
    return out if x == 1 else None


def is_smooth(x: int, P: Sequence[int]) -> bool:
    return factor_over(x, P) is not None


def uvw(e: Sequence[int], P: Sequence[int], N: int) -> Tuple[int, int, int]:
    """(u, v, w) of an exponent vector, exactly."""
    u = 1
    v = 1
    for ei, p in zip(e, P):
        if ei > 0:
            u *= p ** ei
        elif ei < 0:
            v *= p ** (-ei)
    return u, v, u - v * N


def is_fac_relation(e: Sequence[int], P: Sequence[int], N: int) -> bool:
    u, v, w = uvw(e, P, N)
    return w != 0 and is_smooth(w, P)


def relation_key(u: int, w: int) -> Tuple[int, int]:
    """Canonical key for the congruence u = w (mod N) (used to count distinct relations).

    Two triples with the same (u, w) give the same congruence; (u, w) and (w', u') with
    u'=|w|,... are distinct congruences in general, so we keep the ordered pair with sign."""
    return (u, w)


def lambda_vector(u: int, w: int, P: Sequence[int]) -> Optional[Tuple[int, ...]]:
    """Vector of Lambda_N for the congruence u = w (mod N): (exps(u) - exps(|w|), sign bit).

    Sign bit s = 1 iff w < 0 (so that (-1)^s * prod p^{a} = u / w = 1 mod N).  None if not smooth."""
    fu = factor_over(u, P)
    fw = factor_over(w, P)
    if fu is None or fw is None:
        return None
    a = tuple(fu.get(p, 0) - fw.get(p, 0) for p in P)
    return a + ((1 if w < 0 else 0),)


def in_lambda(vec: Sequence[int], P: Sequence[int], N: int) -> bool:
    """Exact membership test for Lambda_N (gcd checks first: p | N would make p non-invertible)."""
    num = 1
    den = 1
    for ai, p in zip(vec[:-1], P):
        if ai > 0:
            num = num * pow(p, ai, N) % N
        elif ai < 0:
            den = den * pow(p, -ai, N) % N
    if gcd(den, N) != 1 or gcd(num, N) != 1:
        return False
    s = vec[-1] % 2
    lhs = (num * (-1 if s else 1)) % N
    return lhs == den % N


def uw_lower_bound_holds(u: int, v: int, w: int, N: int) -> bool:
    """resbarrier §3(h): for integers u, v, w != 0 with u - vN = w and u > 0, v > 0: u*|w| >= N - 1.

    (Proof: if w > 0 then u = vN + w > N so u*w > N.  If w < 0 then u = vN - |w| >= N - |w| > 0, so
    u*|w| >= (N - |w|)*|w| >= N - 1 because 1 <= |w| <= N - 1.)"""
    assert u - v * N == w and u > 0 and v > 0 and w != 0
    return u * abs(w) >= N - 1


def squares_dependency_factor(rel_vectors: List[Tuple[int, ...]], P: Sequence[int], N: int,
                              max_deps: int = 64) -> Tuple[Optional[int], Dict[str, int]]:
    """Linear algebra mod 2 over Lambda_N vectors, then square roots of 1 mod N and gcds.

    rel_vectors are Lambda_N vectors (a, s) with (-1)^s prod p^a = 1 mod N.  A dependency mod 2
    gives y = (sum)/2 with x = (-1)^{s/2} prod p^{y} and x^2 = 1 mod N; gcd(x - 1, N) is tested.
    Returns (factor or None, counters).  Exact (Gaussian elimination over GF(2) on Python ints)."""
    m = len(rel_vectors)
    dim = len(P) + 1
    rows = []
    for i, vec in enumerate(rel_vectors):
        bits = 0
        for j, a in enumerate(vec):
            if a & 1:
                bits |= 1 << j
        rows.append((bits, 1 << i))
    # elimination
    pivots: Dict[int, Tuple[int, int]] = {}
    deps: List[int] = []
    for bits, comb in rows:
        while bits:
            hb = bits.bit_length() - 1
            if hb in pivots:
                pb, pc = pivots[hb]
                bits ^= pb
                comb ^= pc
            else:
                pivots[hb] = (bits, comb)
                break
        if bits == 0:
            deps.append(comb)
    ctr = {"relations": m, "dependencies": len(deps), "gcd_calls": 0, "trivial_roots": 0}
    for comb in deps[:max_deps]:
        tot = [0] * dim
        for i in range(m):
            if comb >> i & 1:
                for j, a in enumerate(rel_vectors[i]):
                    tot[j] += a
        assert all(t % 2 == 0 for t in tot)
        y = [t // 2 for t in tot]
        num, den = 1, 1
        for yj, p in zip(y[:-1], P):
            if yj > 0:
                num = num * pow(p, yj, N) % N
            elif yj < 0:
                den = den * pow(p, -yj, N) % N
        g = gcd(den, N)
        ctr["gcd_calls"] += 1
        if 1 < g < N:
            return g, ctr
        x = num * pow(den, -1, N) % N
        if y[-1] % 2:
            x = (-x) % N
        assert (x * x - 1) % N == 0
        g = gcd(x - 1, N)
        ctr["gcd_calls"] += 1
        if 1 < g < N:
            return g, ctr
        ctr["trivial_roots"] += 1
    return None, ctr
