"""Fourier-Dedekind (lattice-cone) sums, exact.  OPERATIONAL: functions see only the integers they are given.

    s_n(a_1..a_d; N) = (1/N) * sum_{k=1}^{N-1} zeta^{k n} / prod_j (1 - zeta^{k a_j}),   zeta = e(1/N),

for a_j units mod N (Beck-Diaz-Robins; Beck-Robins, "Computing the Continuous Discretely", ch. 8).

Exact evaluation (notes/classE_results.md section 3.1).  For x^N = 1, x != 1:  1/(1 - x) = -(1/N) sum_{r<N} r x^r.
Expanding the product and summing over k (k = 0 term removed) gives

    s_n(a; N) = (-1)^d / N^{d+1} * [ N * C - (N(N-1)/2)^d ],
    C = C_n(a; N) = sum_{r in [0,N)^d,  n + sum_j a_j r_j = 0 (mod N)} r_1 r_2 ... r_d.

So N^d s is an explicit integer combination; everything below computes the integer C exactly.

* d = 1: C = (-n a^{-1} mod N).
* d = 2: after normalising a_2 = 1 (k -> k a_2^{-1}),  C = sum_{r<N} r * ((c0 + c1 r) mod N), a weighted floor sum,
  computed in O(log N) by the Euclid-type recursion `floor_sums` (poly-time; this is the classical d = 2 /
  Dedekind-Rademacher case).
* d = 3: after normalising a_3 = 1,  C = sum_{r1} r1 * G((c0 + c1 r1) mod N) with
  G(c) = sum_{r2} r2 ((c - beta r2) mod N) obeying G(c+1) = G(c) + N(N-1)/2 - N (beta^{-1}(c+1) mod N).
  O(N) time and memory (numpy) -- TOY SIZES ONLY (N <= ~2^24).  A polynomial-time algorithm exists in principle
  (Barvinok, fixed dimension) but is not implemented here.
* `fd_brute` : O(N^{d-1}) direct sum over the hyperplane, any d (tests).
"""
from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import gcd
from typing import Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------------------------------------
# floor sums:  f = sum_{x<n} floor((a x + b)/m),  g = sum_{x<n} x floor((a x + b)/m),  h = sum_{x<n} floor(.)^2
# (classical "Euclid-like" recursion; all exact integers; a, b >= 0, m >= 1, n >= 0)
# ---------------------------------------------------------------------------------------------------------
def floor_sums(n: int, a: int, b: int, m: int) -> Tuple[int, int, int]:
    if n <= 0:
        return 0, 0, 0
    if a == 0:
        t = b // m
        return n * t, t * n * (n - 1) // 2, n * t * t
    if a >= m or b >= m:
        qa, ra = divmod(a, m)
        qb, rb = divmod(b, m)
        f1, g1, h1 = floor_sums(n, ra, rb, m)
        s1 = n * (n - 1) // 2               # sum x
        s2 = (n - 1) * n * (2 * n - 1) // 6  # sum x^2
        f = f1 + qa * s1 + qb * n
        g = g1 + qa * s2 + qb * s1
        h = (h1 + qa * qa * s2 + qb * qb * n + 2 * qa * qb * s1
             + 2 * qa * g1 + 2 * qb * f1)
        return f, g, h
    M = (a * (n - 1) + b) // m
    if M == 0:
        return 0, 0, 0
    # swap: count lattice points; standard identities with (m, a) roles exchanged
    f1, g1, h1 = floor_sums(M, m, m - b - 1, a)
    f = (n - 1) * M - f1
    g = (M * (n - 1) * n - h1 - f1) // 2
    h = (n - 1) * M * (M + 1) - 2 * g1 - 2 * f1 - f
    return f, g, h


def _sum_r_times_mod(N: int, c0: int, c1: int) -> int:
    """sum_{r=0}^{N-1} r * ((c0 + c1 r) mod N), c0, c1 in [0, N)."""
    s1 = N * (N - 1) // 2
    s2 = (N - 1) * N * (2 * N - 1) // 6
    _, g, _ = floor_sums(N, c1, c0, N)
    return c0 * s1 + c1 * s2 - N * g


# ---------------------------------------------------------------------------------------------------------
def _normalise(n: int, a: Sequence[int], N: int):
    """k -> k * a_last^{-1}: returns (m, [alpha_1..alpha_{d-1}]) with last argument 1."""
    for x in a:
        if gcd(x, N) != 1:
            raise ValueError("arguments must be units mod N")
    il = pow(a[-1], -1, N)
    return (n * il) % N, [(x * il) % N for x in a[:-1]]


def fd_C(n: int, a: Sequence[int], N: int) -> int:
    """The integer C_n(a; N) (definition in the module docstring)."""
    d = len(a)
    m, al = _normalise(n, a, N)
    if d == 1:
        return (-m) % N
    if d == 2:
        # r2 = (-m - alpha r1) mod N
        return _sum_r_times_mod(N, (-m) % N, (-al[0]) % N)
    if d == 3:
        return _C3_numpy(N, m, al[0], al[1])
    raise ValueError("d <= 3 (use fd_brute)")


def _C3_numpy(N: int, m: int, alpha: int, beta: int) -> int:
    """C = sum_{r1, r2} r1 r2 ((-m - alpha r1 - beta r2) mod N) in O(N) time, exact.

    G(c) := sum_{r2} r2 ((c - beta r2) mod N) = G(0) + c N(N-1)/2 - N cs(c),  cs(c) = sum_{j=1}^{c} (beta^{-1} j mod N),
    C = sum_{r1} r1 G(c(r1)),  c(r1) = (-m - alpha r1) mod N
      = G(0) S1 + (N(N-1)/2) sum r1 c(r1) - N sum r1 cs(c(r1)).
    int64 arrays hold values <= N^2 (N < 2^31); products are summed in exact Python integers chunk by chunk."""
    if N >= (1 << 31):
        raise ValueError("toy evaluator: N < 2^31")
    ib = pow(beta, -1, N)
    G0 = _sum_r_times_mod(N, 0, (-beta) % N)
    S1 = N * (N - 1) // 2
    r = np.arange(N, dtype=np.int64)
    w = (ib * np.arange(1, N, dtype=np.int64)) % N if N < (1 << 31) else None
    cs = np.empty(N, dtype=np.int64)
    cs[0] = 0
    np.cumsum(w, out=cs[1:])                                   # <= N^2 / 2
    c = (-m - alpha * r) % N if N < (1 << 31) else None        # alpha < N < 2^31: alpha * r < 2^62
    t1 = 0
    t2 = 0
    step = 1 << 14
    for s in range(0, N, step):
        rr = r[s:s + step].astype(object)
        t1 += int(np.dot(rr, c[s:s + step].astype(object)))
        t2 += int(np.dot(rr, cs[c[s:s + step]].astype(object)))
    return G0 * S1 + S1 * t1 - N * t2


def _C3_object(N: int, m: int, alpha: int, beta: int) -> int:
    ib = pow(beta, -1, N)
    G = [0] * N
    G[0] = _sum_r_times_mod(N, 0, (-beta) % N)
    half = N * (N - 1) // 2
    for cc in range(1, N):
        G[cc] = G[cc - 1] + half - N * ((ib * cc) % N)
    return sum(r1 * G[(-m - alpha * r1) % N] for r1 in range(N))


def fd_brute_C(n: int, a: Sequence[int], N: int) -> int:
    """O(N^{d-1}) direct evaluation of C (tests)."""
    d = len(a)
    il = pow(a[-1], -1, N)
    tot = 0
    for rs in product(range(N), repeat=d - 1):
        t = n + sum(x * y for x, y in zip(a[:-1], rs))
        rl = (-t * il) % N
        pr = rl
        for y in rs:
            pr *= y
        tot += pr
    return tot


def fd_exact(n: int, a: Sequence[int], N: int, brute: bool = False) -> Fraction:
    """s_n(a_1..a_d; N) as an exact rational."""
    d = len(a)
    C = fd_brute_C(n, a, N) if brute else fd_C(n, a, N)
    return Fraction((-1) ** d * (N * C - (N * (N - 1) // 2) ** d), N ** (d + 1))


def fd_scaled(n: int, a: Sequence[int], N: int) -> int:
    """N^d s_n(a; N), an integer for odd N (exact)."""
    d = len(a)
    C = fd_C(n, a, N)
    v = (-1) ** d * (N * C - (N * (N - 1) // 2) ** d)   # = N^{d+1} s
    assert v % N == 0
    return v // N


def fd_float(n: int, a: Sequence[int], N: int) -> complex:
    """Direct cyclotomic sum in floating point (tests only; never used for results)."""
    import cmath
    tot = 0j
    for k in range(1, N):
        z = cmath.exp(2j * cmath.pi * k * n / N)
        den = 1
        for x in a:
            den *= 1 - cmath.exp(2j * cmath.pi * k * x / N)
        tot += z / den
    return tot / N
