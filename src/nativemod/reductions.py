"""nativemod OPERATIONAL reductions: factor N = pq from one "fixed-level" datum.

OPERATIONAL: every function receives only N and the stated public datum.  Nothing here uses p, q, phi(N), a
factorisation routine for N, a primality test on N, or PARI.  (Integer-root extraction of an integer polynomial
uses flint's fmpz_poly.roots(), the standard last step of Coppersmith's method; it never touches N.)

Contents
- eis_factor:        Prop EIS-k.  a(N) of E_k(chi, chi') (k >= 2), real characters given as callables.
- crt_trace_factor:  Prop CRT-T.  t = CRT(a_p mod p, a_q mod q) for a fixed elliptic curve (|a_r| <= 2 sqrt r).
- class_pair_factor: Prop CLS.  The unordered pairs {p mod l, q mod l} for primes l <= L (the determinants of the
                     two Frobenius classes of a weight-2 mod-l representation), then CRT on s = p + q.
The lattice code is copied (rule 7) from src/trackE/coppersmith.py (May 2010, delta = 1 lattice) and adapted to
the monic linear polynomial y + c with an unknown divisor b >= L.
"""
from __future__ import annotations

import math
from math import comb, gcd, isqrt
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import flint
import gmpy2


# ---------------------------------------------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------------------------------------------

def _iroot(x: int, k: int) -> Optional[int]:
    """Exact integer k-th root of x >= 0, or None."""
    if x < 0:
        return None
    r, exact = gmpy2.iroot(gmpy2.mpz(x), k)
    return int(r) if exact else None


def _split_from_s(N: int, s: int) -> Optional[int]:
    """Given s = p + q (candidate), return p if s^2 - 4N is a square with the right parity."""
    D = s * s - 4 * N
    if D < 0:
        return None
    r = isqrt(D)
    if r * r != D or (s - r) % 2:
        return None
    p = (s - r) // 2
    if 1 < p < N and N % p == 0:
        return p
    return None


def _trial(N: int, B: int = 1000) -> Optional[int]:
    for d in range(2, B):
        if N % d == 0 and d < N:
            return d
    return None


# ---------------------------------------------------------------------------------------------------------------
# Prop EIS-k
# ---------------------------------------------------------------------------------------------------------------

def eis_factor(N: int, a: int, k: int, chi: Callable[[int], int], chip: Callable[[int], int],
               cond: int = 4) -> Optional[int]:
    """Factor N = pq from a(N) of E_k(chi, chi'), a(n) = sum_{d|n} chi(n/d) chi'(d) d^{k-1}, k >= 2.

    chi, chi' are real (+-1-valued on units) characters of conductor dividing `cond`, given as callables that need
    only n mod cond (e.g. Kronecker symbols).  For N = pq coprime to cond:
        a(N) = chi(N) + chi'(N) N^{k-1} + u p^{k-1} + v q^{k-1},   u = chi(q)chi'(p), v = chi(p)chi'(q) in {+-1}.
    For each (u, v) solve v X^2 - R X + u N^{k-1} = 0 in X = q^{k-1} (R = a - chi(N) - chi'(N) N^{k-1}), take the
    exact (k-1)-th root and test divisibility.  Deterministic, at most 4 patterns."""
    g = gcd(N, cond)
    if 1 < g < N:
        return g
    if k < 2:
        return None
    R = a - chi(N) - chip(N) * N ** (k - 1)
    Nk = N ** (k - 1)
    for u in (1, -1):
        for v in (1, -1):
            # v X^2 - R X + u Nk = 0
            disc = R * R - 4 * v * u * Nk
            if disc < 0:
                continue
            r = isqrt(disc)
            if r * r != disc:
                continue
            for num in (R + r, R - r):
                if num % (2 * v):
                    continue
                X = num // (2 * v)
                if X <= 1:
                    continue
                qq = _iroot(X, k - 1)
                if qq is not None and 1 < qq < N and N % qq == 0:
                    return qq
    return None


# ---------------------------------------------------------------------------------------------------------------
# Prop CRT-T: Coppersmith / Howgrave-Graham for x - t = 0 mod b, b >= L, |x| small
# ---------------------------------------------------------------------------------------------------------------

def _build_basis(N: int, c: int, Y: int, d: int, m: int) -> List[List[int]]:
    """Rows: g_i(yY) = N^{m-i} f(yY)^i (i < m), h_j = (yY)^j f(yY)^m (j < d - m), f(y) = y + c (May 2010)."""
    t = d - m
    fpows = [[comb(kk, i) * c ** (kk - i) for i in range(kk + 1)] for kk in range(m + 1)]
    Ypow = [Y ** j for j in range(d + 1)]
    rows = []
    for i in range(m):
        row = [0] * d
        for j, a in enumerate(fpows[i]):
            row[j] = N ** (m - i) * a * Ypow[j]
        rows.append(row)
    for j in range(t):
        row = [0] * d
        for i2, a in enumerate(fpows[m]):
            row[i2 + j] = a * Ypow[i2 + j]
        rows.append(row)
    return rows


def _rigorous_ok(N: int, L: int, Y: int, d: int, m: int) -> bool:
    det4 = N ** (2 * m * (m + 1)) * Y ** (2 * d * (d - 1))
    return (1 << (d * (d - 1))) * det4 * d ** (2 * d) < L ** (4 * m * d)


def small_root_divisor(N: int, c: int, Y: int, L: int, d: int, m: int, stats: Optional[Dict] = None) -> Optional[int]:
    """Find y0, |y0| <= Y, with b | (y0 + c) for some divisor b >= L of N; return a proper divisor gcd(y0 + c, N)."""
    rows = [[int(a) for a in r] for r in flint.fmpz_mat(_build_basis(N, c, Y, d, m)).lll().tolist()]
    if stats is not None:
        stats["lll_calls"] = stats.get("lll_calls", 0) + 1
    rows.sort(key=lambda v: sum(a * a for a in v))
    for row in rows[:3]:
        coeffs = []
        for j, a in enumerate(row):
            assert a % (Y ** j) == 0
            coeffs.append(a // Y ** j)
        poly = flint.fmpz_poly(coeffs)
        if poly.degree() < 1:
            continue
        for r0, _ in poly.roots():
            y0 = int(r0)
            if abs(y0) > Y:
                continue
            g = gcd(y0 + c, N)
            if 1 < g < N:
                return g
    return None


def crt_trace_factor(N: int, t: int, ratio_bound: Tuple[int, int] = (4, 3), d: int = 9, m: int = 4,
                     stats: Optional[Dict] = None) -> Optional[int]:
    """Factor N = pq (q/p < ratio_bound) from t with t = a_q (mod q), |a_q| <= 2 sqrt(q), q the larger prime.

    Root x0 = a_q of x - t = 0 (mod q), q >= sqrt(N) =: L.  The window |x| <= 2 sqrt(q_max) is covered by
    sub-windows of radius Y, each one Coppersmith call with (d, m) (shifted polynomial y + (x_mid - t)).  For
    the asymptotic statement take d = Theta(log N); then O(1) windows suffice (May 2010, Thm 6 with c_N = O(1)).
    The tiny-window branch enumerates gcd(t - x, N) directly."""
    if stats is None:
        stats = {}
    g = _trial(N, 64)
    if g:
        return g
    L = isqrt(N)
    if L * L < N:
        L += 1
    a, b = ratio_bound
    qmax = isqrt(N * a // b) + 1
    X = 2 * isqrt(qmax) + 2
    # largest Y with the exact rigorous condition for (d, m)
    lo, hi = 1, X
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if _rigorous_ok(N, L, mid, d, m):
            lo = mid
        else:
            hi = mid - 1
    Y = lo
    stats["Y_bits"] = Y.bit_length()
    if Y <= 4:
        for x in range(-X, X + 1):
            g = gcd(t - x, N)
            if 1 < g < N:
                stats["windows"] = "enum"
                return g
        return None
    nwin = 0
    x_mid = -X + Y
    while x_mid - Y <= X:
        nwin += 1
        # f(y) = y + c with y = x - x_mid ... we need (x - t) = 0 mod q, x = x_mid + y  =>  y + (x_mid - t)
        g = small_root_divisor(N, (x_mid - t) % N, Y, L, d, m, stats)
        if g:
            stats["windows"] = nwin
            return g
        x_mid += 2 * Y + 1
    stats["windows"] = nwin
    return None


# ---------------------------------------------------------------------------------------------------------------
# Prop CLS: unordered pairs {p mod l, q mod l} -> s mod l -> CRT -> factor
# ---------------------------------------------------------------------------------------------------------------

def class_pair_factor(N: int, pairs: Dict[int, Tuple[int, int]]) -> Optional[int]:
    """pairs[l] = unordered {p mod l, q mod l} for distinct primes l; needs prod(l) > s_max (s <= N + 1)."""
    M, s = 1, 0
    for l, (u, v) in pairs.items():
        r = (u + v) % l
        # CRT
        inv = pow(M, -1, l) if M % l else None
        if inv is None:
            continue
        s = s + M * (((r - s) * inv) % l)
        M *= l
    if M <= N + 1:
        return None
    return _split_from_s(N, s % M)
