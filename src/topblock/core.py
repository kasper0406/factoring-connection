"""topblock core: exact identities linking the ladder's window queries to divisors / near-divisors.

Notation (Theorem L_j, notes/analytic_results.md §3).  Block start x, offset h >= 0, degree j,
D = x^{j+1}, S_j(h) = sum_{i<=j} (-h)^i x^{j-i}, num_j(h) = N*S_j(h), P_j(h) = num_j(h)/D,
R(h) = N(-h)^{j+1}/(x^{j+1}(x+h)) so that N/(x+h) = P_j(h) + R(h).

Identities (PROVED in notes/topblock_results.md §1; checked here exactly):
  (T1)  (x+h) S_j(h) = x^{j+1} - (-h)^{j+1}.
  (T2)  num_j(h) - q D = S_j(h) (N - q(x+h)) - q(-h)^{j+1}     for every integer q.
  (W)   {P_j(h)} in [a, b]  <=>  exists q in Z with  k := N - q(x+h)  satisfying
        a(x+h) + N(-h)^{j+1}/D  <=  k  <=  b(x+h) + N(-h)^{j+1}/D        (0 <= a <= b < 1).
  (F)   with w_j(h, r) = r - num_j(h) and r_j(h) = -N(-h)^{j+1}/(x+h):
        (x+h) w_j(h, r_j(h)) = -N D,   and   F_j := (x+h) r + N(-h)^{j+1} = (x+h) w_j + N D.

Everything here is HARNESS/DIAGNOSTIC except `locate_divisor_via_window_oracle`, whose only access to
N's factors is through a supplied window oracle (an oracle-relative reduction).  The brute-force
oracles iterate over the block; they test identities and count calls, never costs.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd, isqrt
from typing import Callable, Dict, List, Optional, Tuple


# ---------------------------------------------------------------- polynomials (exact integers)

def S_poly(x: int, j: int, h: int) -> int:
    return sum((-h) ** i * x ** (j - i) for i in range(j + 1))


def num_poly(N: int, x: int, j: int, h: int) -> int:
    return N * S_poly(x, j, h)


def frac_P(N: int, x: int, j: int, h: int) -> Fraction:
    """{P_j(h)} as an exact fraction in [0, 1)."""
    D = x ** (j + 1)
    return Fraction(num_poly(N, x, j, h) % D, D)


def check_T1_T2(x: int, j: int, h: int, N: int, q: int) -> bool:
    D = x ** (j + 1)
    S = S_poly(x, j, h)
    t1 = (x + h) * S == x ** (j + 1) - (-h) ** (j + 1)
    t2 = num_poly(N, x, j, h) - q * D == S * (N - q * (x + h)) - q * (-h) ** (j + 1)
    return t1 and t2


# ---------------------------------------------------------------- (W): window = near-divisor set

def window_set_direct(N: int, x: int, j: int, H: int, a: Fraction, b: Fraction) -> List[int]:
    """h in [0, H) with {P_j(h)} in [a, b]  (brute force)."""
    return [h for h in range(H) if a <= frac_P(N, x, j, h) <= b]


def window_set_neardiv(N: int, x: int, j: int, H: int, a: Fraction, b: Fraction) -> List[int]:
    """h in [0, H) for which some q gives k = N - q(x+h) in the (W) range (brute force over q)."""
    D = x ** (j + 1)
    out = []
    for h in range(H):
        d = x + h
        off = Fraction(N * (-h) ** (j + 1), D)
        lo = a * d + off
        hi = b * d + off
        # k = N - q d  in [lo, hi]  <=>  q in [(N - hi)/d, (N - lo)/d]
        qlo = _ceil_frac(Fraction(N) - hi, d)
        qhi = _floor_frac(Fraction(N) - lo, d)
        if qlo <= qhi:
            out.append(h)
    return out


def _floor_frac(v: Fraction, d: int) -> int:
    return (v.numerator) // (v.denominator * d)


def _ceil_frac(v: Fraction, d: int) -> int:
    return -((-v.numerator) // (v.denominator * d))


def near_divisor_pairs(N: int, x: int, H: int, K: int) -> List[Tuple[int, int]]:
    """All (h, k) with 0 <= h < H, |k| <= K, (x+h) | N - k  (brute force)."""
    out = []
    for h in range(H):
        d = x + h
        r = N % d
        for k in (r, r - d):
            if abs(k) <= K:
                out.append((h, k))
    return sorted(set(out))


# ---------------------------------------------------------------- reduction (iii): oracle => KB

class WindowOracle:
    """HARNESS: brute-force top-block window oracle for one block, counting calls.

    query(h1, h2, a, b) = #{h in [h1, h2) : {P_j(h)} in [a, b]}  (closed window, exact rationals).
    In the ladder this is two floor sums (Lemma W of floorhard plus the eta trick of Theorem L_j)."""

    def __init__(self, N: int, x: int, j: int, H: int):
        self.N, self.x, self.j, self.H = N, x, j, H
        self.D = x ** (j + 1)
        self.vals = [Fraction(num_poly(N, x, j, h) % self.D, self.D) for h in range(H)]
        self.calls = 0

    def query(self, h1: int, h2: int, a: Fraction, b: Fraction) -> int:
        self.calls += 1
        return sum(1 for h in range(h1, h2) if a <= self.vals[h] <= b)


def divisor_window(N: int, x: int, j: int, h1: int, h2: int) -> Tuple[Fraction, Fraction]:
    """Closed window containing {P_j(h)} for every exact divisor x+h, h in [h1, h2) (h1 >= 1).

    For a divisor, {P_j(h)} = -R(h) mod 1 with |R| increasing in h: even j gives |R(h)|, odd j gives
    1 - R(h).  Window = hull over [h1, h2-1]."""
    def Rabs(h):
        return Fraction(N * h ** (j + 1), x ** (j + 1) * (x + h))
    lo, hi = Rabs(h1), Rabs(h2 - 1)
    if j % 2 == 0:
        return lo, hi
    return 1 - hi, 1 - lo


def locate_divisor_via_window_oracle(N: int, x: int, j: int, H: int,
                                     oracle_query: Callable[[int, int, Fraction, Fraction], int]
                                     ) -> Tuple[List[int], int]:
    """ORACLE-RELATIVE (operational given the oracle): all divisors of N in [x, x+H).

    Bisection with windows shrunk to the exact-divisor window of each subinterval (reduction (iii)).
    Returns (divisors, number of oracle calls made by this routine)."""
    calls = 0
    found: List[int] = []
    if x > 1 and N % x == 0:
        found.append(x)
    stack = [(1, H)] if H > 1 else []
    while stack:
        h1, h2 = stack.pop()
        a, b = divisor_window(N, x, j, h1, h2)
        calls += 1
        if oracle_query(h1, h2, a, b) == 0:
            continue
        if h2 - h1 == 1:
            d = x + h1
            g = gcd(d, N)
            if 1 < g < N and N % d == 0:
                found.append(d)
            continue
        mid = (h1 + h2) // 2
        stack.append((mid, h2))
        stack.append((h1, mid))
    return sorted(found), calls


def weighted_neardiv_bound(N: int, x: int, j: int, H: int) -> float:
    """C^w of notes §1.3: 2(log2 H + 1) + 2 sum_{0<|k|<=K'} nu(k) (log2(K'/|k|) + 1), K' = 2(j+1)rho(x+H).

    HARNESS (brute force); float for reporting only."""
    import math
    rho = Fraction(N * H ** (j + 1), x ** (j + 2))
    Kp = 2 * (j + 1) * rho * (x + H)
    Ki = int(Kp)
    tot = 2 * (math.log2(max(H, 2)) + 1)
    for (h, k) in near_divisor_pairs(N, x, H, Ki):
        if k != 0:
            tot += 2 * (math.log2(float(Kp) / abs(k)) + 1)
    return tot


# ---------------------------------------------------------------- reduction (iv): divisors-in-interval => oracle

def divisors_in_interval_bruteforce(M: int, lo: int, hi: int) -> List[int]:
    """HARNESS stand-in for a general 'all divisors of M in [lo, hi)' solver."""
    return [d for d in range(max(lo, 1), hi) if M % d == 0]


def window_count_via_divisor_finder(N: int, x: int, j: int, h1: int, h2: int, a: Fraction, b: Fraction,
                                    finder=divisors_in_interval_bruteforce) -> Tuple[int, int]:
    """Answer #{h in [h1,h2): {P_j(h)} in [a,b]} from a divisors-in-interval finder applied to N - k.

    By (W), h counts iff some k in [a d + off(h), b d + off(h)] has d | N - k, where d = x+h and
    off(h) = N(-h)^{j+1}/D.  All such k lie in [kmin, kmax] below; each k needs one finder call, and
    each returned d is checked against its own exact range.  Returns (count, finder calls)."""
    D = x ** (j + 1)
    offs = [Fraction(N * (-h) ** (j + 1), D) for h in (h1, h2 - 1)]
    dmax = x + h2 - 1
    kmin = _floorF(min(a * x + min(offs), a * dmax + min(offs)))
    kmax = -_floorF(-(max(b * dmax + max(offs), b * x + max(offs))))
    hit = set()
    calls = 0
    for k in range(kmin, kmax + 1):
        M = N - k
        if M <= 0:
            continue
        calls += 1
        for d in finder(M, x + h1, x + h2):
            h = d - x
            off = Fraction(N * (-h) ** (j + 1), D)
            if a * d + off <= k <= b * d + off:
                hit.add(h)
    return len(hit), calls


def _floorF(v: Fraction) -> int:
    return v.numerator // v.denominator


# ---------------------------------------------------------------- block geometry

def ladder_B(N: int, j: int) -> int:
    """B = ceil(N^{1/(j+2)})."""
    B = int(round(N ** (1.0 / (j + 2))))
    while B ** (j + 2) < N:
        B += 1
    while B > 1 and (B - 1) ** (j + 2) >= N:
        B -= 1
    return B


def top_block_containing(N: int, p: int, j: int, c: Fraction = Fraction(1, 2)) -> Tuple[int, int]:
    """HARNESS: a ladder-style block [x, x+H), H = floor(c x / B), with x <= p < x + H."""
    B = ladder_B(N, j)
    H = int(c * p // B)
    x = p - H // 2
    H = int(c * x // B)
    assert x <= p < x + H
    return x, H
