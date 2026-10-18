"""Divisor-cover Lehman (DCL): a several-moduli Lehman/Harvey search.

MODEL + OPERATIONAL-PARAMETER CODE.  Uses only N and public parameters
(never p, q, phi or group orders).  Primes M here are small public moduli
enumerated by a sieve; no factorisation of N is attempted.

Mechanism (notes/multimod_results.md, section 2)
------------------------------------------------
For a weight pair v = (A, B) (gcd 1) and a prime M not dividing 2ABN, let r be
a root of B r^2 == A N (mod M).  If the hidden p satisfies p == r (mod M) then
M | Bp - Aq and

      F := B p + A q  ==  c := B r + A N r^{-1}   (mod M^2)

(the Harvey--Hittmeir Prop. 3.3 lock).  So g = (F - c)/M^2 is an integer in a
window [g_lo, g_hi] fixed by the balanced range p in [sqrt(N/2), sqrt(N)].
One target (v, M, r) covers EVERY integer p == r (mod M) of the range when the
window is the whole range (we always use whole-range windows), so coverage is
a pure sieve: p is covered iff some chosen (v, M) has M | B p^2 - A N.

No residue is guessed: the classes r are determined by N.  The several moduli
share ONE baby list gamma^i, gamma = beta^L, L = prod M^2 (handles are raised
to L/M^2), see operational.py.
"""
from __future__ import annotations

import heapq
from fractions import Fraction
from math import gcd, isqrt, log
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# ----------------------------------------------------------------------------
# small public helpers
# ----------------------------------------------------------------------------
def primes_upto(n: int) -> List[int]:
    if n < 2:
        return []
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def sqrt_mod_prime(a: int, M: int) -> List[int]:
    """All x in [0, M) with x^2 == a (mod M), M an odd PUBLIC prime (Tonelli-Shanks)."""
    a %= M
    if a == 0:
        return [0]
    if pow(a, (M - 1) // 2, M) != 1:
        return []
    if M % 4 == 3:
        x = pow(a, (M + 1) // 4, M)
    else:
        Q, S = M - 1, 0
        while Q % 2 == 0:
            Q //= 2
            S += 1
        z = 2
        while pow(z, (M - 1) // 2, M) != M - 1:
            z += 1
        c, x, t, R = pow(z, Q, M), pow(a, (Q + 1) // 2, M), pow(a, Q, M), S
        while t != 1:
            i, tt = 0, t
            while tt != 1:
                tt = tt * tt % M
                i += 1
            b = pow(c, 1 << (R - i - 1), M)
            x, c, t, R = x * b % M, b * b % M, t * b * b % M, i
    return sorted({x, (M - x) % M})


def isqrt_ceil(n: int) -> int:
    s = isqrt(n)
    return s if s * s == n else s + 1


# ----------------------------------------------------------------------------
# windows
# ----------------------------------------------------------------------------
def f_range(N: int, A: int, B: int, c_num: int = 2, c_den: int = 1) -> Tuple[int, int]:
    """Integer bounds F_lo <= B p + A N/p <= F_hi for real p in [sqrt(N/c), sqrt(N)].
    Exact: F is convex in p; minimum at the stationary point sqrt(AN/B) if it lies
    in the range, else at an endpoint; maximum at an endpoint."""
    # endpoint values squared: p = sqrt(N/c): F^2 = N (B + cA)^2 / c ; p = sqrt N: F^2 = N (A+B)^2
    e1_num, e1_den = N * (B * c_den + c_num * A) ** 2, c_num * c_den  # F^2 = N (B + cA)^2 / c with c = c_num/c_den
    # (B/sqrt(c) + A sqrt(c))^2 N = N (B + cA)^2 / c ; with c = cn/cd: N (B cd + A cn)^2 / (cn cd)
    e2 = N * (A + B) ** 2
    hi = max(isqrt(e1_num // e1_den) + 1, isqrt(e2) + 1)
    # stationary point p* = sqrt(AN/B) in [sqrt(N/c), sqrt N]  <=>  A/B in [1/c, 1]
    if Fraction(A, B) <= 1 and Fraction(A, B) >= Fraction(c_den, c_num):
        lo = isqrt_ceil(4 * A * B * N)
    else:
        lo = min(isqrt(e1_num // e1_den), isqrt(e2))  # floor of the smaller endpoint value (safe)
        lo = max(0, lo - 1)
    return lo, hi


def target_window(N: int, A: int, B: int, M: int, r: int, F_lo: int, F_hi: int) -> Tuple[int, int, int]:
    """(c, g_lo, W): class constant c (mod M^2) and the whole-range window of g."""
    M2 = M * M
    s = N * pow(r, -1, M2) % M2
    c = (B * r + A * s) % M2
    g_lo = -((-(F_lo - c)) // M2)          # ceil
    g_hi = (F_hi - c) // M2               # floor
    return c, g_lo, max(1, g_hi - g_lo + 1)


# ----------------------------------------------------------------------------
# weight pool and greedy design (from N and public parameters only)
# ----------------------------------------------------------------------------
def weight_pool(Hmax: int) -> List[Tuple[int, int]]:
    out = []
    for A in range(1, Hmax + 1):
        for B in range(1, Hmax + 1):
            if gcd(A, B) == 1:
                out.append((A, B))
    return out


def _giants(span: int, M: int, m: int) -> int:
    """Predicted ceil(W/m) for a whole-range window (W <= span/M^2 + 2)."""
    return -(-(span // (M * M) + 2) // m)


def design(N: int, lam: float, m: int, Hmax: int = 12, c_num: int = 2, c_den: int = 1,
           _primes: Optional[List[int]] = None) -> Dict:
    """Greedy choice of (v, M) pairs in decreasing efficiency rho/(M * ceil(W/m))
    (rho = number of roots, known from N) until the expected cover count
    lambda = sum rho/M reaches `lam`.  Efficiency per v is unimodal in M with
    its peak near M_v = sqrt(span_v/m), so candidates are generated lazily by
    a heap over per-v two-pointer walks.  Uses only N and public parameters.
    Returns the target list and the model cost m + sum ceil(W/m)."""
    from bisect import bisect_left
    pool = weight_pool(Hmax)
    spans = {}
    for (A, B) in pool:
        F_lo, F_hi = f_range(N, A, B, c_num, c_den)
        spans[(A, B)] = (F_lo, F_hi, F_hi - F_lo)
    Mpk = max(isqrt(sp // max(1, m)) for (_, _, sp) in spans.values()) + 3
    primes = _primes if _primes is not None else primes_upto(max(1000, 64 * Mpk))
    primes = [M for M in primes if M > 2 and N % M]
    heap = []
    state = {}
    for v, (F_lo, F_hi, sp) in spans.items():
        i0 = bisect_left(primes, isqrt(sp // max(1, m)) + 1)
        # walk down from i0-1 and up from i0
        state[v] = [i0 - 1, i0]
        for idx in (i0 - 1, i0):
            if 0 <= idx < len(primes):
                M = primes[idx]
                heapq.heappush(heap, (M * _giants(sp, M, m), idx, v, 0 if idx < i0 else 1))
    targets = []
    tot_lam, giants = 0.0, 0
    while heap and tot_lam < lam:
        _, idx, v, up = heapq.heappop(heap)
        A, B = v
        F_lo, F_hi, sp = spans[v]
        M = primes[idx]
        nxt = idx + 1 if up else idx - 1
        if 0 <= nxt < len(primes):
            M2 = primes[nxt]
            heapq.heappush(heap, (M2 * _giants(sp, M2, m), nxt, v, up))
        if (A * B) % M == 0:
            continue
        a = A * N * pow(B, -1, M) % M
        roots = [r for r in sqrt_mod_prime(a, M) if r % M]
        if not roots:
            continue
        for r in roots:
            c, g_lo, W = target_window(N, A, B, M, r, F_lo, F_hi)
            k = -(-W // m)
            targets.append({"A": A, "B": B, "M": M, "r": r, "c": c, "g_lo": g_lo, "W": W, "k": k})
            giants += k
        tot_lam += len(roots) / M
    moduli = sorted({t["M"] for t in targets})
    return {"lam": tot_lam, "m": m, "targets": targets, "giants": giants,
            "cost": m + giants, "n_targets": len(targets), "n_moduli": len(moduli),
            "Mmin": moduli[0] if moduli else None, "Mmax": moduli[-1] if moduli else None,
            "Hmax_used": max((max(t["A"], t["B"]) for t in targets), default=0)}


def best_design(N: int, lam: float, Hmax: int = 12, ratio: float = 2 ** 0.25) -> Dict:
    """Optimise the baby length m on a geometric grid (cost = m + sum ceil(W/m))."""
    n = N.bit_length()
    best = None
    m = max(4, int(2 ** (n / 6 - 3)))
    top = int(2 ** (n / 6 + 4))
    # one prime table, large enough for the smallest m
    sp_max = max(f_range(N, A, B)[1] - f_range(N, A, B)[0] for (A, B) in weight_pool(Hmax))
    primes = primes_upto(max(1000, 64 * (isqrt(sp_max // m) + 3)))
    while m <= top:
        d = design(N, lam, m, Hmax, _primes=primes)
        if d["lam"] >= lam and (best is None or d["cost"] < best["cost"]):
            best = d
        m = int(m * ratio) + 1
    return best


# ----------------------------------------------------------------------------
# exact coverage of ALL integer candidates p in the range (verification only;
# cost ~ sum_tau P/M, not part of the algorithm)
# ----------------------------------------------------------------------------
def candidate_range(N: int, c_num: int = 2, c_den: int = 1) -> Tuple[int, int]:
    """Integers p with N/c <= p^2 <= N."""
    lo = isqrt_ceil(-(-N * c_den // c_num))
    while lo * lo * c_num < N * c_den:
        lo += 1
    return lo, isqrt(N)


def uncovered_candidates(N: int, targets: Sequence[Dict], c_num: int = 2, c_den: int = 1,
                         lo: Optional[int] = None, hi: Optional[int] = None):
    """numpy bool array over [lo, hi]: True where no target covers p.
    A target (A,B,M,r,c,g_lo,W) covers p iff p == r (mod M) and the real value
    (B p + A N/p - c)/M^2 lies in [g_lo, g_lo + W).  With whole-range windows
    the second condition always holds; it is nevertheless checked exactly for
    every marked p in `verify_windows`."""
    import numpy as np
    if lo is None:
        lo, hi = candidate_range(N, c_num, c_den)
    unc = np.ones(hi - lo + 1, dtype=bool)
    for t in targets:
        M, r = t["M"], t["r"]
        start = (r - lo) % M
        unc[start::M] = False
    return lo, hi, unc


def verify_windows(N: int, targets: Sequence[Dict], c_num: int = 2, c_den: int = 1,
                   samples_per_target: int = 3) -> int:
    """Exact check (rational arithmetic) that the extreme candidates of every
    target class lie inside its window: g(p) = (B p + A N/p - c)/M^2 is convex,
    so it suffices to check the two endpoints of the class within the range and
    the class members nearest the stationary point.  Returns #violations.
    (For a non-divisor candidate p the value g(p) is not an integer; the
    criterion is that ceil(g) and floor(g) both lie in the window.)"""
    lo, hi = candidate_range(N, c_num, c_den)
    bad = 0
    for t in targets:
        A, B, M, r, c, g_lo, W = t["A"], t["B"], t["M"], t["r"], t["c"], t["g_lo"], t["W"]
        first = lo + ((r - lo) % M)
        last = hi - ((hi - r) % M)
        if first > hi:
            continue
        pts = {first, last}
        ps = isqrt(A * N // B)
        for d in range(-samples_per_target, samples_per_target + 1):
            x = ps + d * M
            x = x - ((x - r) % M)
            for y in (x, x + M):
                if first <= y <= last:
                    pts.add(y)
        for p in pts:
            # if p were a divisor, g(p) would be an integer; the window must
            # contain every integer the real value could round to that way:
            # ceil(g) >= g_lo and floor(g) <= g_lo + W - 1.
            g = Fraction(B * p * p + A * N - c * p, M * M * p)
            gc = -((-g.numerator) // g.denominator)
            gf = g.numerator // g.denominator
            if gc < g_lo or gf > g_lo + W - 1:
                bad += 1
    return bad
