"""Calibration gate: the abstract collision-search cost model reproduces
BSGS-on-(p+q) (1/4), Lehman (1/3), Hittmeir 2021 (2/9), Harvey 2021 (1/5) and
the Harvey--Hittmeir 2022 log-log gain.

MODEL CODE (no factors are used anywhere; N is only a size parameter here).

Abstract model (notes/cycloharvey_results.md, section 2)
-------------------------------------------------------
A *target* tau is (handle h_tau, window [U_tau, U_tau + W_tau)) with the
guarantee h_tau == alpha^{u*} (mod p) for the hidden integer u* = g_tau(p) in
the window, for the branch tau that is "right" for p.  Strategies:

  * square-test (Lehman):     cost = sum_tau W_tau               (no BSGS)
  * per-target BSGS:          cost = sum_tau 2*ceil(sqrt(W_tau))
  * Hittmeir two-tier:        per-target BSGS for k<=xi, one giant step each
                              for xi<k<=eta against a shared baby list
  * Harvey full sharing:      cost = m + sum_tau ceil(W_tau / m)
  * Harvey--Hittmeir:         one giant step per branch (sigma mod M, sigma0
                              bucket), shared baby list of 2*lambda+1 steps
                              in units of M^2.

All window widths are computed with exact integer arithmetic (isqrt);
floating point appears only in the exponent fits.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd, isqrt
from typing import Dict, List, Tuple

import numpy as np


# --------------------------------------------------------------------------
# exact helpers
# --------------------------------------------------------------------------
def ceil_sqrt_ratio(num: int, den: int) -> int:
    """Smallest integer w >= 0 with den*w^2 >= num  (= ceil(sqrt(num/den)))."""
    if num <= 0:
        return 0
    w = isqrt(num // den)
    while den * w * w < num:
        w += 1
    while w > 0 and den * (w - 1) * (w - 1) >= num:
        w -= 1
    return w


def iroot_floor(N: int, k: int) -> int:
    """floor(N^(1/k)) exactly."""
    if N < 2:
        return N
    x = int(round(N ** (1.0 / k)))
    while x ** k > N:
        x -= 1
    while (x + 1) ** k <= N:
        x += 1
    return x


def divisor_counts(r: int) -> np.ndarray:
    """d[k] = number of ordered pairs (a,b) with ab = k, for 0<=k<=r."""
    d = np.zeros(r + 1, dtype=np.int64)
    for a in range(1, r + 1):
        d[a::a] += 1
    return d


def prefix(arr: np.ndarray) -> np.ndarray:
    """P[k] = sum_{j<=k} arr[j] (exact int64)."""
    return np.cumsum(arr, dtype=np.int64)


def sum_ceil_sqrt(N: int, C: int, r: int, P: np.ndarray, floor1: bool = False) -> int:
    """Exact  sum_{k=1}^{r} w(k) * ceil(sqrt(N/(C k)))  where P is the prefix
    sum of the weights w (P[k] = sum_{j<=k} w(j)); with floor1, each ceil is
    replaced by max(1, ceil).  Iterates over k or over the distinct values of
    the ceiling, whichever is fewer (exact integer arithmetic throughout)."""
    if r <= 0:
        return 0
    wmax = ceil_sqrt_ratio(N, C)          # value at k = 1
    wmin = ceil_sqrt_ratio(N, C * r)      # value at k = r
    if r <= (wmax - wmin) + 2:
        tot = 0
        for k in range(1, r + 1):
            v = ceil_sqrt_ratio(N, C * k)
            if floor1:
                v = max(1, v)
            tot += int(P[k] - P[k - 1]) * v
        return tot
    # group: value(k) <= w  iff  C k w^2 >= N  iff k >= ceil(N / (C w^2))
    tot = 0
    hi = r
    w = wmin
    while hi >= 1:
        # k in [lo, hi] all have ceil value == w (value nonincreasing in k)
        if w == 0:
            lo = 1
        else:
            # smallest k with value <= w: k >= N/(C w^2)
            lo = max(1, -(-N // (C * w * w)))
        lo = min(lo, hi)
        v = max(1, w) if floor1 else w
        tot += int(P[hi] - P[lo - 1]) * v
        hi = lo - 1
        if hi >= 1:
            w = ceil_sqrt_ratio(N, C * hi)
    return tot


def lehman_window(N: int, r: int, k: int) -> int:
    """ceil(N^{1/2} / (4 r k^{1/2}))  -- Harvey Lemma 3.3 / Lehman bound."""
    return ceil_sqrt_ratio(N, 16 * r * r * k)


# --------------------------------------------------------------------------
# the five calibrated algorithms (cost = number of list elements / candidates)
# --------------------------------------------------------------------------
def cost_bsgs_sum(N: int, c_num: int = 2, c_den: int = 1) -> Dict:
    """BSGS on s = p+q for q/p <= c: s in [ceil(2 sqrt N), ceil(2 sqrt N)+W].
    W = max_{p in [sqrt(N/c), sqrt(N)]} (p + N/p) - 2 sqrt N, attained at the
    lower end; computed exactly as an integer upper bound."""
    plo = ceil_sqrt_ratio(N * c_den, c_num)          # ceil(sqrt(N/c))
    s0 = ceil_sqrt_ratio(4 * N, 1)                   # ceil(2 sqrt N)
    smax = plo + -(-N // plo)                        # p + ceil(N/p)
    W = max(1, smax - s0 + 1)
    m = ceil_sqrt_ratio(W, 1)
    return {"alg": "bsgs_sum", "W": W, "baby": m, "giant": -(-W // m),
            "cost": m + -(-W // m)}


def cost_lehman(N: int, r: int = None) -> Dict:
    """Lehman: square-test every candidate x with 0 <= x - sqrt(4kN) < W_k,
    k <= r (per k, not per pair; Harvey Remark 3.4).  Default r = N^{1/3}."""
    if r is None:
        r = max(1, iroot_floor(N, 3))
    cand = sum_ceil_sqrt(N, 16 * r * r, r, np.arange(r + 1, dtype=np.int64))
    return {"alg": "lehman", "r": r, "candidates": cand, "k_values": r,
            "cost": cand + r}


def _pairs_by_k(d: np.ndarray, lo: int, hi: int) -> int:
    return int(d[lo + 1:hi + 1].sum()) if hi > lo else 0


def cost_hittmeir(N: int, xi: int, eta: int, d: np.ndarray = None) -> Dict:
    """Hittmeir 2021 Algorithm 5.1 counts (list elements):
    tier 1: for every (a<=b, ab<=xi): BSGS over Lambda_ab = ceil(N^{1/2}/((ab)^{1/2} eta))
            with 2*ceil(sqrt(Lambda)) elements;
    tier 2: baby list of size ceil(N^{1/2}/(xi^{1/2} eta)) plus one giant step
            for each (a<=b, xi<ab<=eta)."""
    if d is None or len(d) <= eta:
        d = divisor_counts(eta)
    U = prefix((d + 1) // 2)                   # unordered pairs a<=b, ab<=k
    tier1 = 0
    for k in range(1, xi + 1):
        lam = ceil_sqrt_ratio(N, k * eta * eta)
        tier1 += int(U[k] - U[k - 1]) * 2 * ceil_sqrt_ratio(lam, 1)
    baby = ceil_sqrt_ratio(N, xi * eta * eta)
    giant = int(U[eta] - U[xi])
    return {"alg": "hittmeir", "xi": xi, "eta": eta, "tier1": tier1,
            "baby": baby, "giant": giant, "cost": tier1 + baby + giant}


def cost_harvey(N: int, r: int, m: int, d: np.ndarray = None) -> Dict:
    """Harvey 2021 Algorithm 4.2 counts: m baby steps, and for every ordered
    pair ab <= r, ceil(W_ab / m) giant steps with W_ab = N^{1/2}/(4 r (ab)^{1/2})
    (Step 2b: 0 <= j < W/m).  Also returns the pair count (handles)."""
    if d is None or len(d) <= r:
        d = divisor_counts(r)
    P = prefix(d)
    # number of j with 0 <= j < N^{1/2}/(4 r m k^{1/2}) is ceil(.) (>= 1)
    giant = sum_ceil_sqrt(N, 16 * r * r * m * m, r, P, floor1=True)
    pairs = int(P[r])
    return {"alg": "harvey", "r": r, "m": m, "baby": m, "giant": giant,
            "pairs": pairs, "cost": m + giant}


def primorial_upto(limit: int) -> List[int]:
    out, M, q = [1], 1, 2
    while True:
        while any(q % s == 0 for s in range(2, isqrt(q) + 1)):
            q += 1
        if M * q > limit:
            return out
        M *= q
        out.append(M)
        q += 1


def phi(n: int) -> int:
    res, x, f = n, n, 2
    while f * f <= x:
        if x % f == 0:
            while x % f == 0:
                x //= f
            res -= res // f
        f += 1
    if x > 1:
        res -= res // x
    return res


def cost_hh(N: int, M: int, m0: int, c_num: int = 2, c_den: int = 1) -> Dict:
    """Harvey--Hittmeir 2022 Algorithm 4.3 counts, restricted to the balanced
    range p in [sqrt(N/c), sqrt(N)]:
      baby  = 2*lambda + 1,  lambda = ceil(4 N^{1/2} / (M m0)^{3/2})
      giant = phi(M) * (#sigma0 buckets with ratio 1 + 1/m0 covering the range).
    """
    lam = ceil_sqrt_ratio(16 * N, (M * m0) ** 3)
    plo = ceil_sqrt_ratio(N * c_den, c_num)
    phi_hi = isqrt(N)
    s0, buckets = plo, 0
    while s0 <= phi_hi:
        buckets += 1
        s0 = -(-(s0 * (m0 + 1)) // m0)
    giant = phi(M) * buckets
    return {"alg": "hh", "M": M, "m0": m0, "baby": 2 * lam + 1,
            "buckets": buckets, "giant": giant, "cost": 2 * lam + 1 + giant}


# --------------------------------------------------------------------------
# optimisation over parameters (geometric grids)
# --------------------------------------------------------------------------
def _geo(lo: int, hi: int, ratio: float = 2 ** 0.25) -> List[int]:
    out, x = [], float(max(1, lo))
    while x <= hi:
        v = int(round(x))
        if not out or v != out[-1]:
            out.append(v)
        x *= ratio
    if not out or out[-1] != hi:
        out.append(hi)
    return out


def best_harvey(N: int, rmax_cap: int = 1 << 22) -> Dict:
    n = N.bit_length()
    r_hi = min(rmax_cap, max(4, int(2 ** (0.30 * n))))
    d = divisor_counts(r_hi)
    P = prefix(d)
    best = None
    for r in _geo(1, r_hi):
        sw = sum_ceil_sqrt(N, 16 * r * r, r, P)          # total window
        m0 = max(1, isqrt(sw))
        for m in _geo(max(1, m0 // 8), 8 * m0 + 1):
            c = cost_harvey(N, r, m, d)
            if best is None or c["cost"] < best["cost"]:
                best = c
    return best


def best_hittmeir(N: int, cap: int = 1 << 22) -> Dict:
    n = N.bit_length()
    e_hi = min(cap, max(4, int(2 ** (0.32 * n))))
    d = divisor_counts(e_hi)
    best = None
    for eta in _geo(1, e_hi):
        for xi in _geo(1, eta):
            c = cost_hittmeir(N, xi, eta, d)
            if best is None or c["cost"] < best["cost"]:
                best = c
    return best


def best_lehman(N: int) -> Dict:
    n = N.bit_length()
    best = None
    for r in _geo(1, max(2, int(2 ** (0.40 * n)))):
        if r > (1 << 22):
            break
        c = cost_lehman(N, r)
        if best is None or c["cost"] < best["cost"]:
            best = c
    return best


def best_hh(N: int, Ms: List[int] = None) -> Dict:
    n = N.bit_length()
    Ms = Ms or primorial_upto(1 << 16)
    best = None
    for M in Ms:
        for m0 in _geo(1, max(2, int(2 ** (0.2 * n + 3)) // M)):
            c = cost_hh(N, M, m0)
            if best is None or c["cost"] < best["cost"]:
                best = c
    return best


# --------------------------------------------------------------------------
# exact continuous-exponent versions (rational grid; exact at 1/360)
# --------------------------------------------------------------------------
def exponent_lp(name: str, den: int = 360) -> Tuple[Fraction, Dict]:
    """Minimise the max of the leading exponents of each algorithm's cost
    terms (log_N of each term, logs dropped) over a rational grid."""
    F = Fraction
    grid = [F(i, den) for i in range(0, den // 2 + 1)]
    best, arg = None, None
    if name == "bsgs_sum":
        return F(1, 4), {}
    if name == "lehman":           # sum_k N^{1/2}/(r k^{1/2}) ~ N^{1/2} r^{-1/2}; + r
        for rho in grid:
            v = max(F(1, 2) - rho / 2, rho)
            if best is None or v < best:
                best, arg = v, {"rho": rho}
    elif name == "hittmeir":       # N^{1/4} xi^{3/4} eta^{-1/2}, N^{1/2} xi^{-1/2} eta^{-1}, eta
        for eta in grid:
            for xi in grid:
                if xi > eta:
                    break
                v = max(F(1, 4) + F(3, 4) * xi - eta / 2,
                        F(1, 2) - xi / 2 - eta, eta)
                if best is None or v < best:
                    best, arg = v, {"xi": xi, "eta": eta}
    elif name == "harvey":         # m, r, N^{1/2} r^{-1/2} m^{-1}
        for rho in grid:
            for mu in grid:
                v = max(mu, rho, F(1, 2) - rho / 2 - mu)
                if best is None or v < best:
                    best, arg = v, {"rho": rho, "mu": mu}
    elif name == "hh":             # phi(M) m0 ~ M m0 = T, N^{1/2} T^{-3/2}
        for tau in grid:
            v = max(tau, F(1, 2) - F(3, 2) * tau)
            if best is None or v < best:
                best, arg = v, {"T": tau}
    else:
        raise ValueError(name)
    return best, arg
