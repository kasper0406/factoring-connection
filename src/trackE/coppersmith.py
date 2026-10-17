"""Univariate Coppersmith / Howgrave-Graham decoder for an *affine hint* on a factor.

OPERATIONAL: every function here receives only N and a public hint; nothing uses
p, q, phi(N) or a factorisation routine for N.  (Integer-root extraction of the
LLL output polynomial uses flint's ``fmpz_poly.roots()``, i.e. factorisation of an
integer polynomial over Z -- the standard final step of Coppersmith's method
[May10, §3 Step 4]; it never touches N.)

Affine hint format  H = (M, r, L, U)
------------------------------------
  "N has a divisor b with  b ≡ r (mod M)  and  L <= b <= U",   L >= 1 public.
Special cases:
  (a) low bits / residue: M = 2^k (or any M), r = b mod M, [L, U] = family promise
      (balanced promise used here:  L = floor(sqrt N)/2, U = floor(sqrt N); the
      smaller factor of a balanced semiprime from src/common/instances lies there).
  (b) interval / high bits: M = 1, r = 0, [L, U] = [A, A+W].
Write b = r + M (x_mid + y), |y| <= Y.  With c = (r M^{-1} + x_mid) mod N, the monic
f(y) = y + c satisfies f(y0) ≡ M^{-1} b ≡ 0 (mod b).

Lattice [May10, proof of Thm 3, δ = 1]:  dimension d = m + t,
  g_i(y) = N^{m-i} f(y)^i       (i = 0..m-1),     h_j(y) = y^j f(y)^m   (j = 0..t-1),
coefficient vectors of g(yY); lower triangular, det = N^{m(m+1)/2} Y^{d(d-1)/2}.
Every vector vanishes mod b^m at y0.

Rigorous success condition (used for parameter choice):
  Howgrave-Graham [May10, Thm 2]: if g(y0) ≡ 0 mod b^m, |y0| <= Y and
  ||g(yY)|| < b^m / sqrt(d), then g(y0) = 0 over Z.
  LLL (δ_LLL = 3/4) [May10, eq. before (1)]: ||v_1|| <= 2^{(d-1)/4} det^{1/d}.
  flint's LLL uses (δ, η) = (0.99, 0.51); its bound (δ-η^2)^{-(d-1)/4} = 1.370^{(d-1)/4}
  is smaller, so the 2^{(d-1)/4} condition is conservative.
  => sufficient:  2^{d(d-1)} * det^4 * d^{2d} < L^{4 m d}      (exact integer test).
Asymptotically (m = (d-1)/2, d odd, b ~ N^{1/2}): Y < N^{(d-1)/(4d)} up to the LLL and
sqrt(d) factors, i.e. Y -> N^{1/4} = N^{β²/δ} [May10, Thm 1].
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass
from math import comb, gcd, isqrt
from typing import List, Optional, Tuple

import flint

try:  # fplll's LLL wrapper (provably (0.99,0.51)-reduced output); ~10x faster than flint here
    from fpylll import IntegerMatrix as _IM, LLL as _FPLLL
except Exception:  # pragma: no cover
    _IM = None

from src.common.metrics import Metrics, NullMetrics
from src.common.ring import NonUnit, inv_mod


@dataclass(frozen=True)
class AffineHint:
    M: int
    r: int
    L: int
    U: int

    def x_range(self) -> Tuple[int, int]:
        """Integer x with L <= r + M x <= U."""
        lo = -((self.r - self.L) // self.M)          # ceil((L - r)/M)
        hi = (self.U - self.r) // self.M             # floor((U - r)/M)
        return lo, hi


def balanced_promise(N: int) -> Tuple[int, int]:
    """Public promise interval for the smaller factor of a balanced semiprime."""
    s = isqrt(N)
    return s // 2, s


# ---------------------------------------------------------------------------
# parameter theory (exact integer inequalities)
# ---------------------------------------------------------------------------

def rigorous_ok(N: int, L: int, Y: int, d: int, m: int) -> bool:
    """Exact check of 2^{(d-1)/4} det^{1/d} sqrt(d) < L^m  (all raised to power 4d)."""
    if not (1 <= m <= d):
        return False
    det4 = N ** (2 * m * (m + 1)) * Y ** (2 * d * (d - 1))   # det^4
    lhs = (1 << (d * (d - 1))) * det4 * d ** (2 * d)
    return lhs < L ** (4 * m * d)


def log2_maxY(n_bits_N: float, log2L: float, d: int, m: int, lll_factor: bool = True,
              hermite: float = 2 ** 0.25) -> float:
    """log2 of the largest Y allowed by the (rigorous or heuristic) condition.

    Condition: hermite^{d-1} det^{1/d} sqrt(d) < L^m with det = N^{m(m+1)/2} Y^{d(d-1)/2}.
    Rigorous: hermite = 2^{1/4}.  Heuristic average LLL: hermite ≈ 1.02 [Nguyen–Stehlé].
    Float output: for plots / threshold tables only.
    """
    if d == 1:
        return float("-inf")
    lh = math.log2(hermite) * (d - 1) if lll_factor else 0.0
    rhs = m * log2L - lh - 0.5 * math.log2(d) - m * (m + 1) / (2 * d) * n_bits_N
    return rhs * 2 / (d - 1)


def best_m(n_bits_N: float, log2L: float, d: int, **kw) -> Tuple[int, float]:
    best = (1, float("-inf"))
    for m in range(1, d + 1):
        v = log2_maxY(n_bits_N, log2L, d, m, **kw)
        if v > best[1]:
            best = (m, v)
    return best


def choose_params(N: int, L: int, Y: int, dmax: int = 200) -> Optional[Tuple[int, int]]:
    """Smallest dimension d (and its m) satisfying the rigorous condition, or None."""
    for d in range(2, dmax + 1):
        m0, _ = best_m(N.bit_length(), math.log2(L), d)
        for m in sorted(range(1, d + 1), key=lambda mm: abs(mm - m0))[:4]:
            if rigorous_ok(N, L, Y, d, m):
                return d, m
    return None


# ---------------------------------------------------------------------------
# lattice + LLL + roots
# ---------------------------------------------------------------------------

def build_basis(N: int, c: int, Y: int, d: int, m: int) -> List[List[int]]:
    """Rows = coefficient vectors of g_i(yY), h_j(yY) (low -> high), f(y) = y + c."""
    t = d - m
    # coefficients of f^k, k <= m
    fpows = [[comb(k, i) * c ** (k - i) for i in range(k + 1)] for k in range(m + 1)]
    Ypow = [1]
    for _ in range(d):
        Ypow.append(Ypow[-1] * Y)
    rows = []
    for i in range(m):
        Ni = N ** (m - i)
        row = [0] * d
        for j, a in enumerate(fpows[i]):
            row[j] = Ni * a * Ypow[j]
        rows.append(row)
    for j in range(t):
        row = [0] * d
        for i2, a in enumerate(fpows[m]):
            row[i2 + j] = a * Ypow[i2 + j]
        rows.append(row)
    return rows


def lll_reduce(B: List[List[int]], backend: str = "fpylll") -> List[List[int]]:
    if backend == "fpylll" and _IM is not None:
        A = _IM.from_matrix(B)
        _FPLLL.reduction(A)
        d1, d2 = A.nrows, A.ncols
        return [[int(A[i, j]) for j in range(d2)] for i in range(d1)]
    return [[int(a) for a in row] for row in flint.fmpz_mat(B).lll().tolist()]


def _int_roots(coeffs: List[int]) -> List[int]:
    poly = flint.fmpz_poly(coeffs)
    if poly.degree() < 1:
        return []
    return [int(r) for r, _ in poly.roots()]


def decode_affine(N: int, hint: AffineHint, d: Optional[int] = None, m: Optional[int] = None,
                  metrics: Optional[Metrics] = None, n_vectors: int = 3,
                  enforce_rigorous: bool = False, backend: str = "fpylll") -> Optional[int]:
    """Return a proper divisor of N consistent with the hint, or None.

    If d is None the smallest rigorous d is used (None returned if > dmax).
    Non-units met while inverting M mod N are returned as factors (rule 3).
    """
    mt = metrics if metrics is not None else NullMetrics()
    M, r = hint.M, hint.r
    mt.inc("gcd")
    gM = gcd(M, N)                      # rule 3: a nonunit modulus is itself output
    if 1 < gM < N:
        mt.set("status", "nonunit_M")
        return gM
    lo, hi = hint.x_range()
    if hi < lo:
        mt.set("status", "empty_range")
        return None
    x_mid = (lo + hi) // 2
    Y = max(1, hi - x_mid, x_mid - lo)
    mt.set("log2Y", round(math.log2(Y), 4))
    # tiny ranges: enumerate directly (cheaper than a lattice)
    if Y <= 2:
        for x in range(lo, hi + 1):
            b = r + M * x
            mt.inc("gcd")
            g = gcd(b, N)
            if 1 < g < N:
                mt.set("status", "enumerated")
                return g
        mt.set("status", "fail_enum")
        return None
    try:
        Minv = inv_mod(M, N, mt) if M % N != 1 else 1
    except NonUnit as e:
        if 1 < e.divisor < N:
            mt.set("status", "nonunit_M")
            return e.divisor
        raise
    c = (r * Minv + x_mid) % N
    if d is None:
        pr = choose_params(N, hint.L, Y)
        if pr is None:
            mt.set("status", "no_rigorous_params")
            return None
        d, m = pr
    elif m is None:
        m = best_m(N.bit_length(), math.log2(hint.L), d)[0]
    rig = rigorous_ok(N, hint.L, Y, d, m)
    mt.set("d", d); mt.set("m", m); mt.set("rigorous_condition", rig)
    if enforce_rigorous and not rig:
        mt.set("status", "below_rigorous")
        return None
    t0 = time.perf_counter()
    B = build_basis(N, c, Y, d, m)
    mt.hi("max_operand_bits", max(abs(v).bit_length() for row in B for v in row))
    mt.set("time_build", time.perf_counter() - t0)
    t0 = time.perf_counter()
    rows = lll_reduce(B, backend)
    mt.set("time_lll", time.perf_counter() - t0)
    mt.set("lll_backend", backend if (backend == "flint" or _IM is not None) else "flint")
    mt.inc("lll_calls")
    t0 = time.perf_counter()
    Ypow = [Y ** j for j in range(d)]
    rows.sort(key=lambda v: sum(int(a) * int(a) for a in v))
    try:
        for idx, row in enumerate(rows[:n_vectors]):
            coeffs = []
            for j, a in enumerate(row):
                a = int(a)
                assert a % Ypow[j] == 0
                coeffs.append(a // Ypow[j])
            mt.inc("root_calls")
            for y0 in _int_roots(coeffs):
                if abs(y0) > Y:
                    continue
                b = r + M * (x_mid + y0)
                mt.inc("gcd")
                g = gcd(b, N)
                if 1 < g < N:
                    mt.set("status", "ok")
                    mt.set("vector_index", idx)
                    return g
    finally:
        mt.set("time_roots", time.perf_counter() - t0)
    mt.set("status", "fail")
    return None


def decode_low_bits(N: int, k: int, r: int, d: Optional[int] = None, m: Optional[int] = None,
                    metrics: Optional[Metrics] = None, **kw) -> Optional[int]:
    """Format (a) with M = 2^k under the balanced promise."""
    L, U = balanced_promise(N)
    return decode_affine(N, AffineHint(1 << k, r % (1 << k), L, U), d=d, m=m, metrics=metrics, **kw)


def decode_residue(N: int, M: int, r: int, d: Optional[int] = None, m: Optional[int] = None,
                   metrics: Optional[Metrics] = None, **kw) -> Optional[int]:
    """Format (a) for an arbitrary modulus M under the balanced promise."""
    L, U = balanced_promise(N)
    return decode_affine(N, AffineHint(M, r % M, L, U), d=d, m=m, metrics=metrics, **kw)


def decode_interval(N: int, A: int, W: int, d: Optional[int] = None, m: Optional[int] = None,
                    metrics: Optional[Metrics] = None, **kw) -> Optional[int]:
    """Format (b): a divisor in [A, A+W]; L = max(A,1) is the HG lower bound."""
    return decode_affine(N, AffineHint(1, 0, max(A, 1), A + W), d=d, m=m, metrics=metrics, **kw)
