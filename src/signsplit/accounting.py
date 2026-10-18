"""Cost accounting for sign-splitting (MODEL code; uses N only as a size/number).

Exact list counts (count_only) for V0/V1/V3 come from hh.hh_search; V2 counts
from hh.v2_counts.  Parameter choice for every variant is a deterministic
function of (N, M, variant) -- frozen in notes/signsplit_prereg.md:

  m0 grid  : integers 1..4 plus round(2^{k/8}) for k >= 16, up to 2^{0.2 n + 4}/M,
             then integer refinement +-3 around the best grid point.
  V2       : (m0, tR, g) minimise the upper-bound formula
             G (2^{t-tR} g + 2^{tR} ceil((2 lambda + 1 + 2t)/g)),
             G = #square classes * #buckets, tR in 1..t, g in powers of 2.

Continuous per-prime merits (section 3 of the notes):  with giant ~ A m0 and
baby ~ B m0^{-3/2}, the optimum is ~ A^{3/5} B^{2/5}; a prime m changes A and B
by factors (x, y) and so the optimum by x^{3/5} y^{2/5}.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd, log2
from typing import Dict, List, Sequence, Tuple

from .hh import buckets, hh_search, lam, square_classes, v2_counts
from .identity import factor_squarefree


def m0_grid(N: int, M: int) -> List[int]:
    n = N.bit_length()
    hi = max(4, int(2 ** (0.2 * n + 4)) // max(1, M))
    out = {1, 2, 3, 4}
    k = 16
    while True:
        v = round(2 ** (k / 8))
        if v > hi:
            break
        out.add(v)
        k += 1
    return sorted(out)


def _count(N, M, m0, variant, anti=()):
    r = hh_search(N, M, m0, variant, anti=anti, count_only=True)["met"]
    return r["giant"] + r["baby"], r


def best_shared(N: int, M: int, variant: str = "V0", anti: Sequence[int] = ()) -> Dict:
    """Own-optimum m0 for V0/V1/V3 (count-only work = giant + baby)."""
    best = None
    for m0 in m0_grid(N, M):
        w, r = _count(N, M, m0, variant, anti)
        if best is None or w < best[0]:
            best = (w, m0)
    w0, m00 = best
    for m0 in range(max(1, m00 - 3), m00 + 4):
        w, r = _count(N, M, m0, variant, anti)
        if w < best[0] or (w == best[0] and m0 < best[1]):
            best = (w, m0)
    w, r = _count(N, M, best[1], variant, anti)
    return {"variant": variant, "M": M, "anti": list(anti), "m0": best[1], "work": w,
            "giant": r["giant"], "baby": r["baby"], "lambda": r["lambda"]}


def v2_formula(N: int, M: int, m0: int, tR: int, g: int) -> int:
    primes = factor_squarefree(M)
    t = sum(1 for m in primes if m % 2)
    G = len(square_classes(N, M)) * len(buckets(N, m0))
    K = 2 * lam(N, M, m0, M) + 1 + 2 * t
    return G * ((1 << (t - tR)) * g + (1 << tR) * (-(-K // g)))


def best_v2(N: int, M: int) -> Dict:
    primes = factor_squarefree(M)
    t = sum(1 for m in primes if m % 2)
    best = None
    for m0 in m0_grid(N, M):
        K = 2 * lam(N, M, m0, M) + 1 + 2 * t
        for tR in range(1, t + 1):
            g = 1
            while g <= 2 * K:
                f = v2_formula(N, M, m0, tR, g)
                if best is None or f < best[0]:
                    best = (f, m0, tR, g)
                g *= 2
    f, m0, tR, g = best
    return {"variant": "V2", "M": M, "m0": m0, "tR": tR, "g": g, "formula": f}


def v2_matched(N: int, M: int, m0: int) -> Dict:
    """V2 at the baseline's (M, m0): tR = floor(t/2) (>= 1), g = 1 (literal claim:
    MITM over the sign vector only)."""
    primes = factor_squarefree(M)
    t = sum(1 for m in primes if m % 2)
    tR = max(1, t // 2)
    return {"variant": "V2", "M": M, "m0": m0, "tR": tR, "g": 1,
            "formula": v2_formula(N, M, m0, tR, 1)}


# --------------------------------------------------------------------------
# continuous merits (exact rationals / floats only for display)
# --------------------------------------------------------------------------
def merit_lock(m: int) -> float:
    """HH lock on prime m: giant x (m-1), baby x m^{-3/2}  ->  ((m-1)/m)^{3/5}."""
    return ((m - 1) / m) ** 0.6


def merit_anti(m: int) -> float:
    """Sign collapse (anti-lock): giant x (m-1)/2, baby x m^{-1/2}."""
    return ((m - 1) / 2) ** 0.6 * m ** (-0.2)


def merit_pm_offsets(m: int) -> float:
    """Force V_m in {+-delta}: giant x (m-1)/2, baby x 2 (m^3/2)^{1/2}/m^2."""
    return ((m - 1) / 2) ** 0.6 * (2 ** 0.5 * m ** -0.5) ** 0.4


def merit_ideal_free(m: int) -> float:
    """Idealised (impossible) free sharing of the sign: HH lock with giant / 2 and
    baby x 2 -> extra factor 2^{-1/5} on top of the lock."""
    return merit_lock(m) * 2 ** -0.2


def merit_table(primes: Sequence[int] = (3, 5, 7, 11, 13, 17, 19, 23)) -> List[Dict]:
    return [{"m": m, "lock_HH": merit_lock(m), "anti": merit_anti(m),
             "pm_offsets": merit_pm_offsets(m), "ideal_free": merit_ideal_free(m),
             "anti_over_lock": merit_anti(m) / merit_lock(m),
             "anti_over_lock_formula": m ** 0.4 * 2 ** -0.6} for m in primes]


def v2_exponent_bound(M_primes: Sequence[int]) -> float:
    """Per-group MITM is per-target BSGS: cost ~ N^{1/4} phi(M) M^{-3/4} 2^{-t/2}
    at m0 = 1 (the constant returned); exponent 1/4 for any M."""
    c = 1.0
    for m in M_primes:
        if m == 2:
            c *= 1 * 2 ** -0.75
        else:
            c *= (m - 1) * m ** -0.75 * 2 ** -0.5
    return c
