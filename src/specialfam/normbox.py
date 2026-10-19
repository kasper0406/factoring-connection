"""Algorithm A (norm boxes) and the generic polynomial-box search (Theorem G).

OPERATIONAL. Inputs: N and public family parameters (eta, box shapes, polynomial
Phi). Never uses the factors. Primality testing is not used. Every nonunit met is
output (rule 3): gcds are checked, and a leading-coefficient nonunit is returned.

Generic search (Theorem G, P1 packaging). Phi(x1, w) is an integer polynomial whose
degree in x1 is d and whose x1-leading coefficient lc(w) is an integer polynomial in
the giant data w. For a box x1 in [x0, x0+A), w in W, write x1 = x0 + m v + j with
0 <= j < m. For each giant point (v, w) the block product
    prod_{j<m} Phi(x0 + m v + j, w) = lc(w)^m * prod_{j<m} chi_{v,w}(j),
    chi_{v,w}(Y) := lc(w)^{-1} Phi(x0 + m v + Y, w)   (monic in Y),
is obtained from src.specialfam.polyprod.block_products with the universal baby
polynomial F(Y) = prod_{j<m}(Y - j). A gcd equal to N is resolved by scanning that
block's m leaves directly (valid when every leaf satisfies 0 < |Phi| < N).

Sum of two squares (Algorithm A). Phi(a, b) = a^2 + b^2, split a, giant b. This is
literally r18's (Z/N)[i] mechanism: (Y + m v)^2 + b^2 is the characteristic
polynomial of -(m v + i b), and Nm_{(Z/N)[i]}(F(m v + i b)) = prod_j ((m v + j)^2 + b^2).
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd, isqrt
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import gmpy2

from src.common.metrics import Metrics, NullMetrics
from src.common.ring import NonUnit
from src.specialfam.polyprod import block_products, monic_from

SMALL_PRIMES = (2, 3, 5, 7)


def small_trial(N: int, mt: Metrics) -> Optional[int]:
    for ell in SMALL_PRIMES:
        mt.inc("gcd")
        if N % ell == 0 and N != ell:
            return ell
    return None


def polybox_search(N: int,
                   shift_coeffs: Callable[[int, object], Sequence[int]],
                   value: Callable[[int, object], int],
                   x0: int, A: int, giants: Sequence[object], m: Optional[int] = None,
                   metrics: Optional[Metrics] = None) -> Dict:
    """Search primes p | N dividing Phi(x1, w), x1 in [x0, x0+A), w in giants.

    shift_coeffs(t, w): integer coefficients (low->high in Y) of Phi(t + Y, w).
    value(x1, w): the integer Phi(x1, w) (used only in the gcd = N scan).
    Returns {"factor": d or None, "block": (v, w) or None, "units": m + d*K, ...}.
    """
    mt = metrics if metrics is not None else NullMetrics()
    W = list(giants)
    if A <= 0 or not W:
        return {"factor": None, "units": 0}
    if m is None:
        m = max(1, isqrt(A * len(W)))
    m = max(1, min(m, A))
    V = -(-A // m)
    S = list(range(m))
    d = len(shift_coeffs(x0, W[0])) - 1
    mt.hi("baby_m", m)
    # batch giants so that each remainder tree has root degree about m
    batch = max(1, m // max(d, 1))
    pts = [(v, w) for w in W for v in range(V)]
    K = len(pts)
    units = m + d * K
    mt.inc("giants", K)
    mt.inc("babies", m)
    for i0 in range(0, K, batch):
        chunk = pts[i0:i0 + batch]
        chis = []
        for (v, w) in chunk:
            try:
                chi, _lc = monic_from(N, list(shift_coeffs(x0 + m * v, w)))
            except NonUnit as e:
                if 1 < e.divisor < N:
                    return {"factor": e.divisor, "block": (v, w), "units": units,
                            "how": "leading_coeff_nonunit", "m": m, "K": K}
                raise
            chis.append(chi)
        vals = block_products(N, S, chis, mt)
        for (v, w), val in zip(chunk, vals):
            mt.inc("gcd")
            g = gcd(val, N)
            if g == 1:
                continue
            if g < N:
                return {"factor": g, "block": (v, w), "units": units, "how": "block_gcd",
                        "m": m, "K": K}
            # gcd = N: scan the m leaves of this block (cost O(m), once)
            mt.inc("full_gcd_blocks")
            for j in S:
                x1 = x0 + m * v + j
                leaf = value(x1, w)
                mt.inc("gcd")
                h = gcd(leaf, N)
                if 1 < h < N:
                    mt.inc("scan_leaves", j + 1)
                    return {"factor": h, "block": (v, w), "units": units + m, "how": "scan",
                            "m": m, "K": K}
            mt.inc("full_gcd_unresolved")
    return {"factor": None, "units": units, "m": m, "K": K}


# ---------------------------------------------------------------------------
# Algorithm A: p = a^2 + b^2, 1 <= b <= p^eta
# ---------------------------------------------------------------------------

def _sq_shift(t: int, b: int) -> List[int]:
    return [t * t + b * b, 2 * t, 1]


def _sq_value(a: int, b: int) -> int:
    return a * a + b * b


def band_box(k: int, eta: Fraction) -> Tuple[int, int]:
    """Public box for band 2^{k-1} < p <= 2^k: a in [0, A], b in [1, B]."""
    A = isqrt(1 << k)
    B = int(gmpy2.iroot(gmpy2.mpz(1) << (k * eta.numerator), eta.denominator)[0])
    return A, max(B, 1)


def normbox_factor(N: int, eta: Fraction, kmin: int = 4, kmax: Optional[int] = None,
                   metrics: Optional[Metrics] = None) -> Dict:
    """Algorithm A. Public input: N and eta. Scans dyadic bands k = kmin, kmin+1, ...

    Guarantee (Theorem A): if N has a prime factor p = a^2 + b^2 with 1 <= b <= p^eta
    and no prime factor below 11, a proper factor is returned by band ceil(log2 p).
    """
    mt = metrics if metrics is not None else NullMetrics()
    t = small_trial(N, mt)
    if t:
        return {"factor": t, "how": "trial", "band": None, "units": 0}
    units = 0
    k = kmin
    while True:
        if kmax is not None and k > kmax:
            break
        A, B = band_box(k, eta)
        if (2 * A + 1) ** 2 + B * B >= N:  # leaves must stay below N (scan validity)
            break
        res = polybox_search(N, _sq_shift, _sq_value, 0, A + 1, range(1, B + 1),
                             m=isqrt((A + 1) * B) + 1, metrics=mt)
        units += res.get("units", 0)
        if res.get("factor"):
            res.update({"band": k, "units_total": units, "A": A, "B": B})
            return res
        k += 1
    return {"factor": None, "units_total": units, "band": None}


# ---------------------------------------------------------------------------
# Other packaged families (Theorem G instances), used by tests / E3
# ---------------------------------------------------------------------------

def disc_box_factor(N: int, a0: int, b0: int, L: int, metrics: Optional[Metrics] = None) -> Dict:
    """p = a^2 + b^2 with (a, b) in [a0, a0+L) x [b0, b0+L) (HKL-type polar box)."""
    return polybox_search(N, _sq_shift, _sq_value, a0, L, range(b0, b0 + L),
                          m=L, metrics=metrics)


def near_square_factor(N: int, Amax: int, Kmax: int, metrics: Optional[Metrics] = None) -> Dict:
    """p = a^2 + k with 0 <= a <= Amax, 1 <= k <= Kmax (not a norm form)."""
    return polybox_search(N, lambda t, k: [t * t + k, 2 * t, 1], lambda a, k: a * a + k,
                          0, Amax + 1, range(1, Kmax + 1), metrics=metrics)


def cubic_form_factor(N: int, X: int, Y: int, metrics: Optional[Metrics] = None) -> Dict:
    """p = x^3 + 2 y^3 with 1 <= x <= X, 1 <= y <= Y (Heath-Brown's form)."""
    return polybox_search(N, lambda t, y: [t ** 3 + 2 * y ** 3, 3 * t * t, 3 * t, 1],
                          lambda x, y: x ** 3 + 2 * y ** 3, 1, X, range(1, Y + 1), metrics=metrics)


def ternary_factor(N: int, X: int, Y: int, Z: int, metrics: Optional[Metrics] = None) -> Dict:
    """p = x^2 + y^2 + z^2, long side x (P1 packaging needs no norm structure)."""
    giants = [(y, z) for y in range(0, Y + 1) for z in range(1, Z + 1)]
    return polybox_search(N, lambda t, w: [t * t + w[0] ** 2 + w[1] ** 2, 2 * t, 1],
                          lambda x, w: x * x + w[0] ** 2 + w[1] ** 2, 0, X + 1, giants,
                          metrics=metrics)
