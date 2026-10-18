"""Cost-model calculations for notes/coverage38_results.md (model level).

Cost model (multimod / below16 M+): a design is a finite set of targets; target tau has a
modulus M, a weight v of height h, a coverage weight pi (probability over the population
S_n that tau covers the true pair) and a window of w g-values.  The cost is
    X >= max(T, sqrt(sum w))      (T = number of (weight, modulus) target groups).
Window hypothesis (WIN): w_tau >= c1 * pi_tau^2 * h_tau * P  (P = sqrt N).
  * full windows:     pi ~ 1/M,              w ~ h P / M^2
  * Lehman windows:   pi ~ Delta/(M h P),    w ~ Delta^2 / (h P M^2)
Both satisfy (WIN) with equality up to constants.

Theorem U:  X >= (c P)^{2/5} * lam^{4/5} * J^{-1/5},   lam = sum pi,  J = #distinct moduli.
All floating point here is model arithmetic (exponents, constants), not factoring output.
"""
from __future__ import annotations

from math import sqrt
from typing import Dict, Iterable, List, Sequence, Tuple


# ---------------------------------------------------------------------------
# exponents
# ---------------------------------------------------------------------------
def exponent_from_theta(theta: float) -> float:
    """Worst-case exponent of N if coverage is provable at lam = J^theta (single band,
    J ~ Y): (1+theta)/(6+theta), capped at the Minkowski value 1/5."""
    return min(0.2, (1.0 + theta) / (6.0 + theta))


def theta_from_c(c: float) -> float:
    """Theorem W's c (H = sqrt(Y) J^{-c})  <->  lam = J^{1-2c}."""
    return 1.0 - 2.0 * c


def gain_factor(lam: float, J: float) -> float:
    """Cost relative to N^{1/5}: (lam^4 / J)^{1/5}  (Theorem U, attained by one band)."""
    return (lam ** 4 / J) ** 0.2


# ---------------------------------------------------------------------------
# designs
# ---------------------------------------------------------------------------
def design_stats(targets: Iterable[Tuple[int, int, float, float]], P: float) -> Dict:
    """targets: (M, h, pi, w).  Returns T, J, lam, sum_w, cost = max(T, sqrt(sum w))."""
    T = 0
    mods = set()
    lam = 0.0
    sw = 0.0
    for M, h, pi, w in targets:
        T += 1
        mods.add(M)
        lam += pi
        sw += w
    return {"T": T, "J": len(mods), "lam": lam, "sum_w": sw, "cost": max(T, sqrt(sw))}


def universal_lower_bound(lam: float, J: int, P: float, c1: float = 1.0, R: int = 1) -> float:
    """Theorem U with explicit constant: X >= (c2 P)^{2/5} lam^{4/5} J^{-1/5},
    c2 = c1 / (2 sqrt 3 R)  (R = max number of root-targets per (weight, modulus))."""
    c2 = c1 / (2.0 * sqrt(3.0) * R)
    return (c2 * P) ** 0.4 * lam ** 0.8 * J ** (-0.2)


def heights_sorted(t: int) -> List[int]:
    """Heights of the t lowest-height signed primitive weights (A != 0, B >= 1)."""
    from math import gcd
    out = []
    h = 1
    while len(out) < t:
        # weights of height exactly h
        cur = []
        for B in range(1, h + 1):
            for A in range(-h, h + 1):
                if A == 0 or max(abs(A), B) != h or gcd(abs(A), B) != 1:
                    continue
                cur.append(h)
        out.extend(cur)
        h += 1
    return out[:t]


def box_design(P: float, Y: float, J: int, H: int) -> List[Tuple[int, int, float, float]]:
    """J moduli 'near Y' (labels Y+1..Y+J; model only), all weights of height <= H,
    full windows: pi = 1/M, w = h P / M^2."""
    hs = heights_sorted_upto(H)
    out = []
    for k in range(J):
        M = int(Y) + 1 + k
        for h in hs:
            out.append((M, h, 1.0 / M, h * P / M ** 2))
    return out


def lehman_design(P: float, Y: float, J: int, H: int, Delta: float):
    hs = heights_sorted_upto(H)
    out = []
    for k in range(J):
        M = int(Y) + 1 + k
        for h in hs:
            d = min(Delta, h * P)
            out.append((M, h, d / (M * h * P), d * d / (h * P * M ** 2)))
    return out


def heights_sorted_upto(H: int) -> List[int]:
    from math import gcd
    return [max(abs(A), B) for B in range(1, H + 1) for A in range(-H, H + 1)
            if A != 0 and gcd(abs(A), B) == 1]


# ---------------------------------------------------------------------------
# closed-form optima of the two single-band families (continuous model)
# ---------------------------------------------------------------------------
def box_cost(P: float, Y: float, J: float, lam: float) -> float:
    """Full windows: T = lam Y, sum w = lam H P / Y, H = sqrt(lam Y / J)."""
    H = sqrt(lam * Y / J)
    return max(lam * Y, sqrt(lam * H * P / Y))


def lehman_cost(P: float, J: float, H: float, lam: float) -> float:
    """Lehman windows (Y free in [J, J H^2/lam]): T = J H^2, sum w = lam^2 P /(J H)."""
    return max(J * H * H, lam * sqrt(P / (J * H)))


def optimise(fn, lo: float, hi: float, iters: int = 200) -> Tuple[float, float]:
    """Golden-section on log scale (unimodal model functions)."""
    import math
    a, b = math.log(lo), math.log(hi)
    g = (sqrt(5) - 1) / 2
    for _ in range(iters):
        c = b - g * (b - a)
        d = a + g * (b - a)
        if fn(math.exp(c)) < fn(math.exp(d)):
            b = d
        else:
            a = c
    x = math.exp((a + b) / 2)
    return x, fn(x)
