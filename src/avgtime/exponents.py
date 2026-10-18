"""Exact (Fraction) exponent bookkeeping for notes/avgtime_results.md.

All exponents are in units of n = log2 N, so "x" means 2^{x n} = N^x, up to
(log N)^{O(1)} factors.  P = N^{1/2}.  Free parameters of one DCL-det stage:

  h   H = N^h            weight height
  y   Y = N^y            band scale
  mu  m = N^mu           baby-list length (m <= mbar is all the success proof needs)
  t   T = N^t            order threshold (analysis only for DCL; detector length if used)

Bounds proved in notes/dcltheorem_results.md (Theorem D', section 4), with Lemma 10
evaluated at the true baby length m (the delta-range of Lemma 7 is (-m, m)):

  cost   = max(mu, 2h + y, 3h + 1/2 - y - mu)        [m + T_k (1 + mbar/m)]
  E_a    = 2t - 1/2                                  [n T^2 / P]
  E_b    = -2 j h                                    [(n/H^2)^j; proved j = 1]
  E_c1   = 3h + mu - t                               [n H^3 m / T]
  E_c2   = 2h + y + mu - 1/2                         [n H^2 Y m / P]
  C(T)   = 3t - 1   (detector failure: ord_p 2 = ord_q 2 <= T)  [T^3 / P^2]
  detector cost = t/2

Side constraints: Y > H^2 (y > 2h, Lemma 9), Y^4 > sqrt2 P (y > 1/8, Lemma 6),
large-sieve level 256 Y^2 <= P (y < 1/4), 0 <= t < 1/2.

Nothing here is an experiment: it is exact linear algebra on these exponents.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from typing import Dict, List, Optional, Sequence, Tuple

ONE_FIFTH = F(1, 5)      # Harvey--Hittmeir fallback exponent
ONE_SIXTH = F(1, 6)


# ---------------------------------------------------------------------------
# 1. Audit of Theorem D' as written: H = n 2^{eps n/2}, Y = 2^{n/6} H^{-1/3},
#    T = 2^{(1/4 - eps) n}, m = mbar.  Each requirement is a + b*eps < 0.
# ---------------------------------------------------------------------------
def dprime_exponents(eps: F) -> Dict[str, F]:
    h = eps / 2
    y = ONE_SIXTH - h / 3
    t = F(1, 4) - eps
    mu = h + F(1, 2) - 2 * y          # mbar = H P / Y^2
    return {"h": h, "y": y, "t": t, "mu": mu,
            "cost": max(mu, 2 * h + y, 3 * h + F(1, 2) - y - mu),
            "E_a": 2 * t - F(1, 2), "E_b": -2 * h,
            "E_c1": 3 * h + mu - t, "E_c2": 2 * h + y + mu - F(1, 2)}


def dprime_requirements() -> Dict[str, Tuple[F, F]]:
    """name -> (a, b): requirement a + b*eps < 0 (strict: polylog factors).

    Separations: Y > H^2 (Lemma 9 distinct classes), Y^4 > sqrt2 P (Lemma 6).
    Exceptional terms must be o(N^{-eps}) * polylog^{-1}: E_x + eps < 0.
    E_b is N^{-eps} exactly by the choice H = n 2^{eps n/2} (no requirement).
    'useful' = D' cheaper than HH (not needed for validity)."""
    def lin(f):
        a = f(F(0))
        b = f(F(1)) - a
        return a, b
    req = {
        "Y>H^2 (Lemma 9)": lin(lambda e: 2 * dprime_exponents(e)["h"] - dprime_exponents(e)["y"]),
        "Y^4>sqrt2 P (Lemma 6)": lin(lambda e: F(1, 2) - 4 * dprime_exponents(e)["y"]),
        "E_a < N^-eps": lin(lambda e: dprime_exponents(e)["E_a"] + e),
        "E_c1 < N^-eps": lin(lambda e: dprime_exponents(e)["E_c1"] + e),
        "E_c2 < N^-eps": lin(lambda e: dprime_exponents(e)["E_c2"] + e),
        "useful: cost < 1/5": lin(lambda e: dprime_exponents(e)["cost"] - ONE_FIFTH),
    }
    return req


def sup_eps(a: F, b: F) -> Optional[F]:
    """sup of eps > 0 with a + b eps < 0 (None = no upper limit; requires a < 0 or b < 0)."""
    if b <= 0:
        return None if a < 0 or b < 0 else F(0)
    return -a / b


def dprime_admissible() -> Tuple[F, str]:
    """Binding validity constraint of D' as written (excluding 'useful')."""
    best, name = None, None
    for k, (a, b) in dprime_requirements().items():
        if k.startswith("useful"):
            continue
        s = sup_eps(a, b)
        if s is not None and (best is None or s < best):
            best, name = s, k
    return best, name


# ---------------------------------------------------------------------------
# 2. A small exact LP solver (two-phase tableau simplex, Bland's rule).
#    maximise c.x  s.t.  A x <= b,  x >= 0.   Fractions throughout.
# ---------------------------------------------------------------------------
def lp_max(c: Sequence[F], A: Sequence[Sequence[F]], b: Sequence[F]
           ) -> Tuple[Optional[F], Optional[List[F]]]:
    m, nv = len(A), len(c)
    # columns: x (nv) | slack/surplus (m) | artificial (m)
    rows, basis, art = [], [], []
    for i in range(m):
        Ai, bi = [F(v) for v in A[i]], F(b[i])
        s = [F(0)] * m
        a = [F(0)] * m
        if bi >= 0:
            s[i] = F(1)
            rows.append(Ai + s + a + [bi])
            basis.append(nv + i)
        else:
            s[i] = F(-1)
            a[i] = F(1)
            rows.append([-v for v in Ai] + s + a + [-bi])
            basis.append(nv + m + i)
            art.append(nv + m + i)
    ncol = nv + 2 * m

    def pivot(r, col):
        pv = rows[r][col]
        rows[r] = [v / pv for v in rows[r]]
        for i in range(m):
            if i != r and rows[i][col] != 0:
                f = rows[i][col]
                rows[i] = [vi - f * vr for vi, vr in zip(rows[i], rows[r])]
        basis[r] = col

    def simplex(obj, allowed):
        # obj: maximise obj.x over columns; reduced costs computed each step
        while True:
            cb = [obj[basis[i]] for i in range(m)]
            enter = None
            for j in allowed:
                if j in basis:
                    continue
                red = obj[j] - sum(cb[i] * rows[i][j] for i in range(m))
                if red > 0:
                    enter = j
                    break
            if enter is None:
                return "opt"
            best, r = None, None
            for i in range(m):
                if rows[i][enter] > 0:
                    ratio = rows[i][-1] / rows[i][enter]
                    if best is None or ratio < best or (ratio == best and basis[i] < basis[r]):
                        best, r = ratio, i
            if r is None:
                return "unbounded"
            pivot(r, enter)

    allcols = list(range(ncol))
    if art:
        obj1 = [F(0)] * ncol
        for j in art:
            obj1[j] = F(-1)
        simplex(obj1, allcols)
        if sum(rows[i][-1] for i in range(m) if basis[i] in art) != 0:
            return None, None                       # infeasible
        # drive remaining artificials out of the basis
        for i in range(m):
            if basis[i] in art:
                for j in range(nv + m):
                    if rows[i][j] != 0:
                        pivot(i, j)
                        break
    obj2 = [F(v) for v in c] + [F(0)] * (2 * m)
    st = simplex(obj2, list(range(nv + m)))
    if st == "unbounded":
        return None, None
    x = [F(0)] * nv
    for i in range(m):
        if basis[i] < nv:
            x[basis[i]] = rows[i][-1]
    return sum(F(ci) * xi for ci, xi in zip(c, x)), x


class LP:
    """Tiny builder: named nonnegative variables, constraints sum coef*var <= rhs."""

    def __init__(self):
        self.names: List[str] = []
        self.cons: List[Tuple[Dict[str, F], F]] = []

    def var(self, name: str) -> str:
        self.names.append(name)
        return name

    def le(self, coefs: Dict[str, F], rhs) -> None:
        self.cons.append(({k: F(v) for k, v in coefs.items()}, F(rhs)))

    def minimise(self, name: str) -> Tuple[F, Dict[str, F]]:
        idx = {v: i for i, v in enumerate(self.names)}
        c = [F(0)] * len(self.names)
        c[idx[name]] = F(-1)
        A, b = [], []
        for coefs, rhs in self.cons:
            row = [F(0)] * len(self.names)
            for k, v in coefs.items():
                row[idx[k]] += v
            A.append(row)
            b.append(rhs)
        val, x = lp_max(c, A, b)
        if val is None:
            raise ValueError("LP infeasible or unbounded")
        return -val, dict(zip(self.names, x))


# ---------------------------------------------------------------------------
# 3. One DCL stage as LP constraints: cost <= C, failure fraction <= N^{-e}.
# ---------------------------------------------------------------------------
def _stage(lp: LP, s: str, C: Dict[str, F], e: str, j: F, detector: bool,
           t_shared: Optional[str] = None, dyadic: bool = False) -> None:
    """Add stage s with variables h,y,mu,(t) ; constraints cost <= C (a linear form
    given as dict, e.g. {'z':1,'e0':1}) and e <= each failure exponent."""
    h, y, mu = lp.var(f"h{s}"), lp.var(f"y{s}"), lp.var(f"mu{s}")
    t = t_shared or lp.var(f"t{s}")
    neg = {k: -v for k, v in C.items()}
    lp.le({mu: 1, **neg}, 0)
    lp.le({h: 2, y: 1, **neg}, 0)
    lp.le({h: 3, y: -1, mu: -1, **neg}, F(-1, 2))
    lp.le({mu: 1, h: -1, y: 2}, F(1, 2))     # m <= mbar  (mu <= h + 1/2 - 2y; avgtime_1 item 8)
    # side constraints
    lp.le({y: -1}, F(-1, 8))                 # Y^4 > sqrt2 P   (infimum)
    lp.le({h: 2, y: -1}, 0)                  # Y > H^2         (infimum)
    lp.le({y: 1}, F(1, 4))                   # large-sieve level
    lp.le({t: 1}, F(1, 2))
    # failure exponent e <= each term's saving
    lp.le({e: 1, h: -2 * j}, 0)                          # E_b
    lp.le({e: 1, h: 2, y: 1, mu: 1}, F(1, 2))            # E_c2
    if dyadic and not detector:
        # Lemma D (dyadic over d = ord_q 2; Opus critic finding 6): no E_a term,
        # E_c1 -> n (A/sqrt P + A^2/P) with A = H^3 m
        lp.le({e: 1, h: 3, mu: 1}, F(1, 4))
        lp.le({e: 1, h: 6, mu: 2}, F(1, 2))
        return
    lp.le({e: 1, t: -1, h: 3, mu: 1}, 0)                 # E_c1
    if detector:
        lp.le({e: 1, t: 3}, 1)                           # C(T) coincidences
    else:
        lp.le({e: 1, t: 2}, F(1, 2))                     # E_a


def single_fallback(j=1, detector: bool = False) -> Tuple[F, Dict[str, F]]:
    """min over all stage parameters of max(cost, 1/5 - e) [D' then HH]."""
    j = F(j)
    lp = LP()
    z, e = lp.var("z"), lp.var("e")
    _stage(lp, "", {z: F(1)}, e, j, detector, t_shared=None)
    lp.le({e: -1, z: -1}, -ONE_FIFTH)                    # 1/5 - e <= z
    if detector:
        lp.le({"t": F(1, 2), z: -1}, 0)                  # detector cost t/2 <= z
    return lp.minimise(z)


def ladder(k: int, j=1, detector: bool = False, dyadic: bool = False) -> Tuple[F, Dict[str, F]]:
    """k DCL stages then HH.  E[T] exponent z with stage i cost <= z + e_{i-1}
    (it is reached on at most N^{-e_{i-1}} of S_n), HH: 1/5 - e_k <= z.
    With the detector, one T is shared (it is run once, cost t/2 <= z)."""
    j = F(j)
    lp = LP()
    z = lp.var("z")
    tsh = lp.var("t") if detector else None
    prev = None
    for i in range(1, k + 1):
        e = lp.var(f"e{i}")
        C = {z: F(1)} if prev is None else {z: F(1), prev: F(1)}
        _stage(lp, str(i), C, e, j, detector, t_shared=tsh, dyadic=dyadic)
        prev = e
    lp.le({prev: -1, z: -1}, -ONE_FIFTH)
    if detector:
        lp.le({tsh: F(1, 2), z: -1}, 0)
    return lp.minimise(z)


# ---------------------------------------------------------------------------
# 4. Closed forms (checked against the LPs in tests/test_avgtime.py)
# ---------------------------------------------------------------------------
def single_fallback_closed(j=1) -> F:
    """1/6 + 1/(30 + 36 j): balance 1/6 + 5h/3 = 1/5 - 2 j h."""
    return ONE_SIXTH + F(1, 30 + 36 * F(j))


LADDER_NO_DETECTOR = F(19, 110)    # inf over k of ladder(k, detector=False): D' one-threshold count
LADDER_NO_DETECTOR_DYADIC = F(17, 100)   # inf over k of ladder(k, detector=False, dyadic=True)


# ---------------------------------------------------------------------------
# 5. The ladder schedule of Theorem AV (integers, n-only)
# ---------------------------------------------------------------------------
def icbrt_like_root(x: int, num: int, den: int) -> int:
    """ceil(x^(num/den)) exactly for integers x >= 1."""
    target = x ** num
    lo, hi = 1, 1
    while hi ** den < target:
        hi *= 2
    while lo < hi:
        mid = (lo + hi) // 2
        if mid ** den >= target:
            hi = mid
        else:
            lo = mid + 1
    return lo


def ladder_schedule(n: int) -> List[int]:
    """H_1 = ceil(sqrt n), H_{i+1} = ceil(H_i^{6/5}), stop at H_max = n^2 * 2^{ceil(n/60)}
    (last step truncated).  All stages satisfy H_{i+1}^5 <= 2^5 H_i^6 (i.e.
    H_{i+1}^{5/3} <= 2^{5/3} H_i^2), which is all the expectation bound uses."""
    Hmax = (n * n) << (-(-n // 60))
    from math import isqrt
    Hs = [isqrt(n - 1) + 1]
    while Hs[-1] < Hmax:
        Hs.append(min(Hmax, icbrt_like_root(Hs[-1], 6, 5)))
    return Hs


def ladder_bound_terms(n: int, c: int = 0) -> List[F]:
    """Per-stage expectation terms / (n^c 2^{n/6}) of the proof's bound (constants 1):
    stage 1: H_1^{5/3};  stage i >= 2: H_i^{5/3} * min(1, n / H_{i-1}^2).
    Returned as exact Fractions of the 5th powers (term^3 to avoid roots)."""
    Hs = ladder_schedule(n)
    out = [F(Hs[0] ** 5)]                       # (H^{5/3})^3
    for a, b in zip(Hs, Hs[1:]):
        out.append(F(b ** 5) * min(F(1), F(n, a * a)) ** 3)
    return out


def min_stage_cost(e, j=1, detector: bool = False, dyadic: bool = False) -> F:
    """LP: least cost exponent of one DCL stage whose proved failure exponent is >= e."""
    lp = LP()
    C, ee = lp.var("C"), lp.var("e")
    _stage(lp, "", {C: F(1)}, ee, F(j), detector, dyadic=dyadic)
    lp.le({ee: -1}, -F(e))
    return lp.minimise(C)[0]


def stage_cost_closed(e) -> F:
    """Closed form of min_stage_cost(e) without detector (third branch: avgtime_1 item 1)."""
    e = F(e)
    if e <= F(1, 46):
        return ONE_SIXTH + 5 * e / 6
    if e <= F(1, 14):
        return F(1, 8) + 11 * e / 4
    return 9 * e / 2


def stage_cost_closed_dyadic(e) -> F:
    """Closed form of min_stage_cost(e, dyadic=True), valid for 0 <= e <= 1/20 (h = e/2)."""
    e = F(e)
    if e <= F(1, 40):
        return ONE_SIXTH + 5 * e / 6
    return F(1, 8) + 5 * e / 2


def dprime_dyadic_sup() -> F:
    """sup eps for which D' (unchanged algorithm) meets O(N^-eps) with Lemma D:
    E_c1' = 3h + mu - 1/4 = -1/12 + 7 eps/3, needs + eps < 0  ->  eps < 1/40;
    E_c1'' = 2(3h + mu) - 1/2 = -1/6 + 14 eps/3 -> eps < 1/34;  E_c2 -> 1/16; Y > H^2 -> 1/7."""
    cands = []
    for f in (lambda e: 3 * dprime_exponents(e)["h"] + dprime_exponents(e)["mu"] - F(1, 4) + e,
              lambda e: 2 * (3 * dprime_exponents(e)["h"] + dprime_exponents(e)["mu"]) - F(1, 2) + e):
        a = f(F(0)); b = f(F(1)) - a
        cands.append(sup_eps(a, b))
    return min(cands + [F(1, 16), F(1, 7)])
