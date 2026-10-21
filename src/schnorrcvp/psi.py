"""Smooth-number counts with rigorous integer brackets, and the heuristic fac-relation count model.

Agent `schnorrcvp`.  MODEL code (predictions), not operational.

Exact part (proved brackets)
----------------------------
For weights L_i^- = floor(K ln p_i) <= K ln p_i <= ceil(K ln p_i) = L_i^+ and an integer threshold,
    C^+(T) := #{e in N^n : sum e_i L_i^+ <= T},   C^-(T) := #{e : sum e_i L_i^- <= T}
satisfy  C^+(floor(K ln x)) <= Psi(x, p_n) <= C^-(ceil(K ln x))  [every smooth m <= x has
sum e_i K ln p_i <= K ln x, hence sum e_i L_i^- <= K ln x; and conversely].  The DP uses Python ints
(exact).  The integers floor/ceil(K ln p) are obtained from mpmath at 60 digits, which determines them
exactly unless K ln p is within 1e-50 of an integer (checked: assertion).

Heuristic part (model)
----------------------
Heuristic H_S: for v, w P-smooth, u = vN + w is P-smooth with probability d(vN), the local density of
P-smooth integers near vN, independently of (v, w).  Then the expected number of fac-relations with
ln v in a set S_v and ln|w| in S_w (|w| < vN) is
    E = sum_{v in S_v} sum_{w in S_w} 2 * d(vN)       (factor 2: both signs of w)
Densities and sums are evaluated in mpmath (floating point used only for these model numbers).
"""
from __future__ import annotations

from typing import List, Sequence, Tuple

import mpmath as mp

mp.mp.dps = 60


def log_weights(P: Sequence[int], K: int) -> Tuple[List[int], List[int]]:
    lo, hi = [], []
    for p in P:
        x = mp.mpf(K) * mp.log(p)
        f = int(mp.floor(x))
        assert abs(x - f) > mp.mpf(10) ** -50 and abs(x - f - 1) > mp.mpf(10) ** -50
        lo.append(f)
        hi.append(f + 1)
    return lo, hi


def exact_count_by_weight(L: Sequence[int], T: int) -> List[int]:
    """c[t] = #{e in N^n : sum e_i L_i = t} for 0 <= t <= T (exact Python ints)."""
    c = [0] * (T + 1)
    c[0] = 1
    for Li in L:
        if Li > T:
            continue
        for t in range(Li, T + 1):
            c[t] += c[t - Li]
    return c


def cumulative(c: Sequence[int]) -> List[int]:
    out = []
    s = 0
    for x in c:
        s += x
        out.append(s)
    return out


class PsiBracket:
    """Exact bracket tables for Psi(x, p_n) with ln x <= lnmax, resolution 1/K."""

    def __init__(self, P: Sequence[int], K: int, lnmax: float):
        self.P = list(P)
        self.K = K
        self.T = int(K * lnmax) + 2
        lo, hi = log_weights(P, K)
        self.c_lo = exact_count_by_weight(lo, self.T)  # floor weights -> upper bounds
        self.c_hi = exact_count_by_weight(hi, self.T)  # ceil weights  -> lower bounds
        self.C_lo = cumulative(self.c_lo)
        self.C_hi = cumulative(self.c_hi)

    def psi_bounds(self, lnx) -> Tuple[int, int]:
        """(lower, upper) integer bounds on Psi(e^lnx, p_n) [proved]."""
        t = mp.mpf(self.K) * mp.mpf(lnx)
        tl = int(mp.floor(t))
        tu = int(mp.ceil(t))
        assert tu <= self.T
        return self.C_hi[max(tl, 0)], self.C_lo[tu]

    def psi_mid(self, lnx):
        a, b = self.psi_bounds(lnx)
        return (mp.mpf(a) + mp.mpf(b)) / 2

    def density(self, lnx, h=0.05, which="mid"):
        """Local density of P-smooth integers on [x, x e^h] (model; brackets propagated if asked)."""
        lnx = mp.mpf(lnx)
        a0, b0 = self.psi_bounds(lnx)
        a1, b1 = self.psi_bounds(lnx + h)
        width = mp.e ** lnx * (mp.e ** h - 1)
        if which == "lower":
            return max(mp.mpf(0), mp.mpf(a1 - b0)) / width
        if which == "upper":
            return mp.mpf(b1 - a0) / width
        return (mp.mpf(a1 + b1) / 2 - mp.mpf(a0 + b0) / 2) / width


def binned_counts(pb: PsiBracket, lnmax: float, h: float, which="mid"):
    """Counts of P-smooth m >= 1 with ln m in [k h, (k+1) h) for k = 0, 1, ...  (bin 0 contains m = 1).

    which = "mid" uses the midpoint of the [proved] bracket at each edge (model value)."""
    def edge(x):
        a, b = pb.psi_bounds(x)
        if which == "lower":
            return mp.mpf(a)
        if which == "upper":
            return mp.mpf(b)
        return (mp.mpf(a) + mp.mpf(b)) / 2
    out = []
    prev = mp.mpf(0)  # Psi(e^0 - ) := 0, so bin 0 counts m = 1 too
    k = 0
    while (k + 1) * h <= lnmax + 1e-12:
        cur = edge((k + 1) * h)
        out.append(max(cur - prev, mp.mpf(0)))
        prev = cur
        k += 1
    return out


def expected_relations(P: Sequence[int], lnN: float, *, K: int = 100, h: float = 0.05,
                       lnv_max: float = None, lnw_max_rel: float = 0.0, close_delta: float = None):
    """Heuristic H_S expected number of fac-relations (u, v, w), u = vN + w, u, v, |w| P-smooth.

    Ranges: ln v in [0, lnv_max]; |w| <= vN * e^{lnw_max_rel} (lnw_max_rel <= 0; 0 means |w| < vN).
    If close_delta is given, restrict to |w| <= vN * close_delta (i.e. |u/(vN) - 1| <= close_delta).
    Returns dict with total E and the profile over ln v.
    """
    lnN = mp.mpf(lnN)
    if lnv_max is None:
        lnv_max = 3 * float(lnN)
    lnmax = float(lnv_max + lnN) + 2.0
    pb = PsiBracket(P, K, lnmax)
    vb = binned_counts(pb, lnv_max, h)
    total = mp.mpf(0)
    prof = []
    for k, nv in enumerate(vb):
        if nv == 0:
            prof.append((k * h, mp.mpf(0)))
            continue
        lnv = (k + mp.mpf(0.5)) * h
        lnvN = lnv + lnN
        cap = lnvN + lnw_max_rel
        if close_delta is not None:
            cap = min(cap, lnvN + mp.log(close_delta))
        if cap < 0:
            prof.append((float(lnv), mp.mpf(0)))
            continue
        nw = pb.psi_mid(cap)
        d = pb.density(lnvN - h / 2, h)
        contrib = nv * 2 * nw * d
        total += contrib
        prof.append((float(lnv), contrib))
    return {"E": total, "profile": prof}


def expected_in_box(N: int, P: Sequence[int], V: int, W: int, *, K: int = 1000, h: float = 0.02):
    """H_S prediction for the E1 box: sum over P-smooth v <= V of 2 * Psi(W) * d(vN).

    Psi(W) and the list of v are exact; d(vN) is the local density from the bracket tables on
    [vN e^{-h/2}, vN e^{h/2}] (midpoint of the proved bracket).  Also returns the bracket-induced
    lower/upper envelope of the prediction."""
    from .enum_rel import smooth_upto
    vs = smooth_upto(V, P)
    lnmax = float(mp.log(mp.mpf(V) * N)) + 1.0
    pb = PsiBracket(P, K, lnmax)
    nw = len(smooth_upto(W, P))
    E = mp.mpf(0)
    Elo = mp.mpf(0)
    Ehi = mp.mpf(0)
    for v in vs:
        x = mp.log(mp.mpf(v) * N) - h / 2
        E += 2 * nw * pb.density(x, h)
        Elo += 2 * nw * pb.density(x, h, which="lower")
        Ehi += 2 * nw * pb.density(x, h, which="upper")
    return {"E": E, "E_lo": Elo, "E_hi": Ehi, "n_v": len(vs), "n_w": nw}


def _prime_support(x: int, P: Sequence[int]):
    return tuple(p for p in P if x % p == 0)


def expected_primitive_in_box(N: int, P: Sequence[int], V: int, W: int, *, K: int = 1000, h: float = 0.02):
    """Refined heuristic H_S' for PRIMITIVE relations (gcd(v, w) = 1, hence u, v, w pairwise coprime).

    For coprime smooth (v, w), u = vN + w is automatically coprime to Q = {p in P : p | vw} and its residues
    mod the other primes are treated as uniform (CRT heuristic).  So
        Pr[u smooth] = dens{P-smooth m near x, gcd(m, Q) = 1} / prod_{q in Q} (1 - 1/q),
    and the numerator is computed by inclusion-exclusion from the P-smooth local density:
        dens{.., gcd(m,Q)=1}(x) = sum_{S subset Q} (-1)^{|S|} d_P(x / prod S) / prod S.
    Model numerics in float64 (a prediction, not an operational routine); w grouped by prime support.
    """
    from collections import Counter
    from itertools import combinations
    from .enum_rel import smooth_upto
    vs = smooth_upto(V, P)
    ws = smooth_upto(W, P)
    lnmax = float(mp.log(mp.mpf(V) * N)) + 1.0
    pb = PsiBracket(P, K, lnmax)
    grid = 1.0 / 500
    cache = {}

    def dens(lnx):
        key = int(round(lnx / grid))
        if key not in cache:
            cache[key] = float(pb.density(mp.mpf(key) * grid - h / 2, h))
        return cache[key]

    supp_w = Counter(_prime_support(w, P) for w in ws)
    E = 0.0
    import math
    for v in vs:
        lnx = math.log(v) + math.log(N)
        sv = set(_prime_support(v, P))
        for sw, cnt in supp_w.items():
            if sv.intersection(sw):
                continue
            Q = tuple(sorted(sv.union(sw)))
            num = 0.0
            for r in range(len(Q) + 1):
                for S in combinations(Q, r):
                    m = 1
                    for q in S:
                        m *= q
                    num += (-1) ** r * dens(lnx - math.log(m)) / m
            den = 1.0
            for q in Q:
                den *= (1 - 1.0 / q)
            E += 2 * cnt * num / den
    return {"E": E, "n_v": len(vs), "n_w": len(ws)}
