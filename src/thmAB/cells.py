"""thmAB: exact checks of the algebraic cell device (notes/thmAB_results.md §§3-5).

Everything is exact (integers / Fractions).  These functions verify, on small instances, the finite statements that
the proof uses:
  - Lemma CL (gcd-class bijection on L' = {x1 = kappa x2 mod Delta'});
  - Lemma CF (cell factorisation: c12rat's Lemma CP with the coupling set {r1 = kappa r2}), in the Parseval form;
  - Lemma PG (parallelogram geometry of the DW pair law: widths, longest chords, fibre masses);
  - Lemma BL(ii) (boundary rectangles of a convex region cover area <= (16h + 32h^2)|Pi|).
No factoring code; nothing here is operational.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from math import gcd
from typing import Dict, List, Sequence, Tuple

Pt = Tuple[Fr, Fr]


# ----------------------------------------------------------------------------- arithmetic helpers
def phi(n: int) -> int:
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def is_smooth(n: int, y: int) -> bool:
    n = abs(n)
    if n == 0:
        return False
    p = 2
    while p <= y and n > 1:
        while n % p == 0:
            n //= p
        p += 1
    return n == 1


# ----------------------------------------------------------------------------- Lemma CL
def class_bijection(Dp: int, kappa: int, X: int) -> Dict[str, int]:
    """Enumerate L' cap [-X, X]^2 and check: gcd(x1, Dp) = gcd(x2, Dp) =: delta, and x -> x/delta is a bijection onto
    {(m1, m2) : (m_i, M) = 1, m1 = kappa m2 (mod M)} cap [-X/delta, X/delta]^2, M = Dp/delta, for every delta | Dp."""
    assert gcd(kappa, Dp) == 1
    images: Dict[int, set] = {}
    n_pts = 0
    for x2 in range(-X, X + 1):
        for x1 in range(-X, X + 1):
            if (x1 - kappa * x2) % Dp:
                continue
            n_pts += 1
            d1, d2 = gcd(x1, Dp), gcd(x2, Dp)
            assert d1 == d2, (x1, x2, d1, d2)
            images.setdefault(d1, set()).add((x1 // d1, x2 // d1))
    for delta in [t for t in range(1, Dp + 1) if Dp % t == 0]:
        M = Dp // delta
        target = set()
        for m2 in range(-(X // delta), X // delta + 1):
            for m1 in range(-(X // delta), X // delta + 1):
                if gcd(m1, M) == 1 and gcd(m2, M) == 1 and (m1 - kappa * m2) % M == 0:
                    target.add((m1, m2))
        assert images.get(delta, set()) == target, (Dp, kappa, delta)
    return {"points": n_pts, "classes": len(images)}


# ----------------------------------------------------------------------------- Lemma CF (Parseval form)
def cf_check(Dp: int, kappa: int, delta: int, I1: Tuple[int, int], I2: Tuple[int, int], y: int) -> Dict[str, Fr]:
    """Rectangle I1 x I2 (closed integer intervals), class delta.  a(m) = 1[delta m in I1, |delta m| y-smooth, (m, M) = 1].
    Checks |#both - A B / phi(M)| <= (1/phi(M)) * sqrt(Ea Eb), with Ea = sum_{chi != chi0} |a^(chi)|^2
    = phi(M) * sum_r (a(r) - A/phi)^2 (Parseval on (Z/M)^x).  Comparison done on squares, exactly."""
    M = Dp // delta
    ph = phi(M)
    units = [r for r in range(M) if gcd(r, M) == 1] if M > 1 else [0]

    def seq(I):
        lo, hi = I
        out = {}
        for x in range(lo, hi + 1):
            if x % delta:
                continue
            m = x // delta
            if gcd(m, M) != 1 or not is_smooth(x, y) or gcd(x, Dp) != delta:
                continue
            out[m] = 1
        return out

    a, b = seq(I1), seq(I2)
    both = sum(1 for m1 in a for m2 in b if (m1 - kappa * m2) % M == 0)
    A, B = len(a), len(b)

    def energy(s):
        cls = {r: 0 for r in units}
        for m in s:
            cls[m % M if M > 1 else 0] += 1
        mean = Fr(len(s), ph)
        return ph * sum((Fr(v) - mean) ** 2 for v in cls.values())

    Ea, Eb = energy(a), energy(b)
    lhs = Fr(both) - Fr(A * B, ph)
    ok = lhs * lhs * ph * ph <= Ea * Eb
    assert ok, (Dp, kappa, delta, I1, I2, lhs, Ea, Eb)
    return {"both": Fr(both), "main": Fr(A * B, ph), "Ea": Ea, "Eb": Eb}


# ----------------------------------------------------------------------------- Lemma PG (parallelogram geometry)
def v_of(P: Tuple[int, int], d: int) -> List[int]:
    a, b = P
    return [a ** (d - 1 - i) * b ** i for i in range(d)]


def mm_choice(P1, P2, d: int) -> Tuple[int, int, int]:
    """(j, k, m_jk) maximising |m_jk| over Lemma MM's five candidate pairs."""
    v1, v2 = v_of(P1, d), v_of(P2, d)
    best = None
    for j, k in [(0, 1), (d - 2, d - 1), (0, d - 1), (0, d - 2), (1, d - 1)]:
        m = v1[j] * v2[k] - v1[k] * v2[j]
        if best is None or abs(m) > abs(best[2]):
            best = (j, k, m)
    return best


def parallelogram(P1, P2, d: int, H: int, G1: int, G2: int, q: int = 1):
    """Vertices of Pi_0 = q^{-1} diag(G1, G2) M_jk [-H, H+1)^2 (closed hull), and the data (j, k, m_jk)."""
    j, k, m = mm_choice(P1, P2, d)
    v1, v2 = v_of(P1, d), v_of(P2, d)
    e = (Fr(G1 * v1[j], q), Fr(G2 * v2[j], q))
    f = (Fr(G1 * v1[k], q), Fr(G2 * v2[k], q))
    vs = []
    for s, t in [(-H, -H), (H + 1, -H), (H + 1, H + 1), (-H, H + 1)]:
        vs.append((s * e[0] + t * f[0], s * e[1] + t * f[1]))
    return vs, (j, k, m)


def area(poly: Sequence[Pt]) -> Fr:
    n = len(poly)
    return abs(sum(poly[i][0] * poly[(i + 1) % n][1] - poly[(i + 1) % n][0] * poly[i][1] for i in range(n))) / 2


def chord(poly: Sequence[Pt], coord: int, level: Fr) -> Fr:
    """Length of {x in poly : x_coord = level} (convex poly), measured along the other coordinate."""
    o = 1 - coord
    pts = []
    n = len(poly)
    for i in range(n):
        p, r = poly[i], poly[(i + 1) % n]
        if p[coord] == r[coord]:
            if p[coord] == level:
                pts += [p[o], r[o]]
            continue
        lo, hi = min(p[coord], r[coord]), max(p[coord], r[coord])
        if lo <= level <= hi:
            t = (level - p[coord]) / (r[coord] - p[coord])
            pts.append(p[o] + t * (r[o] - p[o]))
    return (max(pts) - min(pts)) if pts else Fr(0)


def longest_chord(poly: Sequence[Pt], coord: int) -> Fr:
    """Longest chord along coordinate 1-coord (lines x_coord = const); attained at a vertex level (piecewise linear)."""
    return max(chord(poly, coord, v[coord]) for v in poly)


def width(poly: Sequence[Pt], coord: int) -> Fr:
    return max(v[coord] for v in poly) - min(v[coord] for v in poly)


def pg_check(P1, P2, d: int, H: int, G1: int, G2: int) -> Dict[str, Fr]:
    """Lemma PG (q = 1): (a) W_i >= (2H+1)|G_i||P_i|^{d-1}/2; (b) l_r >= |Pi|/W_2 >= (2H+1)|G1||P1|^{d-1}/(8 kappa),
    kappa := |P1|_inf |P2|_inf / |det(P1, P2)|; symmetric for l_c."""
    poly, (j, k, m) = parallelogram(P1, P2, d, H, G1, G2)
    n1 = max(abs(P1[0]), abs(P1[1]))
    n2 = max(abs(P2[0]), abs(P2[1]))
    Det = abs(P1[0] * P2[1] - P1[1] * P2[0])
    kap = Fr(n1 * n2, Det)
    W1, W2 = width(poly, 0), width(poly, 1)
    A = area(poly)
    lr = longest_chord(poly, 1)        # chords parallel to the x1-axis (lines x2 = const)
    lc = longest_chord(poly, 0)        # chords parallel to the x2-axis (lines x1 = const)
    s1 = Fr((2 * H + 1) * abs(G1) * n1 ** (d - 1))
    s2 = Fr((2 * H + 1) * abs(G2) * n2 ** (d - 1))
    assert abs(m) * 4 >= Det * (n1 * n2) ** (d - 2), "Lemma MM"
    assert W1 >= s1 / 2 and W2 >= s2 / 2, "PG(a)"
    assert lr * W2 >= A and lc * W1 >= A, "area <= chord x extent"
    assert lr >= s1 / (8 * kap) and lc >= s2 / (8 * kap), ("PG(b)", lr, s1 / (8 * kap))
    return {"W1": W1, "W2": W2, "area": A, "lr": lr, "lc": lc, "kappa": kap, "s1": s1, "s2": s2}


# ----------------------------------------------------------------------------- Lemma BL(ii): boundary rectangles
def _inside(poly: Sequence[Pt], p: Pt) -> bool:
    n = len(poly)
    sgn = 0
    for i in range(n):
        a, b = poly[i], poly[(i + 1) % n]
        cr = (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
        if cr != 0:
            s = 1 if cr > 0 else -1
            if sgn == 0:
                sgn = s
            elif s != sgn:
                return False
    return True


def _clip(poly: List[Pt], x0: Fr, x1: Fr, y0: Fr, y1: Fr) -> List[Pt]:
    def clip_edge(pts, inside, inter):
        out = []
        n = len(pts)
        for i in range(n):
            cur, prev = pts[i], pts[i - 1]
            if inside(cur):
                if not inside(prev):
                    out.append(inter(prev, cur))
                out.append(cur)
            elif inside(prev):
                out.append(inter(prev, cur))
        return out

    def ix(c):
        return lambda p, r: (c, p[1] + (r[1] - p[1]) * (c - p[0]) / (r[0] - p[0]))

    def iy(c):
        return lambda p, r: (p[0] + (r[0] - p[0]) * (c - p[1]) / (r[1] - p[1]), c)

    pts = list(poly)
    for ins, inter in [(lambda p: p[0] >= x0, ix(x0)), (lambda p: p[0] <= x1, ix(x1)),
                       (lambda p: p[1] >= y0, iy(y0)), (lambda p: p[1] <= y1, iy(y1))]:
        if not pts:
            return []
        pts = clip_edge(pts, ins, inter)
    return pts


def boundary_area_check(poly: Sequence[Pt], h: Fr) -> Dict[str, Fr]:
    """Tile R^2 by z1 x z2 rectangles (z1 = h l_r, z2 = h l_c), anchored at the origin; rectangles meeting Pi (positive
    intersection area) but not contained in it are boundary rectangles.  Checks: their total area <= (16h + 32h^2)|Pi|."""
    lr, lc = longest_chord(poly, 1), longest_chord(poly, 0)
    z1, z2 = h * lr, h * lc
    xs = [v[0] for v in poly]
    ys = [v[1] for v in poly]
    i0, i1 = int((min(xs) / z1).__floor__()) - 1, int((max(xs) / z1).__ceil__()) + 1
    j0, j1 = int((min(ys) / z2).__floor__()) - 1, int((max(ys) / z2).__ceil__()) + 1
    A = area(poly)
    bd = Fr(0)
    n_in = n_bd = 0
    for i in range(i0, i1):
        for jj in range(j0, j1):
            x0, x1, y0, y1 = i * z1, (i + 1) * z1, jj * z2, (jj + 1) * z2
            corners = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
            if all(_inside(poly, c) for c in corners):
                n_in += 1
                continue
            cl = _clip(list(poly), x0, x1, y0, y1)
            if len(cl) >= 3 and area(cl) > 0:
                n_bd += 1
                bd += z1 * z2
    assert bd <= (16 * h + 32 * h * h) * A, (bd / A, h)
    return {"ratio": bd / A, "bound": 16 * h + 32 * h * h, "inside": Fr(n_in), "boundary": Fr(n_bd)}


# ----------------------------------------------------------------------------- Lemma PG(c): fibre masses
def fibre_mass_check(P1, P2, d: int, H: int, G1: int, G2: int, lam_d: int, Dp: int) -> Dict[str, Fr]:
    """Uniform law on L' cap Pi_0 (L' = {x1 = lam_d x2 mod Dp}, Pi_0 closed hull).  Finite forms of PG(c):
    Z >= I - E2 and max fibre mass <= (l_c/Dp + 1)/(I - E2), and Pr[s | x1] <= (I/s + E1)/(I - E2) for s | Dp, where
    I = |Pi|/Dp, E1 = 2 l_c/Dp + W1/s + 1, E2 = 2 l_c/Dp + W1 + 1 (concavity of the chord function).  As l_c/Dp and W1
    tend to infinity these give max mass <= 3/W1 and Pr[s | x1] <= (1 + o(1))/s."""
    poly, _ = parallelogram(P1, P2, d, H, G1, G2)
    W1 = width(poly, 0)
    lc = longest_chord(poly, 0)
    xs = [v[0] for v in poly]
    lo, hi = int(min(xs).__floor__()), int(max(xs).__ceil__())
    counts: Dict[int, int] = {}
    for x1 in range(lo, hi + 1):
        # x2 range on the fibre
        ys = []
        n = len(poly)
        for i in range(n):
            p, r = poly[i], poly[(i + 1) % n]
            if p[0] == r[0]:
                if p[0] == x1:
                    ys += [p[1], r[1]]
                continue
            a, b = min(p[0], r[0]), max(p[0], r[0])
            if a <= x1 <= b:
                t = (x1 - p[0]) / (r[0] - p[0])
                ys.append(p[1] + t * (r[1] - p[1]))
        if not ys:
            continue
        y0, y1 = min(ys), max(ys)
        c = 0
        # x2 = x2_0 (mod Dp) with lam_d x2_0 = x1
        inv = pow(lam_d, -1, Dp) if Dp > 1 else 0
        r0 = (inv * x1) % Dp if Dp > 1 else 0
        start = int(y0.__ceil__())
        start += (r0 - start) % Dp if Dp > 1 else 0
        x2 = start
        while x2 <= y1:
            c += 1
            x2 += Dp if Dp > 1 else 1
        if c:
            counts[x1] = c
    Z = sum(counts.values())
    A = area(poly)
    I = A / Dp
    E2 = 2 * lc / Dp + W1 + 1
    mx = Fr(max(counts.values()), Z)
    # finite form: each fibre has <= l_c/Dp + 1 points and Z >= I - E2 (concave chord function)
    assert Z >= I - E2, ("Z lower bound", Z, I - E2)
    assert mx <= (lc / Dp + 1) / (I - E2), ("PG(c)", mx)
    out = {"Z": Fr(Z), "max_mass_W1": mx * W1, "I": I}
    for s in [t for t in range(2, Dp + 1) if Dp % t == 0]:
        ps = Fr(sum(v for x, v in counts.items() if x % s == 0), Z)
        E1 = 2 * lc / Dp + W1 / s + 1
        assert ps <= (I / s + E1) / (I - E2), ("local law", s, ps)
        out[f"Pr[{s}|x1]*{s}"] = ps * s
    return out
