"""topblock lattices: does adding the ladder's WINDOW congruences (one block or several) to the
Coppersmith divisor lattice extend its range beyond N^{1/4}?   HARNESS/DIAGNOSTIC measurement code.

Unknown: h0 = p - x (0 <= h0 < X), where x is the block start.  The harness knows p only to plant
the block start x = p - h0 and to report; every lattice is built from (N, x, public parameters)
and the Howgrave-Graham test uses the public bound p >= x.

Lattices (all lower triangular in a monomial order, reduced with fplll LLL (0.99, 0.51)):

U(m, n)       univariate divisor lattice: N^{m-i} f^i (i < m), h^a f^m;  f = x + h.          [May10 Thm 3]
MIX(...)      divisor rows times window rows.  Window b = (Delta_b, j_b): block start x_b = x - Delta_b,
              h_b = h + Delta_b, D_b = x_b^{j_b+1},
              w_b(h, r_b) = r_b - (num_{j_b}(h_b) mod D_b, coefficientwise in h),
              root r_b0 = num_{j_b}(h_b0) - q0 D_b = -q0 (-h_b0)^{j_b+1}   (identity T2 with k = 0),
              |r_b0| <= Y_b := ceil(N/x) * (X + |Delta_b|)^{j_b+1}.
              Rows, for t = (t_b) with t_b <= m2_b and e + sum_b j_b t_b <= E:
                 h^{e - s} f^{s} N^{m1 - s} prod_b w_b^{t_b} D_b^{m2_b - t_b},   s = min(e, m1);
              they vanish mod p^{m1} prod D_b^{m2_b} at the root.
HM(m, E)      Herrmann-May style split p = x0 + H u + v (block index u, offset v):
              u^a f^b N^{(m-b)+}, a + b <= E, f = v + H u + x0.

TRIVIALITY (notes §3, Prop T): substituting r_b -> r_b(h) := -N(-(h+Delta_b))^{j_b+1}/(x+h) and
multiplying by (x+h)^T maps a MIX vector g to a univariate G in prod D_b^{m2_b} * (N, x+h)^{m1+T}.
G == 0 iff g is a multiple of the curve relations F_b = (x+h) r_b + N(-(h+Delta_b))^{j_b+1}; such g
vanish on the whole curve and locate nothing.  A MIX success must come from a vector with G != 0.
"""
from __future__ import annotations

import math
import time
from itertools import product
from math import gcd
from typing import Dict, List, Optional, Sequence, Tuple

import flint
from fpylll import IntegerMatrix, LLL, SVP

from src.topblock.core import S_poly


# ---------------------------------------------------------------- tiny sparse polynomial algebra

Poly = Dict[Tuple[int, ...], int]   # key = (e, t_1, ..., t_k): h^e prod r_b^{t_b}


def pmul(a: Poly, b: Poly) -> Poly:
    out: Poly = {}
    for ka, ca in a.items():
        for kb, cb in b.items():
            k = tuple(x + y for x, y in zip(ka, kb))
            out[k] = out.get(k, 0) + ca * cb
    return {k: v for k, v in out.items() if v}


def ppow(a: Poly, e: int, nvars: int) -> Poly:
    out: Poly = {(0,) * nvars: 1}
    for _ in range(e):
        out = pmul(out, a)
    return out


def pscale(a: Poly, c: int) -> Poly:
    return {k: v * c for k, v in a.items()} if c != 1 else dict(a)


def _num_shifted_coeffs(N: int, xb: int, j: int, Delta: int) -> List[int]:
    """Coefficients (in h) of num_j(h + Delta) for block start xb: N * S_j(xb; h + Delta)."""
    # S_j(xb; y) = sum_i (-y)^i xb^{j-i}; expand y = h + Delta
    coeffs = [0] * (j + 1)
    for i in range(j + 1):
        base = (-1) ** i * xb ** (j - i)
        for a in range(i + 1):
            coeffs[a] += base * math.comb(i, a) * Delta ** (i - a)
    return [N * c for c in coeffs]


# ---------------------------------------------------------------- lattice reduction helpers

def lll_rows(rows: List[List[int]]) -> List[List[int]]:
    A = IntegerMatrix.from_matrix(rows)
    LLL.reduction(A)
    return [[int(A[i, k]) for k in range(A.ncols)] for i in range(A.nrows)]


def sqnorm(v: Sequence[int]) -> int:
    return sum(c * c for c in v)


def log2_int(z: int) -> float:
    if z <= 0:
        return float("-inf")
    b = z.bit_length()
    if b <= 60:
        return math.log2(z)
    return (b - 60) + math.log2(z >> (b - 60))


def exact_shortest(red: List[List[int]]) -> List[int]:
    """Exact SVP by full (unpruned) enumeration on an LLL-reduced basis, mpfr GSO."""
    from fpylll import GSO, Enumeration, FPLLL
    bits = max(max(abs(c) for c in row).bit_length() for row in red)
    FPLLL.set_precision(max(120, 2 * bits + 64))
    A = IntegerMatrix.from_matrix(red)
    M = GSO.Mat(A, float_type="mpfr")
    M.update_gso()
    n = A.nrows
    enum = Enumeration(M)
    radius = M.get_r(0, 0) * 1.0000001
    sols = enum.enumerate(0, n, radius, 0)
    coords = [int(round(c)) for c in sols[0][1]]
    v = [0] * A.ncols
    for i, c in enumerate(coords):
        if c:
            for k in range(A.ncols):
                v[k] += c * red[i][k]
    if sqnorm(v) == 0 or sqnorm(v) > sqnorm(red[0]):
        v = list(red[0])
    return v


# ---------------------------------------------------------------- U(m, n)

def univariate_basis(N: int, x: int, m: int, n: int, X: int) -> Tuple[List[List[int]], List[int]]:
    f = flint.fmpz_poly([x, 1])
    rows, diag = [], []
    for i in range(n):
        if i < m:
            g = f ** i * N ** (m - i)
        else:
            g = flint.fmpz_poly([0] * (i - m) + [1]) * f ** m
        cs = [int(c) for c in g.coeffs()] + [0] * n
        rows.append([cs[k] * X ** k for k in range(n)])
        diag.append(rows[-1][i])
    return rows, diag


def univariate_attempt(N: int, x: int, m: int, n: int, X: int, use_svp: bool = False) -> Dict:
    rows, diag = univariate_basis(N, x, m, n, X)
    red = lll_rows(rows)
    thr2 = x ** (2 * m)
    good = [v for v in red if v and n * sqnorm(v) < thr2]
    found = None
    for v in good[:3]:
        cs = [v[k] // X ** k for k in range(n)]
        assert all(cs[k] * X ** k == v[k] for k in range(n))
        G = flint.fmpz_poly(cs)
        found = _divisor_from_roots(N, x, G, X)
        if found:
            break
    out = {"n": n, "m": m, "good": len(good), "success": bool(found),
           "log2_b1": 0.5 * log2_int(sqnorm(red[0])),
           "log2_thr": m * log2_int(x) - 0.5 * math.log2(n),
           "log2_det_n": sum(log2_int(abs(d)) for d in diag) / n}
    if use_svp:
        sv = exact_shortest(red)
        lam = 0.5 * log2_int(sqnorm(sv))
        out["log2_lambda1"] = lam
        out["svp_success_hg"] = bool(n * sqnorm([int(c) for c in sv]) < thr2)
        out["log2_gh"] = out["log2_det_n"] + 0.5 * math.log2(n / (2 * math.pi * math.e))
    return out


def _divisor_from_roots(N: int, x: int, G, X: int) -> Optional[int]:
    if G.degree() < 1:
        return None
    try:
        rts = G.roots()
    except Exception:
        return None
    for r, _ in rts:
        h = int(r)
        if 0 <= h <= X and x + h > 1:
            d = x + h
            g = gcd(d, N)
            if 1 < g < N:
                return g
    return None


# ---------------------------------------------------------------- MIX

def mixed_basis(N: int, x: int, windows: List[Tuple[int, int, int]], m1: int, E: int, X: int):
    """windows: list of (Delta_b, j_b, m2_b).  Returns (rows, monomials, scales, data)."""
    k = len(windows)
    nv = k + 1
    Ds, ws, Ys = [], [], []
    for b, (Delta, j, m2) in enumerate(windows):
        xb = x - Delta
        D = xb ** (j + 1)
        cs = [c % D for c in _num_shifted_coeffs(N, xb, j, Delta)]
        w: Poly = {}
        key_r = [0] * nv
        key_r[1 + b] = 1
        w[tuple(key_r)] = 1
        for a, c in enumerate(cs):
            if c:
                kk = [0] * nv
                kk[0] = a
                w[tuple(kk)] = w.get(tuple(kk), 0) - c
        Ds.append(D)
        ws.append(w)
        Ys.append(-(-N // x) * (X + abs(Delta)) ** (j + 1))
    f: Poly = {tuple([1] + [0] * k): 1, tuple([0] * nv): x}
    tvecs = [t for t in product(*[range(m2 + 1) for (_, _, m2) in windows])
             if sum(windows[b][1] * t[b] for b in range(k)) <= E]
    monos = []
    for t in tvecs:
        emax = E - sum(windows[b][1] * t[b] for b in range(k))
        for e in range(emax + 1):
            monos.append((e,) + tuple(t))
    monos.sort(key=lambda kk: (sum(kk[1:]), kk[1:], kk[0]))
    idx = {mo: i for i, mo in enumerate(monos)}
    wpow_cache = {}
    rows = []
    for mo in monos:
        e, t = mo[0], mo[1:]
        s = min(e, m1)
        g: Poly = {tuple([e - s] + [0] * k): N ** (m1 - s)}
        g = pmul(g, ppow(f, s, nv))
        for b in range(k):
            if t[b]:
                key = (b, t[b])
                if key not in wpow_cache:
                    wpow_cache[key] = ppow(ws[b], t[b], nv)
                g = pmul(g, wpow_cache[key])
            g = pscale(g, Ds[b] ** (windows[b][2] - t[b]))
        row = [0] * len(monos)
        for kk, c in g.items():
            row[idx[kk]] = c
        rows.append(row)
    scales = []
    for mo in monos:
        sc = X ** mo[0]
        for b in range(k):
            sc *= Ys[b] ** mo[1 + b]
        scales.append(sc)
    srows = [[row[i] * scales[i] for i in range(len(monos))] for row in rows]
    modulus_pub = x ** m1
    for b in range(k):
        modulus_pub *= Ds[b] ** windows[b][2]
    return srows, monos, scales, {"Ds": Ds, "Ys": Ys, "thr": modulus_pub}


def substitute_univariate(N: int, x: int, windows, monos, coeffs) -> "flint.fmpz_poly":
    """G(h) = (x+h)^T g(h, r_b(h)),  r_b(h) = -N(-(h+Delta_b))^{j_b+1}/(x+h),  T = max total r-degree."""
    k = len(windows)
    T = max((sum(mo[1:]) for mo, c in zip(monos, coeffs) if c), default=0)
    xh = flint.fmpz_poly([x, 1])
    rnum = []
    for (Delta, j, _) in windows:
        base = flint.fmpz_poly([Delta, 1])          # h + Delta
        rnum.append(-N * (-base) ** (j + 1))
    G = flint.fmpz_poly([0])
    for mo, c in zip(monos, coeffs):
        if not c:
            continue
        term = flint.fmpz_poly([0] * mo[0] + [c])
        tt = 0
        for b in range(k):
            if mo[1 + b]:
                term = term * rnum[b] ** mo[1 + b]
                tt += mo[1 + b]
        term = term * xh ** (T - tt)
        G = G + term
    return G


def mixed_attempt(N: int, x: int, windows, m1: int, E: int, X: int, root: Optional[Tuple] = None) -> Dict:
    t0 = time.time()
    srows, monos, scales, data = mixed_basis(N, x, windows, m1, E, X)
    n = len(monos)
    red = lll_rows(srows)
    thr2 = data["thr"] ** 2
    good = [v for v in red if any(v) and n * sqnorm(v) < thr2]
    nontriv = 0
    found = None
    for v in good:
        cs = []
        for i in range(n):
            q, rmd = divmod(v[i], scales[i])
            assert rmd == 0
            cs.append(q)
        G = substitute_univariate(N, x, windows, monos, cs)
        if G == 0:
            continue
        nontriv += 1
        if found is None:
            found = _divisor_from_roots(N, x, G, X)
    return {"n": n, "good": len(good), "good_nontrivial": nontriv, "success": bool(found),
            "log2_b1": 0.5 * log2_int(sqnorm(red[0])), "log2_thr": log2_int(data["thr"]) - 0.5 * math.log2(n),
            "secs": time.time() - t0}


def mixed_dim(windows, E: int) -> int:
    k = len(windows)
    cnt = 0
    for t in product(*[range(m2 + 1) for (_, _, m2) in windows]):
        s = sum(windows[b][1] * t[b] for b in range(k))
        if s <= E:
            cnt += E - s + 1
    return cnt


def F_relation_row_is_short(N: int, x: int, j: int, X: int, m1: int, m2: int) -> bool:
    """Diagnostic: is F * N^{m1-1} D^{m2-1} Howgrave-Graham short (it is in MIX when m1, m2 >= 1)?"""
    D = x ** (j + 1)
    Y = -(-N // x) * X ** (j + 1)
    v = [x * Y, X * Y, N * X ** (j + 1)]   # monomials r, h r, h^{j+1}
    c = N ** max(m1 - 1, 0) * D ** max(m2 - 1, 0)
    return 3 * sum((c * a) ** 2 for a in v) < (x ** m1 * D ** m2) ** 2


# ---------------------------------------------------------------- HM split (block index, offset)

def hm_attempt(N: int, x0: int, Hs: int, U: int, m: int, E: int) -> Dict:
    """p = x0 + Hs*u + v, 0 <= u < U, 0 <= v < Hs.  Rows u^a f^b N^{(m-b)+}, a + b <= E."""
    t0 = time.time()
    monos = [(a, b) for b in range(E + 1) for a in range(E + 1 - b)]   # u^a v^b
    monos.sort(key=lambda ab: (ab[1], ab[0]))
    idx = {mo: i for i, mo in enumerate(monos)}
    f = {(0, 1): 1, (1, 0): Hs, (0, 0): x0}
    rows = []
    fp = {0: {(0, 0): 1}}
    for b in range(1, E + 1):
        fp[b] = pmul(fp[b - 1], f)
    for (a, b) in monos:
        g = pmul({(a, 0): N ** max(m - b, 0)}, fp[b])
        row = [0] * len(monos)
        for kk, c in g.items():
            row[idx[kk]] = c
        rows.append(row)
    scales = [U ** a * Hs ** b for (a, b) in monos]
    srows = [[row[i] * scales[i] for i in range(len(monos))] for row in rows]
    n = len(monos)
    red = lll_rows(srows)
    thr2 = x0 ** (2 * m)
    good = [v for v in red if any(v) and n * sqnorm(v) < thr2]
    ctx = flint.fmpz_mpoly_ctx.get(("u", "v"), "lex")
    uu, vv = ctx.gens()
    polys = []
    for v in good:
        g = ctx.from_dict({monos[i]: v[i] // scales[i] for i in range(n) if v[i]})
        polys.append(g)
    found = None
    n_line = 0   # good vectors with a linear factor b(H u + v) + c (i.e. functions of z = H u + v)
    for g in polys:
        try:
            _, facs = g.factor()
        except Exception:
            facs = []
        hit = False
        for (fa, _) in facs:
            dd = fa.to_dict()
            if fa.total_degree() != 1:
                continue
            a = int(dd.get((1, 0), 0)); b = int(dd.get((0, 1), 0)); c = int(dd.get((0, 0), 0))
            if b != 0 and a == Hs * b and c % b == 0:
                hit = True
                d = x0 - c // b
                if 1 < d < N and N % d == 0 and found is None:
                    found = d
        n_line += hit
    for i1 in range(len(polys)):
        for i2 in range(i1 + 1, len(polys)):
            if found:
                break
            res = polys[i1].resultant(polys[i2], "v")
            if res == 0:
                continue
            coeffs = [0] * (res.degrees()[0] + 1)
            for (eu, ev), c in res.to_dict().items():
                coeffs[eu] += int(c)
            R = flint.fmpz_poly(coeffs)
            if R.degree() < 1:
                continue
            for (ur, _) in R.roots():
                u0 = int(ur)
                if not (0 <= u0 <= U):
                    continue
                gv = polys[i1].subs({"u": u0})
                cv = [0] * (gv.degrees()[1] + 1)
                for (eu, ev), c in gv.to_dict().items():
                    cv[ev] += int(c)
                V = flint.fmpz_poly(cv)
                if V.degree() < 1:
                    continue
                for (vr, _) in V.roots():
                    d = x0 + Hs * u0 + int(vr)
                    if 1 < d < N and N % d == 0:
                        found = d
        if found:
            break
    return {"n": n, "good": len(good), "good_line": n_line, "success": bool(found), "secs": time.time() - t0}
