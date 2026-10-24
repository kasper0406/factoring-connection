"""E3: the cell (Cauchy–Schwarz) decomposition of child-box pair covariances, rational block, toy scale.
DIAGNOSTIC (notes/c12rat_prereg.md).

For a label (q, r) and its child box, with g uniform on ALL primitive vectors of [M/2, M]^2 (exact enumeration)
and f(n) = 1[|n| y-smooth]:
  Pair(P', P'') = E_g f(G') f(G''),   S' = E_g f(G').
Cells of a value n (one-point, depending on the pair only through D = |det(P', P'')|):
  (log-interval of |n| with ratio 1 + h)  x  (v_2 <= 3, v_3 <= 2, v_5 <= 1, v_7 <= 1, capped)
  x (min(v_p(n), v_p(D)) for the primes 7 < p <= y dividing D).
psi'(cell) := E[f(G') | cell(G')]; then exactly
  Pair - S'S'' = Cov_nu(psi', psi'') + R,   Cov_nu(psi', psi'') := E_g[psi'(cell(G')) psi''(cell(G''))] - S'S'',
and |Cov_nu| <= sqrt(V' V''), V' := Var_g psi'(cell(G')) (Cauchy–Schwarz; the theorem's one-point bound).
R is the within-cell ("non-principal") part, which the theorem controls only on average over pairs/moduli (BDH).
"""
from __future__ import annotations

from math import log

import numpy as np

from src.pc3 import toy as pc3toy

SMALL = ((2, 3), (3, 2), (5, 1), (7, 1))


def all_primitive_g(M: int) -> np.ndarray:
    lo = (M + 1) // 2
    a = np.arange(lo, M + 1, dtype=np.int64)
    g0, g1 = np.meshgrid(a, a, indexing="ij")
    g0 = g0.ravel()
    g1 = g1.ravel()
    ok = np.gcd(g0, g1) == 1
    return np.stack([g0[ok], g1[ok]], axis=1)


def capped_v(n: np.ndarray, p: int, cap: int) -> np.ndarray:
    v = np.zeros(n.shape, dtype=np.int64)
    m = n.copy()
    for _ in range(cap):
        d = (m % p) == 0
        v += d
        m = np.where(d, m // p, m)
    return v


def base_cells(n: np.ndarray, h: float) -> np.ndarray:
    """interval index (log grid, ratio 1+h) combined with the capped small-prime type."""
    k = np.where(n > 0, np.floor(np.log(np.maximum(n, 1).astype(np.float64)) / log(1.0 + h)).astype(np.int64) + 1, 0)
    cid = k
    for p, cap in SMALL:
        cid = cid * (cap + 1) + capped_v(n, p, cap)
    return cid


def factor_mid(D: int, y: int):
    out = []
    for p in pc3toy.primes_upto(y):
        p = int(p)
        if p <= 7:
            continue
        if D % p == 0:
            k = 0
            m = D
            while m % p == 0:
                m //= p
                k += 1
            out.append((p, k))
    return out


def compact(c: np.ndarray):
    u, inv = np.unique(c, return_inverse=True)
    return inv, len(u)


def cell_mean(cells: np.ndarray, f: np.ndarray, k: int):
    """cells: nonnegative ids < k."""
    cnt = np.bincount(cells, minlength=k)
    s = np.bincount(cells, weights=f, minlength=k)
    return (s / np.maximum(cnt, 1))[cells]


def all_g(M: int) -> np.ndarray:
    """ALL integer g in [M/2, M]^2 (the theory's "version (a)" law; no primitivity)."""
    lo = (M + 1) // 2
    a = np.arange(lo, M + 1, dtype=np.int64)
    g0, g1 = np.meshgrid(a, a, indexing="ij")
    return np.stack([g0.ravel(), g1.ravel()], axis=1)


def local_ratios(smooth: np.ndarray, x: int, y: int) -> dict:
    """Primitivity factor per prime 7 < p <= y (notes/c12rat_prereg.md E3):
    A_p = Pr(smooth | p | n)/Pr(smooth), n uniform in [x/2, 2x];  B_p = (1 - A_p/p)/(1 - 1/p);
    ratio_p = (p+1)(2 A B + (p-1) B^2)/(A + p B)^2 = pair/single^2 for (n1, n2) mod p uniform on F_p^2 minus 0."""
    lo, hi = max(x // 2, 2), min(2 * x, len(smooth) - 1)
    seg = smooth[lo:hi]
    n = np.arange(lo, hi)
    base = seg.mean()
    out = {}
    for p in pc3toy.primes_upto(y):
        p = int(p)
        if p <= 7:
            continue
        A = seg[n % p == 0].mean() / base
        B = (1 - A / p) / (1 - 1 / p)
        out[p] = float(np.log((p + 1) * (2 * A * B + (p - 1) * B * B) / (A + p * B) ** 2))
    return out


def label_run(q: int, r: int, S: int, g: np.ndarray, smooth: np.ndarray, y: int, h: float, Kthin: float,
              prim_model: bool = False):
    e1, e2 = pc3toy.label_basis(q, r)
    box = pc3toy.child_box(e1, e2, S)
    ab = np.array([(b[2], b[3]) for b in box], dtype=np.int64)
    npts = len(ab)
    vals = np.abs(np.outer(g[:, 1], ab[:, 0]) - np.outer(g[:, 0], ab[:, 1]))  # (ng, npts)
    if vals.max() >= len(smooth):
        raise ValueError("table too small")
    F = smooth[vals].astype(np.float64)
    ng = len(g)
    Sv = F.mean(axis=0)
    Pair = (F.T @ F) / ng
    base = np.empty(vals.shape, dtype=np.int64)
    nbase = np.empty(npts, dtype=np.int64)
    for k in range(npts):
        base[:, k], nbase[k] = compact(base_cells(vals[:, k], h))
    lr = local_ratios(smooth, int(np.median(vals)), y) if prim_model else None
    m = Sv.sum()
    tot = dict(ex=0.0, cov=0.0, B=0.0, R=0.0, absR=0.0, thin=0.0, prim=0.0)
    rr, cc = [], []
    norms = np.sqrt((ab.astype(float) ** 2).sum(axis=1))
    for a in range(npts):
        for b in range(a + 1, npts):
            D = abs(int(ab[a, 0] * ab[b, 1] - ab[a, 1] * ab[b, 0]))
            mids = factor_mid(D, y)
            c1 = base[:, a]
            c2 = base[:, b]
            k1, k2 = int(nbase[a]), int(nbase[b])
            for p, k in mids:
                c1 = c1 * (k + 1) + capped_v(vals[:, a], p, k)
                c2 = c2 * (k + 1) + capped_v(vals[:, b], p, k)
                k1 *= k + 1
                k2 *= k + 1
            p1 = cell_mean(c1, F[:, a], k1)
            p2 = cell_mean(c2, F[:, b], k2)
            ss = Sv[a] * Sv[b]
            cov_nu = float((p1 * p2).mean() - ss)
            V1 = float((p1 * p1).mean() - Sv[a] ** 2)
            V2 = float((p2 * p2).mean() - Sv[b] ** 2)
            ex = float(Pair[a, b] - ss)
            R = ex - cov_nu
            w = 2.0  # ordered pairs
            tot["ex"] += w * ex
            tot["cov"] += w * cov_nu
            tot["B"] += w * np.sqrt(max(V1, 0) * max(V2, 0))
            tot["R"] += w * R
            tot["absR"] += w * abs(R)
            if norms[a] * norms[b] / D > Kthin:
                tot["thin"] += w * ex
            if lr is not None:
                mid_ps = {p for p, _ in mids}
                lp = sum(v for p, v in lr.items() if p not in mid_ps)
                tot["prim"] += w * ss * (np.exp(lp) - 1.0)
            if ss > 0:
                rr.append(R / ss)
                cc.append(ex / ss)
    m2 = m * m
    return dict(q=q, r=r, npts=npts, ng=ng, m=float(m), phi_ex=tot["ex"] / m2, C_nu=tot["cov"] / m2,
                B1=tot["B"] / m2, R_lab=tot["R"] / m2, absR=tot["absR"] / m2, thin_ex=tot["thin"] / m2,
                R_prim=(tot["prim"] / m2 if lr is not None else None),
                sd_pair_R=float(np.std(rr)), sd_pair_cov=float(np.std(cc)), n_pairs=len(rr))
