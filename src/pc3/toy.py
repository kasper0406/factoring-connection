"""pc3 toy model: the RATIONAL block of Hyp PC3 (multidescent §3.4) at toy scale.

DIAGNOSTIC / HARNESS code (notes/pc3_prereg.md).  No N is involved: the rational block is a statement about
smooth values of the linear forms G_g(P) = g1*a - g0*b, g uniform on primitive vectors of [M/2, M]^2.
Smoothness is decided EXACTLY (all primes <= y and their powers are divided out of every n <= X).

Objects (mirroring multidescent §3.2 / this note §1):
  label q = (q, r): lattice L_q = {(a, b): a = r b mod q}, Gauss-reduced basis e1, e2;
  child box: primitive P' = i e1 + j e2, max(|i|,|j|) <= S, modulo +-1;
  parents of the label: primitive P = i0 e1 + j0 e2 with K <= |(i0, j0)|_inf < 2K (a random sample);
  a parent is accepted for g iff |G_g(P)| = l * s with l prime in (Y, 2Y] and s y-smooth (l = the special-q);
  a child point is good for g iff |G_g(P')| is y-smooth.
N(g) = number of good child points, Z(g) = number of accepted parents (the label's parent count).
"""
from __future__ import annotations

from math import gcd, isqrt
from typing import Dict, List, Tuple

import numpy as np


# ------------------------------------------------------------------ exact sieves
def primes_upto(n: int) -> np.ndarray:
    if n < 2:
        return np.zeros(0, dtype=np.int64)
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)


def build_tables(X: int, y: int, Y: int) -> Tuple[np.ndarray, np.ndarray]:
    """Exact tables for 0 <= n <= X:  smooth[n] = (n >= 1 and P+(n) <= y);
    parent[n] = (n = l*s, l prime in (Y, 2Y], s y-smooth).  Requires Y > y.  Peak memory ~ 6 X bytes."""
    assert Y > y
    rem = np.arange(X + 1, dtype=np.uint32)  # 4X bytes
    for p in primes_upto(y):
        p = int(p)
        pk = p
        while pk <= X:
            rem[pk::pk] //= p
            pk *= p
    smooth = rem == 1
    smooth[0] = False
    big = primes_upto(2 * Y)
    isbig = np.zeros(2 * Y + 1, dtype=bool)
    isbig[big[big > Y]] = True
    parent = np.zeros(X + 1, dtype=bool)
    step = 5_000_000
    for c0 in range(0, X + 1, step):
        blk = rem[c0:c0 + step]
        idx = np.nonzero(blk <= 2 * Y)[0]
        parent[c0 + idx] = isbig[blk[idx]]
    del rem
    return smooth, parent


# ------------------------------------------------------------------ lattices and boxes
def gauss_reduce(u: Tuple[int, int], v: Tuple[int, int]) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]

    if n2(u) > n2(v):
        u, v = v, u
    while True:
        k = (u[0] * v[0] + u[1] * v[1])
        m = round(k / n2(u))
        v = (v[0] - m * u[0], v[1] - m * u[1])
        if n2(v) >= n2(u):
            return u, v
        u, v = v, u


def label_basis(q: int, r: int):
    return gauss_reduce((q, 0), (r % q, 1))


def canon(P: Tuple[int, int]) -> Tuple[int, int]:
    a, b = P
    if b < 0 or (b == 0 and a < 0):
        return (-a, -b)
    return (a, b)


def child_box(e1, e2, S: int) -> List[Tuple[int, int, int, int]]:
    """Primitive points P' = i e1 + j e2, max(|i|,|j|) <= S, gcd(i,j)=1, P' primitive in Z^2, one per +-pair.
    Returns (i, j, a', b')."""
    out = {}
    for i in range(-S, S + 1):
        for j in range(-S, S + 1):
            if (i, j) == (0, 0) or gcd(i, j) != 1:
                continue
            a = i * e1[0] + j * e2[0]
            b = i * e1[1] + j * e2[1]
            if gcd(a, b) != 1:
                continue
            key = canon((a, b))
            if key not in out:
                out[key] = (i, j, key[0], key[1])
    return list(out.values())


def sample_parents(e1, e2, K: int, n: int, rng: np.random.Generator) -> List[Tuple[int, int, int, int]]:
    out = {}
    tries = 0
    while len(out) < n and tries < 50 * n:
        tries += 1
        i0 = int(rng.integers(-2 * K + 1, 2 * K))
        j0 = int(rng.integers(-2 * K + 1, 2 * K))
        if max(abs(i0), abs(j0)) < K or gcd(i0, j0) != 1:
            continue
        a = i0 * e1[0] + j0 * e2[0]
        b = i0 * e1[1] + j0 * e2[1]
        if gcd(a, b) != 1:
            continue
        key = canon((a, b))
        out.setdefault(key, (i0, j0, key[0], key[1]))
    return list(out.values())


def sample_g(M: int, n: int, rng: np.random.Generator) -> np.ndarray:
    """n primitive g = (g0, g1) uniform on [M/2, M]^2 (rejection on gcd)."""
    lo = (M + 1) // 2
    gs = []
    need = n
    while need > 0:
        cand = rng.integers(lo, M + 1, size=(int(need * 1.7) + 16, 2))
        g0, g1 = cand[:, 0], cand[:, 1]
        ok = np.gcd(g0, g1) == 1
        cand = cand[ok][:need]
        gs.append(cand)
        need -= len(cand)
    return np.concatenate(gs).astype(np.int64)


def values(g: np.ndarray, pts: np.ndarray) -> np.ndarray:
    """|G_g(P)| = |g1 a - g0 b| for all (g, P): shape (len(g), len(pts)), int64."""
    v = np.outer(g[:, 1], pts[:, 0])
    v -= np.outer(g[:, 0], pts[:, 1])
    np.abs(v, out=v)
    return v


def indicator(table: np.ndarray, vals: np.ndarray) -> np.ndarray:
    if vals.max() >= len(table):
        raise ValueError("value exceeds table size: %d >= %d" % (vals.max(), len(table)))
    return table[vals]


# ------------------------------------------------------------------ local model (E1)
def vp(n: int, p: int) -> int:
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


def factor_small(n: int) -> Dict[int, int]:
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def tau_table(smooth: np.ndarray, lo: int, hi: int, y: int, kmax: int = 6) -> Dict[Tuple[int, int], float]:
    """tau_p(k) = Pr(n smooth | v_p(n) = k) / Pr(n smooth), n uniform in [lo, hi), p <= y, k < kmax
    (k = kmax means v_p >= kmax).  Exact counts."""
    seg = smooth[lo:hi]
    n = np.arange(lo, hi, dtype=np.int64)
    base = seg.mean()
    out = {}
    for p in primes_upto(y):
        p = int(p)
        v = np.zeros(len(n), dtype=np.int64)
        m = n.copy()
        for _ in range(kmax):
            div = (m % p) == 0
            v += div
            m = np.where(div, m // p, m)
        for k in range(kmax + 1):
            sel = v == k
            if sel.any():
                out[(p, k)] = float(seg[sel].mean() / base)
    return out


def local_model(D: int, tau: Dict[Tuple[int, int], float], y: int, kmax: int = 6) -> float:
    """C_loc(D) = E[tau(gcd(n,D))^2]/E[tau(gcd(n,D))]^2 - 1 for n uniform, tau multiplicative over p:
    on the pair lattice gcd(n1, D) = gcd(n2, D) (n1 = lam n2 mod D, lam a unit), so the principal joint is
    E[tau1 tau2] with tau1 = tau2 evaluated at the common gcd (this note §2.3).  Primes p > y dividing D:
    tau = 0 on p | n (smooth numbers avoid p), 1/(1-1/p) otherwise (normalised)."""
    fac = factor_small(D)
    ratio = 1.0
    for p, e in fac.items():
        if p <= y:
            probs, taus = [], []
            for k in range(0, e + 1):
                pk = (1 - 1 / p) * p ** (-k) if k < e else p ** (-e)
                if k < e:
                    t = tau.get((p, min(k, kmax)), 1.0)
                else:
                    # v_p(n) >= e: average tau over k >= e with geometric weights (tail lumped at kmax)
                    if e >= kmax:
                        t = tau.get((p, kmax), 0.0)
                    else:
                        num, den = 0.0, 0.0
                        for kk in range(e, kmax + 1):
                            w = (1 - 1 / p) * p ** (-(kk - e)) if kk < kmax else p ** (-(kmax - e))
                            num += w * tau.get((p, kk), 0.0)
                            den += w
                        t = num / den
                probs.append(pk)
                taus.append(t)
            probs = np.array(probs)
            taus = np.array(taus)
            ratio *= float((probs * taus ** 2).sum() / (probs * taus).sum() ** 2)
        else:
            # P(p | n) = 1/p, tau = 0 there
            ratio *= 1.0 / (1 - 1 / p)
    return ratio - 1.0


# ------------------------------------------------------------------ principal proxy (E1, after the smoke test)
SMALL = ((2, 6), (3, 4), (5, 3), (7, 2), (11, 1), (13, 1))  # (p, cap)


def capped_vals(n: np.ndarray, p: int, cap: int) -> np.ndarray:
    v = np.zeros(n.shape, dtype=np.int64)
    m = n.copy()
    for _ in range(cap):
        div = (m % p) == 0
        v += div
        m = np.where(div, m // p, m)
    return v


def class_id(n: np.ndarray) -> np.ndarray:
    cid = np.zeros(n.shape, dtype=np.int64)
    for p, cap in SMALL:
        cid = cid * (cap + 1) + capped_vals(n, p, cap)
    return cid


def n_classes() -> int:
    k = 1
    for _, cap in SMALL:
        k *= cap + 1
    return k


def principal_tables(smooth: np.ndarray, lo: int, hi: int, y: int, kmax: int = 3):
    """Exact empirical local factors on n in [lo, hi): tau_cls[c] = Pr(smooth | class c)/Pr(smooth) for the
    small-prime valuation class (p <= 13, capped), and tau_p(k) for 17 <= p <= y (k capped at kmax)."""
    n = np.arange(lo, hi, dtype=np.int64)
    sm = smooth[lo:hi]
    base = sm.mean()
    cid = class_id(n)
    cnt = np.bincount(cid, minlength=n_classes()).astype(float)
    hit = np.bincount(cid, weights=sm.astype(float), minlength=n_classes())
    tau_cls = np.where(cnt > 0, hit / np.maximum(cnt, 1) / base, 1.0)
    taup = {}
    for p in primes_upto(y):
        p = int(p)
        if p <= 13:
            continue
        v = capped_vals(n, p, kmax)
        for k in range(kmax + 1):
            sel = v == k
            taup[(p, k)] = float(sm[sel].mean() / base) if sel.any() else 1.0
    return tau_cls, taup


def principal_proxy(vals: np.ndarray, tau_cls, taup, kmax: int = 3) -> np.ndarray:
    """f(n) = tau_cls(class(n)) * prod_{17<=p<=y} tau_p(v_p(n)); the archimedean factor is applied separately."""
    f = tau_cls[class_id(vals)]
    for (p, k), t in taup.items():
        if k == 0:
            v = capped_vals(vals, p, kmax)
            arr = np.array([taup[(p, kk)] for kk in range(kmax + 1)])
            f = f * arr[v]
    return f
