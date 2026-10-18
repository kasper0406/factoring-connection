"""HARNESS / DIAGNOSTIC code for sieve2d: coverage of (p, q) pairs by the
Theorem D' targets (l | Bp - Aq), computed from p and q.  Never imported by
operational code.

Notation (notes/sieve2d_results.md):
  V      weight set of D' (coprime (A,B) in [1,H]^2, A != B)
  band   a list of primes l
  M(p,q) = #{(l, v) in band x V : l | Bp - Aq}       (cover multiplicity)
  uncovered  <=>  M = 0.

For a prime l, (p, q) mod l determines M_l := #{v : l | Bp - Aq}:
  (p,q) == (0,0)            -> |V|
  q != 0, r = p/q mod l     -> C_l[r] = #{v : B r == A (mod l)}
  q == 0, p != 0            -> C_l[l] = #{v : B == 0 (mod l)}   (point at infinity)
"""
from __future__ import annotations

from dataclasses import dataclass
from math import exp, log
from typing import Dict, List, Sequence, Tuple

import numpy as np


@dataclass
class BandTables:
    l: int
    C: np.ndarray        # int64, size l+1: projective-line class counts (index l = infinity)
    inv: np.ndarray      # int64, size l: inverse table (inv[0] = 0 unused)
    nV: int


def modpow_arr(x: np.ndarray, e: int, l: int) -> np.ndarray:
    """x^e mod l elementwise (l < 2^31)."""
    res = np.ones(x.shape, dtype=np.int64)
    b = x % l
    while e:
        if e & 1:
            res = (res * b) % l
        b = (b * b) % l
        e >>= 1
    return res


def inverse_table(l: int) -> np.ndarray:
    """inv[i] = i^(-1) mod l for prime l (inv[0] = 0)."""
    return modpow_arr(np.arange(l, dtype=np.int64), l - 2, l)


def class_table(l: int, V: Sequence[Tuple[int, int]]) -> np.ndarray:
    C = np.zeros(l + 1, dtype=np.int64)
    for (A, B) in V:
        b = B % l
        if b == 0:
            C[l] += 1
        else:
            C[(A * pow(b, -1, l)) % l] += 1
    return C


def band_tables(band: Sequence[int], V: Sequence[Tuple[int, int]]) -> List[BandTables]:
    out = []
    for l in band:
        C = np.zeros(l + 1, dtype=np.int64)
        for (A, B) in V:
            b = B % l
            if b == 0:
                C[l] += 1          # Bp - Aq == -Aq: covers q == 0 (A != 0 mod l since gcd)
            else:
                C[(A * pow(b, -1, l)) % l] += 1
        inv = inverse_table(l)
        out.append(BandTables(l=l, C=C, inv=inv, nV=len(V)))
    return out


def multiplicity(p: np.ndarray, q: np.ndarray, tabs: List[BandTables]) -> np.ndarray:
    """M(p,q) for int64 arrays p, q (values < 2^62)."""
    M = np.zeros(p.shape, dtype=np.int64)
    for t in tabs:
        l = t.l
        pm = p % l
        qm = q % l
        qz = qm == 0
        r = np.where(qz, l, (pm * t.inv[qm]) % l)
        c = t.C[r]
        both = qz & (pm == 0)
        if both.any():
            c = np.where(both, t.nV, c)
        M += c
    return M


def covered_by_count(p: np.ndarray, q: np.ndarray, tabs: List[BandTables]) -> np.ndarray:
    """Number of band primes l with M_l > 0 (the 'N' of the note)."""
    Nn = np.zeros(p.shape, dtype=np.int64)
    for t in tabs:
        l = t.l
        pm = p % l
        qm = q % l
        qz = qm == 0
        r = np.where(qz, l, (pm * t.inv[qm]) % l)
        c = (t.C[r] > 0) | (qz & (pm == 0))
        Nn += c
    return Nn


# ---------------------------------------------------------------------------
# Independent-residue models
# ---------------------------------------------------------------------------
def model_prime(tabs: List[BandTables]) -> Dict:
    """Model for prime pairs (p, q not divisible by any band prime): p/q mod l uniform
    on (Z/l)^*, independent over l.  Returns the uncovered probability, lambda's and
    the exact mean / variance / 4th central moment of M under the model."""
    logf = 0.0
    lam_eff = 0.0
    mu = var = m4c = m3 = 0.0
    for t in tabs:
        l = t.l
        Cs = t.C[1:l].astype(float)          # nonzero residues r = 1..l-1
        rho = float((Cs > 0).sum())
        pc = rho / (l - 1)
        lam_eff += pc
        logf += log(1.0 - pc) if pc < 1 else float("-inf")
        m1 = Cs.mean()
        d = Cs - m1
        mu += m1
        v = (d ** 2).mean()
        var += v
        m3 += (d ** 3).mean()
        m4c += (d ** 4).mean() - 3 * v * v
    m4 = m4c + 3 * var * var
    return {"f_model": exp(logf), "lam_eff": lam_eff, "mu": mu, "var": var, "m3": m3, "m4": m4}


def model_integer(tabs: List[BandTables]) -> Dict:
    """Model for integer pairs: (p, q) uniform on (Z/l)^2, independent over l."""
    logf = 0.0
    mu = var = m4c = m3 = 0.0
    for t in tabs:
        l = t.l
        # values over the l^2 points: (0,0) -> nV; each projective class r (l+1 of them)
        # carries l-1 points with value C[r]
        vals = np.concatenate([[t.nV], t.C.astype(float)])
        wts = np.concatenate([[1.0], np.full(l + 1, l - 1.0)]) / (l * l)
        pc = float(wts[vals > 0].sum())
        logf += log(1.0 - pc) if pc < 1 else float("-inf")
        m1 = float((wts * vals).sum())
        d = vals - m1
        v = float((wts * d ** 2).sum())
        mu += m1
        var += v
        m3 += float((wts * d ** 3).sum())
        m4c += float((wts * d ** 4).sum()) - 3 * v * v
    return {"f_model": exp(logf), "mu": mu, "var": var, "m3": m3, "m4": m4c + 3 * var * var}


def lam_nominal(band: Sequence[int], nV: int) -> float:
    return nV * sum(1.0 / l for l in band)


# ---------------------------------------------------------------------------
# Populations
# ---------------------------------------------------------------------------
def primes_in(lo: int, hi: int) -> np.ndarray:
    """Primes in [lo, hi] (hi <= ~1e8), numpy sieve."""
    s = np.ones(hi + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(hi ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    idx = np.nonzero(s)[0]
    return idx[idx >= lo].astype(np.int64)


def sn_bounds(n: int):
    from math import isqrt
    plo = isqrt(1 << (n - 2))
    phi = isqrt((1 << n) - 1) + 1
    qlo = isqrt(1 << (n - 1))
    qhi = isqrt(1 << (n + 1)) + 1
    return plo, phi, qlo, qhi


def exhaustive_prime_pairs(n: int) -> Tuple[np.ndarray, np.ndarray]:
    """All (p, q) in S_n = {p < q <= 2p primes, 2^(n-1) <= pq < 2^n}."""
    plo, phi, qlo, qhi = sn_bounds(n)
    pr = primes_in(2, max(phi, qhi) + 2)
    lo, hi = 1 << (n - 1), 1 << n
    ps, qs = [], []
    P = pr[(pr >= plo) & (pr <= phi)]
    for p in P.tolist():
        a = max(p + 1, -(-lo // p))
        b = min(2 * p, (hi - 1) // p)
        if a > b:
            continue
        i0 = np.searchsorted(pr, a, "left")
        i1 = np.searchsorted(pr, b, "right")
        if i1 > i0:
            qs.append(pr[i0:i1])
            ps.append(np.full(i1 - i0, p, dtype=np.int64))
    return np.concatenate(ps), np.concatenate(qs)


def exhaustive_integer_pairs(n: int) -> Tuple[np.ndarray, np.ndarray]:
    """All integer (x, y) with x < y <= 2x, 2^(n-1) <= xy < 2^n (n <= 28)."""
    plo, phi, qlo, qhi = sn_bounds(n)
    lo, hi = 1 << (n - 1), 1 << n
    xs, ys = [], []
    for x in range(plo, phi + 1):
        a = max(x + 1, -(-lo // x))
        b = min(2 * x, (hi - 1) // x)
        if a > b:
            continue
        ys.append(np.arange(a, b + 1, dtype=np.int64))
        xs.append(np.full(b - a + 1, x, dtype=np.int64))
    return np.concatenate(xs), np.concatenate(ys)


def degenerate_mask(x: np.ndarray, y: np.ndarray, V) -> np.ndarray:
    """Points with Bx == Ay for some v in V (never prime pairs)."""
    m = np.zeros(x.shape, dtype=bool)
    for (A, B) in V:
        m |= (B * x == A * y)
    return m
