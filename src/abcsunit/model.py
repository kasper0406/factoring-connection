"""Counting models (floating point is used ONLY to evaluate heuristic predictions; no operational use).

Local-factor model for doubly-smooth congruent pairs (the N-dependent analogue of the
Lagarias–Soundararajan singular series):

  pairs(N, X, y) ≈ (N/|H|) · Σ_{k=1}^{⌊X/N⌋} 𝔖_y(kN) · Σ_{u smooth, u ≤ X − kN} d(u + kN),

with d(t) the local density of y-smooth integers near t, and the shift singular series
𝔖_y(h) = Π_{ℓ ≤ y} σ_ℓ(v_ℓ(h)) · Π_{y < ℓ, ℓ ∤ N} σ_ℓ^{>y}(h), where for ℓ ≤ y the local weight of a
y-smooth integer n with v_ℓ(n) = j is  w_ℓ(j) = ℓ^{j(1−α)}(1 − ℓ^{−α})/(1 − ℓ^{−1})  (so that Haar-mean
of w_ℓ is 1; α = α(X, y) the Hildebrand–Tenenbaum saddle point), and
σ_ℓ(m) = E_{u Haar}[ w_ℓ(v_ℓ(u)) · w_ℓ(v_ℓ(u + h)) ] for v_ℓ(h) = m.
For ℓ > y (ℓ ∤ N): σ = 1 − 1/(ℓ − 1)^2 if ℓ ∤ h, and ℓ/(ℓ − 1) if ℓ | h.  The factors at the primes of N
combine with equidistribution inside H = <S> ⊂ (Z/N)^* into the prefactor N/|H|.
"""
from __future__ import annotations

import math
from functools import lru_cache
from typing import Sequence

import numpy as np

from src.abcsunit import core


def saddle_alpha(logX: float, primes: Sequence[int]) -> float:
    """Solve Σ_{p ≤ y} log p / (p^α − 1) = log X for α > 0 (bisection)."""
    lo, hi = 1e-6, 2.0

    def f(a):
        return sum(math.log(p) / (p ** a - 1.0) for p in primes) - logX

    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _w(l: int, j: int, a: float) -> float:
    C = (1.0 - l ** (-a)) / (1.0 - 1.0 / l)
    return C * l ** (j * (1.0 - a))


def sigma_small(l: int, m: int, a: float, jmax: int = 200, av: float = None) -> float:
    """σ_ℓ(m) for ℓ ≤ y and a shift h with v_ℓ(h) = m: E_{u Haar}[w_ℓ(v(u); a) · w_ℓ(v(u+h); av)],
    av defaulting to a (one saddle point for both members).  Series truncated at jmax."""
    if av is None:
        av = a
    tot = 0.0
    for j in range(0, jmax):
        aj = l ** (-j) * (1.0 - 1.0 / l)
        if j < m:
            tot += aj * _w(l, j, a) * _w(l, j, av)
        elif j > m:
            tot += aj * _w(l, j, a) * _w(l, m, av)
        else:
            # v(u) = m: v(u+h) = m w.p. (ℓ−2)/(ℓ−1); else m + t (t ≥ 1) w.p. (1/(ℓ−1))·ℓ^{−(t−1)}(1 − 1/ℓ)
            inner = (l - 2) / (l - 1) * _w(l, m, av)
            for t in range(1, jmax):
                inner += (1.0 / (l - 1)) * l ** (-(t - 1)) * (1.0 - 1.0 / l) * _w(l, m + t, av)
            tot += aj * _w(l, m, a) * inner
    return tot


def shift_series(k: int, primes: Sequence[int], a: float, N: int = 1, ylimit_tail: int = 100000,
                 av: float = None) -> float:
    """𝔖_y(kN)/(N/φ(N)) restricted to ℓ ∤ N: product over ℓ ≤ y of σ_ℓ(v_ℓ(k)) times the ℓ > y factors
    (truncated at ylimit_tail; the tail product is 1 − O(1/ylimit_tail))."""
    s = 1.0
    for l in primes:
        m = 0
        kk = k
        while kk % l == 0:
            kk //= l
            m += 1
        s *= sigma_small(l, m, a, av=av)
    y = primes[-1]
    for l in core.primes_upto(ylimit_tail):
        if l <= y or N % l == 0:
            continue
        if k % l == 0:
            s *= l / (l - 1.0)
        else:
            s *= 1.0 - 1.0 / (l - 1.0) ** 2
    return s


def local_density(vals: np.ndarray, t: np.ndarray, delta: float = 0.1, Xmax: int = None) -> np.ndarray:
    """d(t) = #smooth in a window around t / window length, from the exact sorted smooth list `vals`
    (complete up to Xmax).  Symmetric window [t(1−δ), t(1+δ)] when t(1+δ) ≤ Xmax, else the left window
    [t(1−2δ), t]."""
    if Xmax is None:
        Xmax = int(vals[-1])
    sym = t * (1 + delta) <= Xmax
    hi = np.where(sym, np.floor(t * (1 + delta)), np.floor(t)).astype(np.int64)
    lo = np.where(sym, np.ceil(t * (1 - delta)), np.ceil(t * (1 - 2 * delta))).astype(np.int64)
    c = np.searchsorted(vals, hi, side="right") - np.searchsorted(vals, lo, side="left")
    return c / (hi - lo + 1)


def pair_model(vals: np.ndarray, N: int, X: int, H: int, primes: Sequence[int], vals_ext: np.ndarray,
               with_series: bool = True) -> float:
    """Local-factor prediction for #{u < v ≤ X smooth, u ≡ v (mod N)}.  `vals_ext` must list the smooth
    numbers up to at least 1.1·X (for the density window)."""
    a = saddle_alpha(math.log(X), primes)
    kmax = X // N
    tot = 0.0
    X_ext = X if vals_ext is vals else int(vals_ext[-1])
    for k in range(1, kmax + 1):
        h = k * N
        us = vals[vals <= X - h].astype(np.float64)
        if us.size == 0:
            continue
        d = local_density(vals_ext, us + h, Xmax=X_ext)
        s = shift_series(k, primes, a, N) if with_series else 1.0
        tot += s * float(d.sum())
    return tot * N / H



def pair_model_fast(vals: np.ndarray, N: int, X: int, H: int, primes: Sequence[int], vals_ext: np.ndarray,
                    with_series: bool = True, usample: int = 20000, seed: int = 12345,
                    k_filter=None) -> float:
    """Registered model M1 (= pair_model with one saddle point α(X, y)), evaluated with a seeded
    subsample of at most `usample` values of u (weight Ψ/usample) and a factored ℓ > y tail."""
    a = saddle_alpha(math.log(X), primes)
    kmax = X // N
    y = primes[-1]
    X_ext = X if vals_ext is vals else int(vals_ext[-1])
    rng = np.random.default_rng(seed)
    if vals.size > usample:
        us_all = np.sort(rng.choice(vals, size=usample, replace=False)).astype(np.float64)
        wt = vals.size / usample
    else:
        us_all = vals.astype(np.float64)
        wt = 1.0
    tailp = [l for l in core.primes_upto(100000) if l > y and N % l != 0]
    T0 = 1.0
    for l in tailp:
        T0 *= 1.0 - 1.0 / (l - 1.0) ** 2
    small_cache = {}
    tot = 0.0
    for k in range(1, kmax + 1):
        if k_filter is not None and not k_filter(k):
            continue
        h = k * N
        us = us_all[us_all <= X - h]
        if us.size == 0:
            continue
        d = local_density(vals_ext, us + h, Xmax=X_ext)
        if with_series:
            sig = []
            kk = k
            s = 1.0
            for l in primes:
                m = 0
                while kk % l == 0:
                    kk //= l
                    m += 1
                key = (l, m)
                if key not in small_cache:
                    small_cache[key] = sigma_small(l, m, a)
                s *= small_cache[key]
            s *= T0
            for l in tailp:
                if l > kk:
                    break
                if kk % l == 0:
                    s *= (l / (l - 1.0)) / (1.0 - 1.0 / (l - 1.0) ** 2)
        else:
            s = 1.0
        tot += s * float(d.sum()) * wt
    return tot * N / H
