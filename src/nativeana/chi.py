"""Kronecker-symbol tables and estimators of L(1, chi_D).  OPERATIONAL: sees D (or N, m) only.

Floating point is used only for L-value approximations (allowed for fits/plots by rule 1); every character value is
exact (gmpy2.kronecker + complete multiplicativity).
"""
from __future__ import annotations

import math
import random
from typing import Dict, List, Sequence

import gmpy2
import numpy as np


def primes_upto(x: int) -> np.ndarray:
    """All primes <= x (simple numpy sieve)."""
    if x < 2:
        return np.zeros(0, dtype=np.int64)
    s = np.ones(x + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(x ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)


def chi_table(D: int, T: int) -> np.ndarray:
    """chi_D(n) for 0 <= n <= T as int8, via chi_D at primes and complete multiplicativity (exact)."""
    chi = np.ones(T + 1, dtype=np.int8)
    chi[0] = 0
    for ell in primes_upto(T).tolist():
        v = gmpy2.kronecker(D, ell)
        if v == 1:
            continue
        if v == 0:
            chi[ell::ell] = 0
            continue
        pe = ell
        while pe <= T:
            chi[pe::pe] *= -1
            pe *= ell
    return chi


def kron(D: int, n: int) -> int:
    return int(gmpy2.kronecker(D, n))


# ---------------------------------------------------------------- estimators of h(D), D < -4 fundamental
def euler_product_h(D: int, chi: np.ndarray, xs: Sequence[int]) -> Dict[int, float]:
    """h_hat = sqrt|D|/pi * prod_{l <= x} (1 - chi(l)/l)^{-1}, for each x in xs."""
    Q = abs(D)
    pr = primes_upto(max(xs))
    terms = -np.log1p(-chi[pr].astype(np.float64) / pr)
    cum = np.cumsum(terms)
    out = {}
    for x in xs:
        k = int(np.searchsorted(pr, x, side="right"))
        out[x] = math.sqrt(Q) / math.pi * math.exp(cum[k - 1] if k else 0.0)
    return out


def partial_sum_h(D: int, chi: np.ndarray, Ts: Sequence[int]) -> Dict[int, float]:
    Q = abs(D)
    n = np.arange(1, max(Ts) + 1, dtype=np.float64)
    cum = np.cumsum(chi[1:max(Ts) + 1] / n)
    return {T: math.sqrt(Q) / math.pi * float(cum[T - 1]) for T in Ts}


def smooth_sum_h(D: int, chi: np.ndarray, Ts: Sequence[int]) -> Dict[int, float]:
    """sum_n chi(n) exp(-(n/T)^2)/n, n <= 6T."""
    Q = abs(D)
    out = {}
    for T in Ts:
        n = np.arange(1, 6 * T + 1, dtype=np.float64)
        out[T] = math.sqrt(Q) / math.pi * float(np.sum(chi[1:6 * T + 1] * np.exp(-(n / T) ** 2) / n))
    return out


def afe_h(D: int, chi: np.ndarray, Ts: Sequence[int], chunk: int = 1 << 20) -> Dict[int, float]:
    """Exact formula h = sum_n chi(n)[erfc(n sqrt(pi/Q)) + sqrt(Q)/(pi n) exp(-pi n^2/Q)] (D < -4 fundamental),
    truncated at n <= T.  Terms beyond len(chi)-1 are dropped (caller sizes chi so that they are negligible)."""
    Q = abs(D)
    a = math.sqrt(math.pi / Q)
    b = math.sqrt(Q) / math.pi
    Tmax = min(max(Ts), len(chi) - 1)
    targets = sorted(set(min(T, Tmax) for T in Ts))
    out_by = {}
    total = 0.0
    ti = 0
    erfc = np.frompyfunc(math.erfc, 1, 1)
    for start in range(1, Tmax + 1, chunk):
        stop = min(start + chunk, Tmax + 1)
        n = np.arange(start, stop, dtype=np.float64)
        c = chi[start:stop].astype(np.float64)
        term = c * (erfc(n * a).astype(np.float64) + b / n * np.exp(-math.pi * n * n / Q))
        cum = np.cumsum(term) + total
        while ti < len(targets) and targets[ti] < stop:
            out_by[targets[ti]] = float(cum[targets[ti] - start])
            ti += 1
        total = float(cum[-1])
    return {T: out_by[min(T, Tmax)] for T in Ts}


def random_summation_h(D: int, Ks: Sequence[int], seed_key: str, cutoff_factor: int = 1 << 12) -> Dict[str, object]:
    """Srinivasan-type unbiased estimator of sum_{n < 2^B} chi(n)/n, 2^B >= cutoff_factor*|D|.
    One stream of max(Ks) samples; estimates use prefixes."""
    Q = abs(D)
    B = (cutoff_factor * Q).bit_length()
    rng = random.Random(seed_key)
    s = 0.0
    out = {}
    Ks = sorted(Ks)
    ki = 0
    for i in range(1, Ks[-1] + 1):
        j = rng.randrange(B)
        n = (1 << j) + rng.randrange(1 << j)
        s += B * (1 << j) * kron(D, n) / n
        if i == Ks[ki]:
            out[i] = math.sqrt(Q) / math.pi * s / i
            ki += 1
    return {"B": B, "est": out}
