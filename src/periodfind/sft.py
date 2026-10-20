"""OPERATIONAL sparse Fourier transform on Z_Q, Q = 2^m (Kushilevitz-Mansour style, low-bit buckets).

Bucket (l, alpha) = {k in Z_Q : k = alpha mod 2^l}. Its weight
  W(l, alpha) = sum_{k = alpha (2^l)} |F^(k)|^2 = E_x |E_{u in H_l} F(x+u) e(-alpha u/Q)|^2,
with H_l = 2^{m-l} Z_Q (size 2^l). Unbiased estimator: Re F(x+u) conj(F(x+u')) e(-alpha(u-u')/Q) for
uniform x, u, u' -- two queries per sample, bounded by ||F||_inf^2 <= 1.

The algorithm sees only an oracle F (query access), Q, tau, a seed; it never sees r. The order is
recovered from an output frequency k by continued fractions of k/Q and checked with a^{r'} = 1.
"""
from __future__ import annotations

import math
from typing import Callable, Dict, List

import numpy as np


class CountingOracle:
    def __init__(self, f: Callable[[np.ndarray], np.ndarray]):
        self.f = f
        self.queries = 0

    def __call__(self, x: np.ndarray) -> np.ndarray:
        self.queries += len(x)
        return self.f(x)


def hoeffding_samples(eps: float, delta: float) -> int:
    # real part in [-1, 1]: P(|mean - mu| > eps) <= 2 exp(-s eps^2 / 2)
    return int(math.ceil(2.0 * math.log(2.0 / delta) / eps ** 2))


def sft(oracle: CountingOracle, m: int, tau: float, seed: int, delta: float = 0.05,
        max_queries: int = 10 ** 8) -> Dict:
    """Return candidate tau-heavy frequencies (relative to the estimated ||F||^2)."""
    rng = np.random.default_rng(seed)
    mask = (np.uint64(1) << np.uint64(m)) - np.uint64(1)
    cap = int(math.ceil(4 / tau))
    s = hoeffding_samples(tau / 4, delta / (2 * m * cap))
    # estimate ||F||^2
    x = rng.integers(0, 2 ** 62, size=s, dtype=np.uint64) & mask
    norm2 = float(np.mean(np.abs(oracle(x)) ** 2))
    if norm2 == 0.0:
        return {"out": [], "norm2_est": 0.0, "samples_per_est": s, "queries": oracle.queries,
                "aborted": False}
    theta = tau * norm2
    survivors: List[int] = [0]
    for l in range(1, m + 1):
        shift = np.uint64(m - l)
        lowmask = (np.uint64(1) << np.uint64(l)) - np.uint64(1)
        nxt = []
        for alpha0 in survivors:
            for alpha in (alpha0, alpha0 + (1 << (l - 1))):
                if oracle.queries > max_queries:
                    return {"out": [], "norm2_est": norm2, "samples_per_est": s,
                            "queries": oracle.queries, "aborted": True}
                x = rng.integers(0, 2 ** 62, size=s, dtype=np.uint64) & mask
                v1 = rng.integers(0, 2 ** 62, size=s, dtype=np.uint64) & lowmask
                v2 = rng.integers(0, 2 ** 62, size=s, dtype=np.uint64) & lowmask
                f1 = oracle((x + (v1 << shift)) & mask)
                f2 = oracle((x + (v2 << shift)) & mask)
                dv = (v1 - v2) & lowmask  # exact mod 2^l (uint64 wraparound)
                with np.errstate(over="ignore"):
                    ph = (np.uint64(alpha) * dv) & lowmask
                phase = np.exp(-2j * np.pi * (ph.astype(np.float64) / float(2 ** l)))
                est = float(np.mean((f1 * np.conj(f2) * phase).real))
                if est >= theta / 2:
                    nxt.append((est, alpha))
        nxt.sort(reverse=True)
        survivors = [al for _, al in nxt[:cap]]
        if not survivors:
            break
    return {"out": survivors, "norm2_est": norm2, "samples_per_est": s,
            "queries": oracle.queries, "aborted": False}


def cf_denominators(k: int, Q: int, bound: int) -> List[int]:
    """Denominators of continued-fraction convergents of k/Q up to `bound`."""
    out = []
    a, b = k, Q
    h0, h1, k0, k1 = 0, 1, 1, 0
    while b:
        t = a // b
        a, b = b, a - t * b
        h0, h1 = h1, t * h1 + h0
        k0, k1 = k1, t * k1 + k0
        if k1 > bound:
            break
        out.append(k1)
    return out


def recover_order(freqs: List[int], Q: int, a: int, N: int) -> int:
    """Smallest r' (multiple-closure of CF denominators) with a^{r'} = 1, or 0."""
    cands = set()
    for k in freqs:
        if k == 0:
            continue
        for d in cf_denominators(k, Q, N):
            cands.add(d)
    best = 0
    cl = sorted(cands)
    # lcm of pairs too (one sample gives r/gcd(j, r))
    pool = set(cl)
    for i in range(len(cl)):
        for j in range(i + 1, len(cl)):
            pool.add(cl[i] * cl[j] // math.gcd(cl[i], cl[j]))
    for d in sorted(pool):
        if d > 1 and pow(a, d, N) == 1:
            return d
    return best
