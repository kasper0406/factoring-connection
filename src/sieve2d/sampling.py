"""HARNESS: seeded uniform samplers for S_n (prime pairs) and for the integer
region R_n = {x < y <= 2x, 2^(n-1) <= xy < 2^n} at sizes where exhaustion is
impossible.  Uniformity: independent uniform primes in the p- and q-boxes,
then rejection to the region, is exactly uniform over S_n (the region is a
subset of the product box)."""
from __future__ import annotations

import gmpy2
import numpy as np

from .coverage import sn_bounds

_SMALL = np.array([3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47], dtype=np.int64)


def uniform_primes(lo: int, hi: int, count: int, rng: np.random.Generator) -> np.ndarray:
    out = []
    got = 0
    while got < count:
        x = rng.integers(lo, hi + 1, size=max(4 * count, 4096), dtype=np.int64)
        x = x[(x & 1) == 1]
        keep = np.ones(x.shape, dtype=bool)
        for s in _SMALL:
            keep &= (x % s) != 0
        # numbers removed by the prefilter are composite (x > 47), so this is exact
        cand = x[keep]
        ok = np.fromiter((gmpy2.is_prime(int(v), 40) for v in cand.tolist()), dtype=bool,
                         count=cand.size)
        pr = cand[ok]
        out.append(pr)
        got += pr.size
    return np.concatenate(out)[:count]


def in_region(n: int, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    lo, hi = 1 << (n - 1), 1 << n
    # products may exceed int64 for n > 62: use Python ints via object arrays then
    if n <= 60:
        prod = x * y
        return (x < y) & (y <= 2 * x) & (prod >= lo) & (prod < hi)
    xo = x.astype(object)
    yo = y.astype(object)
    prod = xo * yo
    return (x < y) & (y <= 2 * x) & np.array([(lo <= v < hi) for v in prod], dtype=bool)


def sample_prime_pairs(n: int, count: int, seed: int):
    assert n <= 122
    rng = np.random.default_rng([0x51E7E2D, n, seed])
    plo, phi, qlo, qhi = sn_bounds(n)
    ps, qs = [], []
    got = 0
    while got < count:
        m = count
        p = uniform_primes(plo, phi, m, rng)
        q = uniform_primes(qlo, qhi, m, rng)
        ok = in_region(n, p, q)
        ps.append(p[ok])
        qs.append(q[ok])
        got += int(ok.sum())
    return np.concatenate(ps)[:count], np.concatenate(qs)[:count]


def sample_integer_pairs(n: int, count: int, seed: int):
    assert n <= 122
    rng = np.random.default_rng([0x1A7E6E5, n, seed])
    plo, phi, qlo, qhi = sn_bounds(n)
    xs, ys = [], []
    got = 0
    while got < count:
        x = rng.integers(plo, phi + 1, size=2 * count, dtype=np.int64)
        y = rng.integers(qlo, qhi + 1, size=2 * count, dtype=np.int64)
        ok = in_region(n, x, y)
        xs.append(x[ok])
        ys.append(y[ok])
        got += int(ok.sum())
    return np.concatenate(xs)[:count], np.concatenate(ys)[:count]
