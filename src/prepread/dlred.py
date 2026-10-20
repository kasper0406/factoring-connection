"""Prop O1-DL (operational part): discrete logs from an observable correlated with a character.

Given p, a base g, a target y = g^x, query access to f: F_p^* -> unit disk (vectorised: f(array)),
a frequency j and the coefficient Fj = E_t f(g^t) e(-jt/(p-1)) (or an estimate of it), recover x.

Steps (results note, Prop O1-DL):
  1. est(z) = mean_s f(z g^s) e(-j s/(p-1)) / Fj  ~  e(j log_g(z)/(p-1)), m samples per estimate.
  2. theta_k = arg est(y^(2^k)) / 2pi ~ frac(2^k theta), theta = j x/(p-1), k = 0..K.
  3. Backward bootstrapping: theta_K, then theta_k in {theta_{k+1}/2, (theta_{k+1}+1)/2}, pick the
     candidate closest to the measurement; after K steps theta is known to < 1/(2(p-1)).
  4. j x mod (p-1) = round(theta (p-1)); x mod (p-1)/gcd(j, p-1); brute-force the cofactor.
The procedure never computes a discrete log itself; it uses only f-queries, products and powers.
"""
from __future__ import annotations
import math
import numpy as np


def _circ(a, b):
    d = (a - b) % 1.0
    return min(d, 1.0 - d)


def recover_dlog(p, g, y, f, j, Fj, rng, m=None, max_cofactor=10 ** 5, gpow=None, chunk=1 << 18):
    """gpow: optional table gpow[s] = g^s mod p (a speed shortcut; it is just repeated
    multiplication by g and leaks nothing)."""
    n = p - 1
    if m is None:
        m = int(math.ceil(200.0 / max(abs(Fj), 1e-9) ** 2))
    K = int(math.ceil(math.log2(n))) + 3
    queries = 0
    thetas = []
    z = y % p
    for k in range(K + 1):
        acc = 0j
        left = m
        while left > 0:
            mm = min(left, chunk)
            s = rng.integers(0, n, size=mm)
            if gpow is not None:
                gs = gpow[s]
            else:
                gs = np.array([pow(g, int(t), p) for t in s], dtype=np.int64)
            vals = f((z * gs) % p)
            acc += np.sum(vals * np.exp(-2j * np.pi * ((j * s) % n) / n))
            left -= mm
        queries += m
        est = acc / m / Fj
        thetas.append((np.angle(est) / (2 * np.pi)) % 1.0)
        z = z * z % p
    th = thetas[K]
    for k in range(K - 1, -1, -1):
        c1, c2 = th / 2.0, (th + 1.0) / 2.0
        th = c1 if _circ(c1, thetas[k]) <= _circ(c2, thetas[k]) else c2
    jx = int(round(th * n)) % n
    gg = math.gcd(j, n)
    if jx % gg:
        return None, queries
    n2 = n // gg
    x0 = (jx // gg) * pow(j // gg, -1, n2) % n2 if n2 > 1 else 0
    if gg > max_cofactor:
        return None, queries
    for t in range(gg):
        x = x0 + t * n2
        if pow(g, x, p) == y % p:
            return x, queries
    return None, queries
