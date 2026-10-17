"""Transparent baselines: trial division and Pollard rho (Brent), instrumented.

Operational: input N (and a seed); no factor information.  Counters use the same
keys as ring.py: ``rn_mul`` (multiplications mod N), ``rn_add``, ``gcd``, plus
``trial_div`` for trial divisions.
Each returns a proper divisor d (1 < d < N) or None.
"""
from __future__ import annotations

import random
from math import gcd, isqrt
from typing import Optional

from .metrics import Metrics, NullMetrics


def trial_division(N: int, m: Optional[Metrics] = None, bound: Optional[int] = None) -> Optional[int]:
    m = m if m is not None else NullMetrics()
    lim = isqrt(N) if bound is None else min(bound, isqrt(N))
    m.inc("trial_div")
    if N % 2 == 0 and N > 2:
        return 2
    d = 3
    while d <= lim:
        m.inc("trial_div")
        if N % d == 0:
            return d
        d += 2
    return None


def pollard_rho_brent(N: int, seed: int = 0, m: Optional[Metrics] = None,
                      max_iters: int = 1 << 26, batch: int = 128) -> Optional[int]:
    """Brent's variant of Pollard rho with batched gcds, f(x)=x^2+c.

    Restarts with a new c if the cycle closes with gcd == N.  Returns a proper
    divisor or None after `max_iters` total iterations.
    """
    m = m if m is not None else NullMetrics()
    if N % 2 == 0:
        return 2 if N > 2 else None
    rng = random.Random(f"rho|{N}|{seed}")
    iters = 0
    while iters < max_iters:
        y = rng.randrange(1, N)
        c = rng.randrange(1, N)
        g = r = q = 1
        x = ys = y
        while g == 1 and iters < max_iters:
            x = y
            for _ in range(r):
                y = (y * y + c) % N
            m.inc("rn_mul", r)
            m.inc("rn_add", r)
            iters += r
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(batch, r - k)):
                    y = (y * y + c) % N
                    q = q * abs(x - y) % N
                m.inc("rn_mul", 2 * min(batch, r - k))
                m.inc("rn_add", 2 * min(batch, r - k))
                iters += min(batch, r - k)
                g = gcd(q, N)
                m.inc("gcd")
                k += batch
            r *= 2
        if g == N:
            g = 1
            while g == 1:
                ys = (ys * ys + c) % N
                m.inc("rn_mul")
                g = gcd(abs(x - ys), N)
                m.inc("gcd")
        if 1 < g < N:
            m.set("rho_iters", iters)
            return g
        m.inc("rho_restarts")
    return None
