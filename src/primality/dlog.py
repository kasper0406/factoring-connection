"""OPERATIONAL: recover s = p+q (or q-p) from a primality-test residue by an
interval discrete logarithm, then factor.  Inputs: N, a seed, public interval
conventions.  No factor, order or phi is used.  Cost ~ sqrt(interval width) = N^(1/4).

Prior art: Pollard 1978 (Math. Comp. 32, Sec. 3) proposed exactly this
(kangaroo on a^(2x) = a^(N+1), x = (p+q)/2); Shanks' BSGS is the deterministic
version; McKee-Pinch 1996 and Hittmeir 2018 refine it.  This module is a
baseline demonstration, not a new algorithm.
"""
from __future__ import annotations

import random
from math import gcd, isqrt
from typing import Callable, Optional, Tuple

import gmpy2

from src.common.metrics import Metrics, NullMetrics
from src.primality.ops import QuadRing, selfridge

_NULL = NullMetrics()


class Mult:
    """Abstract group: Z/N^* (ints) or the norm-one torus (pairs) with counted multiplications."""

    def __init__(self, N: int, D: Optional[int] = None, metrics=_NULL):
        self.N, self.D, self.m = N, D, metrics
        if D is not None:
            self.R = QuadRing(N, 0, -D, metrics)  # w^2 = D

    def one(self):
        return 1 if self.D is None else (1, 0)

    def mul(self, x, y):
        self.m.inc("grp_mul")
        if self.D is None:
            return x * y % self.N
        return self.R.mul(x, y)

    def pow(self, x, e: int):
        if e < 0:
            x, e = self.inv(x), -e
        r = self.one()
        while e:
            if e & 1:
                r = self.mul(r, x)
            e >>= 1
            if e:
                x = self.mul(x, x)
        return r

    def inv(self, x):
        if self.D is None:
            from src.common.ring import inv_mod
            self.m.inc("gcd")
            return inv_mod(x, self.N)  # raises NonUnit (a factor) if not invertible
        # norm-one elements: inverse = conjugate
        return (x[0], (-x[1]) % self.N)

    def key(self, x) -> int:
        return x if self.D is None else x[0] * self.N + x[1]


def bsgs(G: Mult, g, y, lo: int, hi: int) -> Optional[int]:
    """Least x in [lo, hi] with g^x = y (BSGS, ~2 sqrt(hi-lo) multiplications, sqrt memory)."""
    W = hi - lo + 1
    M = isqrt(W - 1) + 1
    table = {}
    cur = G.one()
    for j in range(M):
        table.setdefault(G.key(cur), j)
        cur = G.mul(cur, g)
    G.m.hi("stored", M)
    step = G.inv(G.pow(g, M))
    z = G.mul(y, G.inv(G.pow(g, lo)))  # y g^-lo = g^(x-lo)
    for i in range(M + 1):
        j = table.get(G.key(z))
        if j is not None:
            x = lo + i * M + j
            if x <= hi:
                return x
        z = G.mul(z, step)
    return None


def kangaroo(G: Mult, g, y, lo: int, hi: int, rng: random.Random, theta: Optional[int] = None,
             max_factor: int = 40) -> Optional[int]:
    """Pollard lambda with distinguished points (van Oorschot-Wiener, one tame, one wild).
    Returns x in [lo, hi] with g^x = y, or None after ~max_factor*sqrt(W) steps."""
    W = hi - lo
    sW = max(4, isqrt(W))
    k = 1
    while (2 ** k - 1) / k < sW / 2:
        k += 1
    jumps = [2 ** i for i in range(k)]
    salt = rng.getrandbits(32) | 1
    jpow = [G.pow(g, s) for s in jumps]
    if theta is None:
        theta = max(1, (sW.bit_length() - 4) // 2)
    mask = (1 << theta) - 1
    G.m.hi("jumps", k)

    def h(x):
        return ((G.key(x) * salt) >> 7) % k

    tame_pos = lo + W // 2
    tame = G.pow(g, tame_pos)
    wild, wild_pos = y, 0  # wild = g^(x + wild_pos)
    dp = {}
    steps, limit = 0, max_factor * sW + 64
    while steps < limit:
        for which in (0, 1):
            cur = tame if which == 0 else wild
            j = h(cur)
            if which == 0:
                tame = G.mul(tame, jpow[j]); tame_pos += jumps[j]; cur = tame; pos = tame_pos
            else:
                wild = G.mul(wild, jpow[j]); wild_pos += jumps[j]; cur = wild; pos = wild_pos
            steps += 1
            kk = G.key(cur)
            if (kk * salt >> 3) & mask == 0:
                prev = dp.get(kk)
                if prev is None:
                    dp[kk] = (which, pos)
                elif prev[0] != which:
                    tpos, wpos = (pos, prev[1]) if which == 0 else (prev[1], pos)
                    x = tpos - wpos
                    if lo <= x <= hi and G.pow(g, x) == y:
                        G.m.hi("stored", len(dp))
                        return x
    G.m.hi("stored", len(dp))
    return None


def factor_from_s(N: int, s: int) -> Optional[Tuple[int, int]]:
    disc = s * s - 4 * N
    if disc < 0 or not gmpy2.is_square(disc):
        return None
    r = isqrt(disc)
    p, q = (s - r) // 2, (s + r) // 2
    return (p, q) if p * q == N and 1 < p < N else None


def factor_from_diff(N: int, d: int) -> Optional[Tuple[int, int]]:
    disc = d * d + 4 * N
    if not gmpy2.is_square(disc):
        return None
    r = isqrt(disc)
    p, q = (r - abs(d)) // 2, (r + abs(d)) // 2
    return (p, q) if p * q == N and 1 < p < N else None


def fermat_interval(N: int) -> Tuple[int, int]:
    """s - 2 for balanced semiprimes with q/p <= ((3+sqrt5)/2)^2: s in [2 sqrt N, 3 sqrt N]."""
    return isqrt(4 * N) - 2, isqrt(9 * N) + 1 - 2


def factor_via_fermat_residue(N: int, seed: int, method: str = "kangaroo", metrics=None,
                              max_tries: int = 8) -> dict:
    """Compute the Fermat residue y = a^(N-1) (what the test outputs), recover x = s-2 from y alone by an
    interval DLP on [2 sqrt N - 2, 3 sqrt N - 2], then factor via s^2 - 4N = (q-p)^2."""
    m = metrics or Metrics()
    rng = random.Random(seed)
    G = Mult(N, None, m)
    lo, hi = fermat_interval(N)
    for attempt in range(max_tries):
        a = rng.randrange(2, N - 1)
        m.inc("gcd")
        d = gcd(a, N)
        if d != 1:
            return {"factors": (d, N // d), "class": "ACCIDENTAL", "attempts": attempt + 1, "metrics": m.to_record()}
        y = pow(a, N - 1, N)  # the Fermat test's residue
        m.inc("rn_mul", 2 * N.bit_length())
        x = bsgs(G, a, y, lo, hi) if method == "bsgs" else kangaroo(G, a, y, lo, hi, rng)
        if x is not None:
            f = factor_from_s(N, x + 2)
            if f:
                return {"factors": f, "class": "STRUCTURAL", "attempts": attempt + 1, "x": x,
                        "metrics": m.to_record()}
    return {"factors": None, "attempts": max_tries, "metrics": m.to_record()}


def factor_via_lucas_residue(N: int, seed: int, method: str = "kangaroo", metrics=None, max_tries: int = 8) -> dict:
    """Selfridge D ((D/N) = -1, so exactly one prime is inert): for a random norm-one z,
    z^(N+1) = z^(c) with c = +-(q-p) + 2 (Thm 5).  Interval |c - 2| <= sqrt(5N)."""
    m = metrics or Metrics()
    rng = random.Random(seed)
    sel = selfridge(N)
    if "D" not in sel or "P" not in sel:
        return {"factors": None, "selfridge": sel}
    D = sel["D"]
    G = Mult(N, D, m)
    W = isqrt(5 * N) + 1
    lo, hi = 2 - W, 2 + W
    for attempt in range(max_tries):
        a0, a1 = rng.randrange(N), rng.randrange(1, N)
        nu = (a0 * a0 - D * a1 * a1) % N
        m.inc("gcd")
        d = gcd(nu, N)
        if d != 1:
            if 1 < d < N:
                return {"factors": (d, N // d), "class": "ACCIDENTAL", "metrics": m.to_record()}
            continue
        inv = pow(nu, -1, N)
        a2 = G.R.mul((a0, a1), (a0, a1))
        z = (a2[0] * inv % N, a2[1] * inv % N)  # alpha / sigma(alpha), norm one
        y = G.pow(z, N + 1)  # Lucas-test residue
        x = bsgs(G, z, y, lo, hi) if method == "bsgs" else kangaroo(G, z, y, lo, hi, rng)
        if x is not None:
            f = factor_from_diff(N, x - 2)
            if f:
                return {"factors": f, "class": "STRUCTURAL", "attempts": attempt + 1, "x": x,
                        "metrics": m.to_record()}
    return {"factors": None, "attempts": max_tries, "metrics": m.to_record()}


def split_with_order_multiple(N: int, a: int, M: int, metrics=_NULL) -> Optional[int]:
    """Given M > 0 with a^M = 1 (mod N) (e.g. M = E - x from a collision a^x = a^E with x != E),
    walk the 2-power ladder a^(M'/2^i), M' = M / 2^v2(M) * 2^j, and gcd every value with N
    (Miller/Long splitting restricted to the single base a).  Returns a proper divisor or None."""
    if M <= 0:
        return None
    t = (M & -M).bit_length() - 1
    m = M >> t
    b = pow(a, m, N)
    for _ in range(t + 1):
        for z in (b - 1, b + 1):
            metrics.inc("gcd")
            d = gcd(z % N, N)
            if 1 < d < N:
                return d
        b = b * b % N
    return None
