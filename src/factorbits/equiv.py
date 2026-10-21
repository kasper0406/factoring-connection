"""OPERATIONAL reductions between prime-exponent addition and factoring (notes/factorbits_results.md §1).

A "vector" is a dict {prime: exponent}; value(v) = prod p**e.  The empty dict is 1.
Every reduction receives its oracle as an argument and counts the calls it makes.  The reductions
themselves only multiply, compare, take integer roots, run trial division up to an explicit bound
and run a probable-prime test (gmpy2.is_prime, BPSW; AKS in theory).  They never call a factoring
routine on N.  Oracles live in src/factorbits/harness.py (HARNESS).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, Optional

import gmpy2

Vec = Dict[int, int]


# ----------------------------------------------------------------------------- vector helpers
def value(v: Vec) -> int:
    out = 1
    for p, e in v.items():
        out *= p ** e
    return out


def vmul(a: Vec, b: Vec) -> Vec:
    """Multiplication = vector addition."""
    out = dict(a)
    for p, e in b.items():
        out[p] = out.get(p, 0) + e
        if out[p] == 0:
            del out[p]
    return out


def vdiv(a: Vec, b: Vec) -> Vec:
    out = dict(a)
    for p, e in b.items():
        out[p] = out.get(p, 0) - e
        if out[p] < 0:
            raise ValueError("not a divisor")
        if out[p] == 0:
            del out[p]
    return out


def is_prime(x: int) -> bool:
    return x >= 2 and bool(gmpy2.is_prime(x))


def trial_vec(x: int, bound: int) -> Optional[Vec]:
    """Full factorisation of x >= 1 by trial division up to `bound`, plus a primality test on the
    cofactor.  Returns None if the cofactor is composite (x not fully factored)."""
    if x < 1:
        raise ValueError
    v: Vec = {}
    d = 2
    while d <= bound and d * d <= x:
        while x % d == 0:
            v[d] = v.get(d, 0) + 1
            x //= d
        d += 1 if d == 2 else 2
    if x > 1:
        if x <= bound or is_prime(x) or (d * d > x):
            v[x] = v.get(x, 0) + 1
        else:
            return None
    return v


@dataclass
class Counter:
    calls: int = 0
    log: list = field(default_factory=list)


# ----------------------------------------------------------------------------- Prop E1: ladders
def factor_via_successor(N: int, succ: Callable[[Vec], Vec], ctr: Counter) -> Vec:
    """Horner ladder with free doubling: x <- 2x (vector op), x <- x+1 (one oracle call).
    Exactly popcount(N) - 1 calls.  The oracle is only ever queried on even inputs."""
    bits = bin(N)[2:]
    x: Vec = {}                                  # value 1 = leading bit
    for b in bits[1:]:
        x = vmul(x, {2: 1})
        if b == "1":
            ctr.calls += 1
            x = succ(x)
    assert value(x) == N
    return x


def factor_via_addition_window(N: int, add: Callable[[Vec, Vec], Vec], w: int, ctr: Counter) -> Vec:
    """Base-2^w Horner: x <- 2^w x (free), x <- x + d with d < 2^w factored by trial division.
    At most ceil(log2(N)/w) calls (one per nonzero digit after the leading one)."""
    digits = []
    M = N
    while M:
        digits.append(M % (1 << w))
        M >>= w
    digits.reverse()
    x = trial_vec(digits[0], 1 << w)
    assert x is not None
    for d in digits[1:]:
        x = vmul(x, {2: w})
        if d:
            vd = trial_vec(d, 1 << w)
            ctr.calls += 1
            x = add(x, vd)
    assert value(x) == N
    return x


def factor_via_addition_prime(N: int, add: Callable[[Vec, Vec], Vec], ctr: Counter,
                              max_gap: int = 10 ** 6) -> Optional[Vec]:
    """One call (heuristically): find the largest prime v < N with u = N - v fully factored by trial
    division (u <= max_gap), then ADD(u, v).  Prime gaps are O(log^2 N) under Cramer's model, so
    u is polynomially small and the search is polynomial heuristically (BHP: unconditionally
    u <= N^0.525 for the true gap)."""
    if N <= 3:
        return {N: 1} if N > 1 else {}
    v = N - 1
    while v > 1 and N - v <= max_gap:
        if is_prime(v):
            vu = trial_vec(N - v, max_gap)
            if vu is not None:
                ctr.calls += 1
                return add(vu, {v: 1})
        v -= 1
    return None


def factor_via_successor_multiplier(N: int, succ: Callable[[Vec], Vec], ctr: Counter,
                                    tmax: int = 10 ** 5) -> Optional[Vec]:
    """One successor call (heuristically): find t <= tmax with t*N - 1 prime (vector trivial),
    call succ on it to get the vector of t*N, divide out t (trial-factored)."""
    for t in range(1, tmax + 1):
        if is_prime(t * N - 1):
            vt = trial_vec(t, tmax)
            ctr.calls += 1
            vtn = succ({t * N - 1: 1})
            return vdiv(vtn, vt)
    return None


# ----------------------------------------------------------------------------- Prop E2: converse
def add_via_factoring(va: Vec, vb: Vec, fact: Callable[[int], Vec], ctr: Counter) -> Vec:
    """One complete-factorisation call (input size: log2(a+b) <= n, the convention of Prop E2)."""
    ctr.calls += 1
    return fact(value(va) + value(vb))


def add_via_splitting(va: Vec, vb: Vec, split: Callable[[int], int], ctr: Counter) -> Vec:
    """A splitting oracle returns a nontrivial divisor of a composite input.  Complete factorisation
    of S = a+b uses at most Omega(S) - 1 <= log2(S) - 1 calls (plus primality tests)."""
    def rec(x: int) -> Vec:
        if x == 1:
            return {}
        if is_prime(x):
            return {x: 1}
        ctr.calls += 1
        d = split(x)
        assert 1 < d < x and x % d == 0
        return vmul(rec(d), rec(x // d))
    return rec(value(va) + value(vb))


# ----------------------------------------------------------------------------- Prop E3: partial outputs
def factor_via_support(N: int, supp: Callable[[Vec, Vec], set], ctr: Counter) -> Vec:
    """Oracle returns only the set of primes dividing a+b.  Exponents by repeated division.
    Same ladder as factor_via_successor (popcount(N) - 1 calls)."""
    def full(x: int, primes: set) -> Vec:
        v: Vec = {}
        for p in primes:
            while x % p == 0:
                v[p] = v.get(p, 0) + 1
                x //= p
        assert x == 1
        return v
    bits = bin(N)[2:]
    x: Vec = {}
    for b in bits[1:]:
        x = vmul(x, {2: 1})
        if b == "1":
            ctr.calls += 1
            S = supp(x, {})
            x = full(value(x) + 1, S)
    return x


def iroot(x: int, m: int) -> Optional[int]:
    r, exact = gmpy2.iroot(x, m)
    return int(r) if exact else None


def factor_via_expmod(N: int, expmod: Callable[[Vec, Vec], Vec], m: int, ctr: Counter,
                      trial: int = 1000) -> Vec:
    """Oracle returns {p: e_p mod m} (zeros dropped) for a+b.  The primes with m | e_p are hidden;
    their product is an m-th power R = r^m.  r is trial-divided up to `trial`, tested for primality,
    and otherwise factored recursively by the same ladder (recursion depth <= log_m log2 N)."""
    def complete(x: int, part: Vec) -> Vec:
        v: Vec = {}
        for p in part:
            while x % p == 0:
                v[p] = v.get(p, 0) + 1
                x //= p
        if x == 1:
            return v
        r = iroot(x, m)
        assert r is not None, "hidden part must be an m-th power"
        vr = trial_vec(r, trial)
        if vr is None:
            vr = ladder(r)
        return vmul(v, {p: e * m for p, e in vr.items()})

    def ladder(M: int) -> Vec:
        bits = bin(M)[2:]
        x: Vec = {}
        for b in bits[1:]:
            x = vmul(x, {2: 1})
            if b == "1":
                ctr.calls += 1
                part = expmod(x, {})
                x = complete(value(x) + 1, part)
        return x
    return ladder(N)


# ----------------------------------------------------------------------------- Prop Z: Zech logarithms
def dlog_via_zech(h: int, p: int, g: int, zech: Callable[[int], int], ctr: Counter) -> int:
    """Discrete log in F_p^* from a Zech oracle Z(k) = log_g(1 + g^k) (k mod p-1; 1+g^k != 0).
    log 2 = Z(0) costs one call; Horner on the bits of h (1 <= h < p) costs popcount(h) - 1 more.
    Total <= floor(log2 p) + 1 calls; 1 + g^k never vanishes along the ladder (all values < p)."""
    assert 1 <= h < p
    if h == 1:
        return 0
    ctr.calls += 1
    L2 = zech(0)
    bits = bin(h)[2:]
    x = 0                                       # log of the current prefix (prefix 1)
    for b in bits[1:]:
        x = (x + L2) % (p - 1)
        if b == "1":
            ctr.calls += 1
            x = zech(x)
    return x


def zech_via_dlog(k: int, p: int, g: int, dlog: Callable[[int], int], ctr: Counter) -> int:
    """One DL call."""
    t = (1 + pow(g, k, p)) % p
    assert t != 0
    ctr.calls += 1
    return dlog(t)


# ----------------------------------------------------------------------------- randomized short-interval reductions
def _random_prime_below(x: int, theta: float, rng) -> int:
    w = max(2, 1 << int(theta * x.bit_length()))
    while True:
        v = rng.randrange(x - w, x)
        if is_prime(v):
            return v


def factor_via_addition_shortint(N: int, add: Callable[[Vec, Vec], Vec], ctr: Counter, rng,
                                 theta: float = 0.6, base: int = 1 << 20) -> Vec:
    """Las Vegas (with a certified primality test): X = u + v, v a random prime in [X - X^theta, X),
    u factored recursively; O(log log N) oracle calls (Huxley's short-interval PNT, theta > 7/12)."""
    if N <= base:
        return trial_vec(N, base)
    v = _random_prime_below(N, theta, rng)
    vu = factor_via_addition_shortint(N - v, add, ctr, rng, theta, base)
    ctr.calls += 1
    return add(vu, {v: 1})


def factor_via_expmod_shortint(N: int, expmod: Callable[[Vec, Vec], Vec], m: int, ctr: Counter, rng,
                               theta: float = 0.6, base: int = 1 << 20) -> Vec:
    """Codex's randomized polynomial reduction for the fixed-m exponent oracle: query (u, v) as above,
    strip the returned primes, the rest is r^m with r <= X^{1/m}; recurse on r.
    C(b) <= C(theta b) + C(b/m) + 1 oracle calls, i.e. polynomial."""
    if N <= base:
        return trial_vec(N, base)
    v = _random_prime_below(N, theta, rng)
    vu = factor_via_expmod_shortint(N - v, expmod, m, ctr, rng, theta, base)
    ctr.calls += 1
    part = expmod(vu, {v: 1})
    x, out = N, {}
    for p in part:
        while x % p == 0:
            out[p] = out.get(p, 0) + 1
            x //= p
    if x > 1:
        r = iroot(x, m)
        assert r is not None
        vr = factor_via_expmod_shortint(r, expmod, m, ctr, rng, theta, base)
        out = vmul(out, {p: e * m for p, e in vr.items()})
    return out
