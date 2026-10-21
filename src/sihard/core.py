"""Exact helpers for sihard: primes, exact smoothness sieves on uint64 windows, Dickman rho (design only),
certificate counting for Johnson-radius decoders, two-list merge counts.

All smoothness decisions are exact integer computations. Floating point appears only in `dickman_rho`
(used to choose interval lengths, a design parameter) and in reported ratios.
"""
from __future__ import annotations

import math
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


def primes_upto(n: int) -> List[int]:
    if n < 2:
        return []
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i:: i] = bytearray(len(s[i * i:: i]))
    return [i for i in range(n + 1) if s[i]]


def smooth_mask(lo: int, length: int, primes: Sequence[int]) -> np.ndarray:
    """Exact: mask[i] is True iff lo + i is y-smooth (all prime factors in `primes`, which must be all primes <= y).

    Requires lo + length < 2^63. Repeated exact division on a uint64 cofactor array.
    """
    if lo < 1 or lo + length >= (1 << 63):
        raise ValueError("window out of uint64 range")
    cof = np.arange(lo, lo + length, dtype=np.uint64)
    for p in primes:
        p64 = np.uint64(p)
        start = (-lo) % p
        if start >= length:
            continue
        view = cof[start::p]
        # divide out p repeatedly (exact)
        while True:
            q = view // p64
            m = (q * p64) == view
            if not m.any():
                break
            view[m] = q[m]
        cof[start::p] = view
    return cof == 1


def is_smooth_int(z: int, y: int) -> bool:
    """Exact trial division (small y only)."""
    for p in primes_upto(y):
        while z % p == 0:
            z //= p
    return z == 1


def dickman_rho(u: float, step: float = 1e-4) -> float:
    """Numerical Dickman rho (design parameter only; not used for any claim)."""
    if u <= 1:
        return 1.0
    n = int(math.ceil(u / step))
    h = u / n
    m1 = int(round(1.0 / h))
    rho = np.ones(n + 1)
    # u*rho'(u) = -rho(u-1); integrate with trapezoid on rho(t-1)/t
    for i in range(m1 + 1, n + 1):
        t0, t1 = (i - 1) * h, i * h
        f0 = rho[i - 1 - m1] / t0
        f1 = rho[i - m1] / t1
        rho[i] = rho[i - 1] - 0.5 * h * (f0 + f1)
    return float(rho[n])


# ---------------------------------------------------------------- Johnson-radius certificate counting (E2)

def johnson_S(V: int, P: int) -> int:
    """Smallest integer S with (log S)^2 >= log V * log P (idealised unweighted Johnson condition).

    Numerical with 80-digit mpmath (not a certified computation); exact when V == P (then S = V).
    """
    import mpmath
    if V == P:
        return V
    with mpmath.workdps(80):
        s = mpmath.sqrt(mpmath.log(V) * mpmath.log(P))
        S = int(mpmath.ceil(mpmath.exp(s)))
        while S > 1 and mpmath.log(S - 1) >= s:
            S -= 1
        while mpmath.log(S) < s:
            S += 1
    return S


def num_divisors_at_least(g: int, S: int, fac: Dict[int, int]) -> int:
    """#{d | g : d >= S}, where `fac` is the factorisation of g (small primes)."""
    divs = [1]
    for p, e in fac.items():
        divs = [d * p ** k for d in divs for k in range(e + 1)]
    return sum(1 for d in divs if d >= S)


def factor_small(z: int, primes: Sequence[int]) -> Dict[int, int]:
    out: Dict[int, int] = {}
    for p in primes:
        if z % p == 0:
            e = 0
            while z % p == 0:
                z //= p
                e += 1
            out[p] = e
    if z != 1:
        raise ValueError("not smooth over given primes")
    return out


def johnson_counts(X: int, M: int, V: int, P: int, S: int, smooth_z: Sequence[int],
                   primes: Sequence[int]) -> Dict[str, float]:
    """Exact E2 quantities for one (P, S, V).

    smooth_z: all y-smooth integers in [X, X + M + V].
    A  = (1/M) sum_{x in [X, X+M)} #{z in [x, x+V] : gcd(z,P) >= S, z smooth}
    B  = ((V+1)/M) #{z in [X, X+M+V] : gcd(z,P) >= S, z smooth}
    Bd = ((V+1)/M) sum_{z smooth} #{d | gcd(z,P) : d >= S}
    """
    A_num = 0
    cnt = 0
    dsum = 0
    for z in smooth_z:
        g = math.gcd(z, P)
        if g >= S:
            # number of x in [X, X+M) with x <= z <= x + V
            lo = max(X, z - V)
            hi = min(X + M - 1, z)
            if hi >= lo:
                A_num += hi - lo + 1
            cnt += 1
            dsum += num_divisors_at_least(g, S, factor_small(g, primes))
    return {"A": A_num / M, "B": (V + 1) * cnt / M, "Bd": (V + 1) * dsum / M, "n_cert_smooth": cnt}


# ---------------------------------------------------------------- two-list merge (E3)

def merge_interval_average(A: Sequence[int], B: Sequence[int], X: int, h: int) -> Tuple[float, float]:
    """Exact average over x in [X, 2X) of #{(a,b) : ab in [x, x+h]}, and the bound |A||B|(h+1)/X."""
    a = np.array(A, dtype=np.uint64)
    b = np.array(B, dtype=np.uint64)
    w = np.multiply.outer(a, b).ravel().astype(np.int64)  # caller guarantees products < 2^62
    lo = np.maximum(np.int64(X), w - np.int64(h))
    hi = np.minimum(np.int64(2 * X - 1), w)
    c = np.clip(hi - lo + 1, 0, None)
    total = int(c.sum())
    return total / X, len(A) * len(B) * (h + 1) / X


def merge_bohr_counts(A: Sequence[int], B: Sequence[int], Ns: Iterable[int], h: int) -> List[int]:
    """For each N: #{(a,b) : ab mod N in [0, h], ab >= N}."""
    a = np.array(A, dtype=np.uint64)
    b = np.array(B, dtype=np.uint64)
    prods = np.multiply.outer(a, b).ravel()
    out = []
    for N in Ns:
        r = prods % np.uint64(N)
        out.append(int(((r <= np.uint64(h)) & (prods >= np.uint64(N))).sum()))
    return out


def smooth_in_range(lo: int, hi: int, y: int) -> List[int]:
    """All y-smooth integers in [lo, hi) (exact)."""
    m = smooth_mask(lo, hi - lo, primes_upto(y))
    return [lo + int(i) for i in np.nonzero(m)[0]]


def lemma_fp(x: int, h: int, y: int) -> int:
    """Lemma FP: if h >= 2x/(y-1) (or x <= y), return a y-smooth z = 2^a * ceil(x / 2^a) in [x, x+h] (poly time).

    y = 2 is allowed (then h >= 2x and z is a power of 2)."""
    if x <= y:
        return x
    if y < 2 or h * (y - 1) < 2 * x:
        raise ValueError("needs h >= 2x/(y-1)")
    if y == 2:
        return 1 << (x - 1).bit_length()  # least power of 2 >= x, <= 2x - 1 <= x + h
    a = 0
    while (1 << (a + 1)) * (y - 1) <= 2 * x:  # largest 2^a with 2^a (y-1) <= 2x, so 2^a >= x/(y-1)
        a += 1
    d = 1 << a
    return d * (-(-x // d))
