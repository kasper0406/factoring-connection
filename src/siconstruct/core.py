"""Core exact helpers for agent `siconstruct` (phase 20): constructing y-smooth u in [kN, kN + h].

Operational code here sees only N and public parameters.  The only factorisations are trial divisions by the
public factor base (the smoothness test that defines a relation).  Floats are used only to choose toy parameters
(the L-path) and in the prediction models (`ht_local_density`), never in a decision that an output is valid.
"""
from __future__ import annotations

import math
import random
from bisect import bisect_left, bisect_right
from math import isqrt
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

KAPPA_STAR = (2.0 / 9.0) ** (1.0 / 3.0)          # (2/9)^{1/3} = 0.6057
NU_STAR = 2.0 / (3.0 * KAPPA_STAR)               # existence threshold 2/(3 kappa*) = 1.1006


# ----------------------------------------------------------------------------------------------- primes / smoothness
def primes_upto(B: int) -> List[int]:
    if B < 2:
        return []
    s = bytearray([1]) * (B + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, isqrt(B) + 1):
        if s[i]:
            s[i * i :: i] = bytearray(len(s[i * i :: i]))
    return [i for i in range(B + 1) if s[i]]


def factor_over(x: int, P: Sequence[int]) -> Optional[Dict[int, int]]:
    """Exponents of |x| over P if |x| is P-smooth (x != 0), else None.  Trial division only (exact)."""
    if x == 0:
        return None
    x = abs(x)
    out: Dict[int, int] = {}
    for p in P:
        if x == 1:
            break
        if x % p == 0:
            k = 0
            while x % p == 0:
                x //= p
                k += 1
            out[p] = k
    return out if x == 1 else None


def is_smooth(x: int, P: Sequence[int]) -> bool:
    return factor_over(x, P) is not None


# ----------------------------------------------------------------------------------------------- toy L-path
def lpath_params(bits: int, kappa: float = KAPPA_STAR, nu_factor: float = 1.25) -> Dict[str, float]:
    """Toy parameters on the Prop SI path with the o(1) terms dropped: y = L[1/3, kappa], h = L[2/3, nu],
    nu = nu_factor * 2/(3 kappa).  Returns ints y, h and the float diagnostics."""
    x = bits * math.log(2.0)
    lx = math.log(x)
    F13 = x ** (1 / 3) * lx ** (2 / 3)
    F23 = x ** (2 / 3) * lx ** (1 / 3)
    nu = nu_factor * 2.0 / (3.0 * kappa)
    y = int(round(math.exp(kappa * F13)))
    h = int(round(math.exp(nu * F23)))
    return {"bits": bits, "y": y, "h": h, "nu": nu, "kappa": kappa, "x": x,
            "log2h": math.log2(h), "u_N": x / math.log(y)}


# ----------------------------------------------------------------------------------------------- exact sieve
def _lg16(p: int) -> int:
    """floor(16 log2 p), exactly: the largest t with 2^t <= p^16."""
    return (p ** 16).bit_length() - 1


def sieve_window(X: int, W: int, P: Sequence[int], chunk: int = 1 << 22) -> Tuple[List[int], int]:
    """All P-smooth integers in [X, X + W), exactly (log sieve with integer weights + trial-division check).

    Completeness: a P-smooth z accumulates sum_p v_p(z) floor(16 log2 p) >= 16 log2 z - Omega(z)
    >= 16 (bitlen(X) - 1) - bitlen(X + W) =: T, so thresholding at T misses nothing.  Returns (list, points sieved).
    """
    out: List[int] = []
    T = 16 * (X.bit_length() - 1) - (X + W).bit_length()
    top = X + W
    pw: List[Tuple[int, int]] = []
    for p in P:
        lg = _lg16(p)
        q = p
        while q <= top:
            pw.append((q, lg))
            q *= p
    start = X
    while start < X + W:
        L = min(chunk, X + W - start)
        acc = np.zeros(L, dtype=np.int32)
        for q, lg in pw:
            off = (-start) % q
            if off < L:
                acc[off::q] += lg
        idx = np.nonzero(acc >= T)[0]
        for i in idx.tolist():
            z = start + i
            if is_smooth(z, P):
                out.append(z)
        start += L
    return out, W


# ----------------------------------------------------------------------------------------------- random smooth lists
def random_smooth_list(lo: int, hi: int, P: Sequence[int], count: int, rng: random.Random,
                       max_tries: int = 10 ** 8) -> List[int]:
    """`count` distinct P-smooth integers in [lo, hi], built as products of uniformly chosen primes of P.

    The distribution is not uniform on smooth numbers (it is N-oblivious, which is what the MITM model needs).
    Each accepted element costs O(number of prime draws) = L^{o(1)}."""
    seen = set()
    out: List[int] = []
    tries = 0
    plist = list(P)
    while len(out) < count:
        tries += 1
        if tries > max_tries:
            raise RuntimeError("random_smooth_list: too many tries")
        z = 1
        while z < lo:
            z *= rng.choice(plist)
        if z <= hi and z not in seen:
            seen.add(z)
            out.append(z)
    return out


# ----------------------------------------------------------------------------------------------- prediction model
def ht_local_density(X: float, y: int, P: Optional[Sequence[int]] = None) -> float:
    """Hildebrand–Tenenbaum saddle-point model of the local density of y-smooth integers near X
    (prediction only; floats).  dPsi/dX ~ alpha * Psi(X, y) / X, Psi ~ X^a zeta(a, y) / (a sqrt(2 pi phi2))."""
    if P is None:
        P = primes_upto(y)
    lX = math.log(X)
    lp = [math.log(p) for p in P]

    def phi1(a):  # sum log p / (p^a - 1)
        return sum(l / (math.exp(a * l) - 1.0) for l in lp)

    lo, hi = 1e-6, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if phi1(mid) > lX:
            lo = mid
        else:
            hi = mid
    a = (lo + hi) / 2
    logzeta = -sum(math.log1p(-math.exp(-a * l)) for l in lp)
    phi2 = sum(math.exp(a * l) * l * l / (math.exp(a * l) - 1.0) ** 2 for l in lp)
    logpsi = a * lX + logzeta - math.log(a) - 0.5 * math.log(2 * math.pi * phi2)
    return math.exp(logpsi - lX) * a
