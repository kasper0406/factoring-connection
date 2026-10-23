"""Uniform-r zero sets of PRIM nodes at a prime p, their coset structure, and the KS bound.  DIAGNOSTIC (uses p).

Z_v(p) = {r in Z/(p-1) : node v = 0 mod p when every affine gate exponent a*N + b is replaced by a*r + b}.
At the true partner q, r = N mod (p-1) = q mod (p-1)  (Theorem U of notes/prim_results.md).
"""
from __future__ import annotations

from functools import lru_cache
from math import ceil, e as EULER, floor, log2
from typing import List, Tuple

import numpy as np
import sympy

from .prog import Ex


def primes_b(b: int) -> List[int]:
    """b-bit primes with the top two bits set: [3*2^(b-2), 2^b)."""
    return list(sympy.primerange(3 << (b - 2), 1 << b))


@lru_cache(maxsize=None)
def divisors(n: int) -> Tuple[int, ...]:
    return tuple(sympy.divisors(n))


def zero_set(node: Ex, p: int) -> np.ndarray:
    r = np.arange(p - 1, dtype=np.int64)
    return node.vec(p, r) == 0


def coset_union(Z: np.ndarray, L: float) -> np.ndarray:
    """Union of all cosets r0 + m*Z/n (m | n, coset size n/m >= L) contained in Z (boolean array of length n)."""
    n = Z.size
    V = np.zeros(n, dtype=bool)
    if not Z.any():
        return V
    for m in divisors(n):
        if n // m < L:
            continue
        full = Z.reshape(n // m, m).all(axis=0)          # residues r0 mod m whose whole coset lies in Z
        if full.any():
            V |= np.tile(full, n // m)
    return V


def max_vanishing_coset(Z: np.ndarray) -> int:
    """size of the largest coset (of a subgroup of Z/n) contained in Z; 0 if Z is empty."""
    n = Z.size
    if not Z.any():
        return 0
    best = 1
    for m in divisors(n):
        if n // m <= best:
            continue
        if Z.reshape(n // m, m).all(axis=0).any():
            best = n // m
    return best


def ks_threshold(n: int) -> float:
    """coset-size threshold of Lemma KS: n^{1/e}."""
    return n ** (1.0 / EULER)


def ks_bound(n: int, t: int) -> float:
    """Lemma KS residual bound: 2*lam*n^{1 - 1/(e*(1+lam)*(t-1))}, lam = floor(log2 t) (t >= 2)."""
    lam = int(floor(log2(t)))
    return 2 * lam * n ** (1.0 - 1.0 / (EULER * (1 + lam) * (t - 1)))


def kelley_bound(n: int, t: int, C: int) -> float:
    """Kelley 2016 Thm 2.3 on a cyclic group of order n: R <= 2 n^{1-1/(t-1)} C^{1/(t-1)}."""
    return 2 * n ** (1 - 1 / (t - 1)) * max(C, 1) ** (1 / (t - 1))


def L_cube(n: int) -> int:
    return int(ceil(n ** (1.0 / 3.0)))


def order_mod(a: int, p: int) -> int:
    return int(sympy.n_order(a, p))


def gate_subexprs(ex) -> list:
    """all affine gate / pure-power subexpressions of a node (post-order)."""
    from .prog import Add, G, Mul, PPow, Sub
    out = []
    if isinstance(ex, (Add, Sub, Mul)):
        out += gate_subexprs(ex.x) + gate_subexprs(ex.y)
    elif isinstance(ex, (G, PPow)):
        out += gate_subexprs(ex.x)
        out.append(ex)
    return out


def frozen_set(ex, p: int, T: int = 64) -> np.ndarray:
    """POST HOC class FROZEN: r such that every gate output at r equals its output at some small integer t,
    |t| <= T (t taken mod p-1).  On such r the node equals the integer tower evaluated at t, mod p."""
    n = p - 1
    r = np.arange(n, dtype=np.int64)
    gates = [g for g in gate_subexprs(ex)]
    if not gates:
        return np.zeros(n, dtype=bool)
    outs = [g.vec(p, r) for g in gates]
    F = np.zeros(n, dtype=bool)
    for t in range(-T, T + 1):
        tt = t % n
        m = np.ones(n, dtype=bool)
        for o in outs:
            m &= o == o[tt]
        F |= m
    return F


def ks_bound_compress(n: int, t: int, L: float) -> float:
    """Lemma KS (codex prim_1 (12), degree compression): residual roots outside zero cosets of order >= L
    number at most 2n / (floor(n/L) + 1)^{1/(t-1)}."""
    return 2 * n / (int(n // L) + 1) ** (1.0 / (t - 1))
