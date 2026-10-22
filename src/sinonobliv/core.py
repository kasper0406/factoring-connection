"""Exact helpers for sinonobliv.

- primes, exact smoothness masks for arbitrary uint64 arrays (repeated exact division);
- exact empirical smooth-density curves (used only to form model predictions; floats appear only in ratios);
- polynomial identity families with forced factorisations (route 2);
- divisors of lcm(1..B) and certificate thresholds (route 3).

Nothing here sees a factor of N. Harness-only code (none needed) would be labelled as such.
"""
from __future__ import annotations

import math
from functools import lru_cache
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np


# ------------------------------------------------------------------------------------------------ primes / masks

@lru_cache(maxsize=None)
def primes_upto(n: int) -> Tuple[int, ...]:
    """All primes <= n (copied from sihard.core.primes_upto; owner rule 7)."""
    if n < 2:
        return ()
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i:: i] = bytearray(len(s[i * i:: i]))
    return tuple(i for i in range(n + 1) if s[i])


def smooth_mask_values(values: np.ndarray, y: int) -> np.ndarray:
    """Exact: mask[i] True iff values[i] (a positive integer < 2^63) is y-smooth.

    Repeated exact division of a uint64 cofactor array by every prime <= y.
    """
    v = np.asarray(values)
    if v.size == 0:
        return np.zeros(0, dtype=bool)
    if (v < 1).any():
        raise ValueError("values must be positive")
    if int(v.max()) >= (1 << 63):
        raise ValueError("values out of uint64-safe range")
    cof = v.astype(np.uint64).copy()
    for p in primes_upto(y):
        p64 = np.uint64(p)
        while True:
            q = cof // p64
            m = (q * p64) == cof
            if not m.any():
                break
            cof[m] = q[m]
    return cof == 1


def smooth_mask_window(lo: int, length: int, y: int) -> np.ndarray:
    """Exact mask for the consecutive window [lo, lo+length) (strided division, faster)."""
    if lo < 1 or lo + length >= (1 << 63):
        raise ValueError("window out of uint64 range")
    cof = np.arange(lo, lo + length, dtype=np.uint64)
    for p in primes_upto(y):
        p64 = np.uint64(p)
        start = (-lo) % p
        if start >= length:
            continue
        view = cof[start::p]
        while True:
            q = view // p64
            m = (q * p64) == view
            if not m.any():
                break
            view[m] = q[m]
        cof[start::p] = view
    return cof == 1


def is_smooth_int(z: int, y: int) -> bool:
    """Exact trial division for Python ints (any size)."""
    if z < 1:
        raise ValueError
    for p in primes_upto(y):
        while z % p == 0:
            z //= p
        if z == 1:
            return True
    return z == 1


# ------------------------------------------------------------------------------------- density curve (model)

class DensityCurve:
    """Empirical y-smooth density near 2^b on a half-bit grid, measured exactly on windows of `width` integers
    starting at round(2^b); log-linear interpolation in log2(size). Used only to form *predictions*.
    Sizes below the grid are clamped to the first grid point (densities there are close to 1 when b is small).
    """

    def __init__(self, y: int, b_lo: float, b_hi: float, width: int = 1 << 16, step: float = 0.5):
        self.y = y
        self.grid: List[float] = []
        self.dens: List[float] = []
        self.widths: List[int] = []
        b = b_lo
        while b <= b_hi + 1e-9:
            lo = max(int(round(2.0 ** b)), 1)
            w = min(width, max(16, lo))  # keep the window local: [2^b, 2^b + min(width, 2^b))
            m = smooth_mask_window(lo, w, y)
            c = int(m.sum())
            self.grid.append(b)
            self.dens.append(max(c, 0.5) / w)  # 0.5 floor avoids log(0); flagged by `zero_cells`
            self.widths.append(w)
            b += step
        self.zero_cells = [g for g, d, w in zip(self.grid, self.dens, self.widths) if d * w < 1]
        self._lg = np.array(self.grid)
        self._ld = np.log(np.array(self.dens))

    def __call__(self, size) -> np.ndarray:
        b = np.log2(np.asarray(size, dtype=float))
        b = np.clip(b, self._lg[0], self._lg[-1])
        return np.exp(np.interp(b, self._lg, self._ld))


# -------------------------------------------------------------------------------- route 2: identity families

def family_pieces(name: str, t: np.ndarray) -> List[np.ndarray]:
    """Forced factorisation of the member F(t) into integer pieces (all positive for t >= 2)."""
    t = np.asarray(t, dtype=object)
    if name == "cp2":            # t(t+1)
        return [t, t + 1]
    if name == "cp3":            # t(t+1)(t+2)
        return [t, t + 1, t + 2]
    if name == "aur":            # 4t^4 + 1 = (2t^2+2t+1)(2t^2-2t+1)   (Sophie Germain / Aurifeuillian)
        return [2 * t * t + 2 * t + 1, 2 * t * t - 2 * t + 1]
    if name == "sq1":            # t^2 + 1: irreducible over Z (control: no forced factorisation)
        return [t * t + 1]
    raise ValueError(name)


def family_value(name: str, t: np.ndarray) -> np.ndarray:
    v = np.ones(len(t), dtype=object)
    for p in family_pieces(name, t):
        v = v * p
    return v


FAMILY_DEGREE = {"cp2": 2, "cp3": 3, "aur": 4, "sq1": 2}
FAMILY_PIECES = {"cp2": 2, "cp3": 3, "aur": 2, "sq1": 1}


def family_t_for_size(name: str, X: int) -> int:
    """Least t >= 2 with F(t) >= X (exact integer search)."""
    lo, hi = 2, 2
    while int(family_value(name, np.array([hi], dtype=object))[0]) < X:
        hi *= 2
    while lo < hi:
        mid = (lo + hi) // 2
        if int(family_value(name, np.array([mid], dtype=object))[0]) < X:
            lo = mid + 1
        else:
            hi = mid
    return lo


def count_members_in_window(name: str, a: int, b: int) -> int:
    """Exact number of t >= 2 with a <= F(t) <= b (F increasing for t >= 2)."""
    t1 = family_t_for_size(name, a)
    t2 = family_t_for_size(name, b + 1)
    return t2 - t1


# ------------------------------------------------------------------------- route 3: divisors and thresholds

def lcm_upto_factored(B: int) -> Dict[int, int]:
    """lcm(1..B) as {p: e} with p^e <= B < p^{e+1}."""
    f = {}
    for p in primes_upto(B):
        e = 1
        while p ** (e + 1) <= B:
            e += 1
        f[p] = e
    return f


def divisors_from_factored(f: Dict[int, int]) -> np.ndarray:
    ds = [1]
    for p, e in f.items():
        new = []
        pk = 1
        for _ in range(e + 1):
            new.extend(d * pk for d in ds)
            pk *= p
        ds = new
    return np.array(sorted(ds), dtype=object)


def value_of_factored(f: Dict[int, int]) -> int:
    v = 1
    for p, e in f.items():
        v *= p ** e
    return v


def johnson_log_threshold(t: float, pi: float) -> float:
    """Unweighted Johnson / Coppersmith–Howgrave-Graham: certified iff ln gcd >= sqrt(ln V * ln P)."""
    return math.sqrt(t * pi)


def ch_log_threshold(t: float, pi: float, m: int) -> float:
    """Cohn–Heninger PACD region with one exact multiple P and m approximate samples, error R = e^t:
    divisor d = P^beta recoverable if R <= P^{beta^{(m+1)/m}}, i.e. ln d >= t^{m/(m+1)} pi^{1/(m+1)}.
    m = 1 gives the Johnson threshold. (Asymptotic region; granted at unit cost, as Theorem J grants Johnson.)"""
    return t ** (m / (m + 1.0)) * pi ** (1.0 / (m + 1.0))
