"""OPERATIONAL observables g: (Z/N)^* -> C, evaluated on numpy arrays of residues y in [0, N).

Each takes (y: np.ndarray[int64], N: int) and returns a complex or float array. Only y and N are used
(no factors, no group orders). P11 uses a public smoothness bound only.
"""
from __future__ import annotations

import math

import gmpy2
import numpy as np

SMALL_PRIMES_50 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]


def jacobi_vec(y: np.ndarray, N: int) -> np.ndarray:
    return np.fromiter((gmpy2.jacobi(int(v), N) for v in y), dtype=np.float64, count=len(y))


def obs_comb(y, N, c):
    return (y == c).astype(np.float64)


def obs_additive(y, N, h):
    # e(h y / N); h*y < 2^63 for N < 2^31
    return np.exp(2j * np.pi * ((h * y) % N) / N)


def obs_interval(y, N):
    return np.where(2 * y < N, 1.0, -1.0)


def obs_lsb(y, N):
    return np.where(y & 1, -1.0, 1.0)


def obs_jacobi(y, N):
    return jacobi_vec(y, N)


def obs_shifted_jacobi(y, N):
    return jacobi_vec((y + 1) % N, N)


def obs_jacobi_interval(y, N):
    return jacobi_vec(y, N) * obs_interval(y, N)


def obs_lsb_square(y, N):
    return obs_lsb((y * y) % N, N)


def obs_lsb_inverse(y, N):
    inv = np.fromiter((pow(int(v), -1, N) for v in y), dtype=np.int64, count=len(y))
    return obs_lsb(inv, N)


def obs_smooth(y, N, B=50, centre=True):
    """Indicator (centred by its sample mean unless centre=False) that the integer y in [1, N) is B-smooth."""
    z = y.copy()
    for ell in (p for p in SMALL_PRIMES_50 if p <= B):
        while True:
            m = (z % ell == 0) & (z > 1)
            if not m.any():
                break
            z[m] //= ell
    ind = (z == 1).astype(np.float64)
    return ind - ind.mean() if centre else ind


def lcm_upto(B: int) -> int:
    M = 1
    for k in range(2, B + 1):
        M = M * k // math.gcd(M, k)
    return M


def obs_oa(y, N, B1=30):
    """[y^M = 1 mod N], M = lcm(1..B1): an order-annihilation (p-1 style) test."""
    M = lcm_upto(B1)
    return np.fromiter((1.0 if pow(int(v), M, N) == 1 else 0.0 for v in y), dtype=np.float64,
                       count=len(y))


def obs_hash(y, N, salt=0x9E3779B97F4A7C15):
    """Pseudorandom +-1 function of y (splitmix64 finaliser): random-function baseline."""
    z = (y.astype(np.uint64) + np.uint64(salt))
    with np.errstate(over="ignore"):
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        z = z ^ (z >> np.uint64(31))
    return np.where(z & np.uint64(1), -1.0, 1.0)


PANEL = ["P1_comb", "P2_additive", "P3_interval", "P4_lsb", "P5_jacobi", "P6_shifted_jacobi",
         "P7_jacobi_interval", "P8_lsb_square", "P9_lsb_inverse", "P10_smooth", "P11_oa", "P12_hash"]


def evaluate(name: str, y: np.ndarray, N: int, rng) -> np.ndarray:
    if name == "P1_comb":
        return obs_comb(y, N, int(y[rng.randrange(len(y))]))
    if name == "P2_additive":
        while True:
            h = rng.randrange(1, N)
            if math.gcd(h, N) == 1:
                return obs_additive(y, N, h)
    f = {"P3_interval": obs_interval, "P4_lsb": obs_lsb, "P5_jacobi": obs_jacobi,
         "P6_shifted_jacobi": obs_shifted_jacobi, "P7_jacobi_interval": obs_jacobi_interval,
         "P8_lsb_square": obs_lsb_square, "P9_lsb_inverse": obs_lsb_inverse,
         "P10_smooth": obs_smooth, "P11_oa": obs_oa, "P12_hash": obs_hash}[name]
    return f(y, N)
