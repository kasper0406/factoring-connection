"""OPERATIONAL factor-blind evaluators: input N only. No factorisation of N, no p, q.
a(N, Delta^k) mod 2 for k = 9 (x^2+8y^2) and k = 11 (x^2+2y^2+8z^2), all variables odd positive."""
from math import isqrt
import numpy as np
import cypari2
pari = cypari2.Pari()

def eval_d9(N):
    """parity of #{x,y odd > 0 : x^2 + 8 y^2 = N}; cost ~ sqrt(N/8)/2 square tests."""
    ys = np.arange(1, isqrt(N // 8) + 1, 2, dtype=object)
    c = 0
    for y in range(1, isqrt(N // 8) + 1, 2):
        m = N - 8 * y * y
        r = isqrt(m)
        if r * r == m and r & 1: c ^= 1
    return c

def _half_sum_chi8(m):
    """(1/2) sum_{d|m} (-8/d) for odd m = 3 mod 8 (number of x,y odd >0 with x^2+2y^2=m). Factors m (not N)."""
    f = pari.factor(m)
    tot = 1
    for l, e in zip(f[0], f[1]):
        l, e = int(l), int(e)
        chi = int(pari.kronecker(-8, l))
        tot *= sum(chi ** i for i in range(e + 1))
    return tot // 2

def eval_d11(N):
    """parity of #{x,y,z odd>0: x^2+2y^2+8z^2 = N}; loops z, factors N - 8z^2 (never N). cost ~ sqrt(N/8)/2 factorings."""
    c = 0
    z = 1
    while 8 * z * z < N:
        m = N - 8 * z * z
        if m % 8 == 3: c ^= _half_sum_chi8(m) & 1
        z += 2
    return c
