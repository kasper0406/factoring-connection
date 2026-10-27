"""nativemod HARNESS (uses PARI via cypari2 and the factors).  Never imported by operational code.

- Eichler-Selberg Hurwitz moments C_k(N) and level-1 Hecke traces.
- Hecke matrices on M_2(Gamma_0(M)).
- Theta series of fixed quaternary lattices (level 11 genus, class-number-one lattices).
- Merel's Heilbronn set sizes.
- a_p of fixed curves, division polynomials, and the E2/E3 data.
"""
from __future__ import annotations

import math
from fractions import Fraction
from functools import lru_cache
from typing import Dict, List, Sequence, Tuple

import cypari2

pari = cypari2.Pari()
try:
    pari.allocatemem(512 * 10 ** 6, silent=True)
except TypeError:  # older cypari2
    pari.allocatemem(512 * 10 ** 6)
pari.default("nbthreads", 1)

CURVES = {"37a1": [0, 0, 1, -1, 0], "11a1": [0, -1, 1, -10, -20], "32a2": [0, 0, 0, -1, 0]}

# level-11 genus of even quaternary lattices of determinant 121 (Gram matrices of 2Q); theta_1 - theta_2 = 4 f_11a1
L11_1 = [[2, 1, 0, 0], [1, 6, 0, 0], [0, 0, 2, 1], [0, 0, 1, 6]]
L11_2 = [[4, -1, -1, 2], [-1, 4, 2, -1], [-1, 2, 4, 0], [2, -1, 0, 4]]
# post hoc (E1c(ii')): a third lattice found by a wider search; with L11_1, L11_2 it carries Siegel's genus average
L11_3 = [[2, -1, -1, -1], [-1, 2, 0, 0], [-1, 0, 8, -3], [-1, 0, -3, 8]]
# class-number-one lattices (Gram of 2Q for even ones; Z^4 is odd, Gram of Q)
D4_HURWITZ = [[2, 0, 0, 1], [0, 2, 0, 1], [0, 0, 2, 1], [1, 1, 1, 2]]   # 2*nrd on the Hurwitz order (D4)
A2A2 = [[2, 1, 0, 0], [1, 2, 0, 0], [0, 0, 2, 1], [0, 0, 1, 2]]


def sigma(n: int, k: int = 1) -> int:
    return int(pari.sigma(n, k))


def hclassno(n: int) -> Fraction:
    """Hurwitz class number H(n) (H(0) = -1/12, H(3) = 1/3, H(4) = 1/2)."""
    v = pari(f"qfbhclassno({n})")
    return Fraction(int(pari.numerator(v)), int(pari.denominator(v)))


def P_k(t: int, N: int, k: int) -> int:
    """[z^{k-2}] (1 - t z + N z^2)^{-1}."""
    c0, c1 = 1, t
    if k == 2:
        return 1
    for _ in range(k - 3):
        c0, c1 = c1, t * c1 - N * c0
    return c1


def hurwitz_moment(N: int, k: int) -> Fraction:
    tot = Fraction(0)
    t = 0
    while t * t < 4 * N:
        H = hclassno(4 * N - t * t)
        w = 1 if t == 0 else 2
        tot += w * P_k(t, N, k) * H
        t += 1
    return tot


@lru_cache(maxsize=None)
def _mf_level1(k: int):
    return pari(f"mfinit([1,{k}],1)")


def trace_TN_level1(N: int, k: int) -> int:
    mf = _mf_level1(k)
    if int(pari.mfdim(mf)) == 0:
        return 0
    return int(pari.trace(pari.mfheckemat(mf, N)))


@lru_cache(maxsize=None)
def _mf_full(M: int, k: int = 2):
    return pari(f"mfinit([{M},{k}],4)")


def hecke_full(M: int, n: int, k: int = 2):
    return pari.mfheckemat(_mf_full(M, k), n)


def theta_counts(G: Sequence[Sequence[int]], B: int, even: bool = True) -> List[int]:
    """r_L(n), n = 0..B, for the lattice with Gram matrix G; even=True means Q = x^T G x / 2."""
    M = pari.matrix(len(G), len(G), [x for row in G for x in row])
    if even:
        v = pari.qfrep(M, B, 1)
    else:
        v = pari.qfrep(M, B, 0)
    return [1] + [2 * int(x) for x in v]


def aut_order(G: Sequence[Sequence[int]]) -> int:
    M = pari.matrix(len(G), len(G), [x for row in G for x in row])
    return int(pari.qfauto(M)[0])


def coeffs_11a1(B: int) -> List[int]:
    return [int(c) for c in pari(f"mfcoefs(mfeigenbasis(mfinit([11,2],0))[1], {B})")]


def merel_size(n: int) -> int:
    """#{(a b; c d): ad - bc = n, a > b >= 0, d > c >= 0} (Merel 1994), O(n^2) loop over (a, d)... counted by b."""
    cnt = 0
    for a in range(1, n + 1):
        for b in range(0, a):
            # for each d >= 1 with c = (a d - n)/b in [0, d)
            if b == 0:
                if n % a == 0:
                    cnt += n // a          # d = n/a, c in [0, d)
                continue
            # need a d - n = b c, 0 <= c < d  <=>  n/a <= d  and  a d - n < b d  <=>  d < n/(a-b)
            dlo = -(-n // a)
            dhi = (n - 1) // (a - b)   # d (a-b) < n
            # among d in [dlo, dhi] count those with a d = n mod b
            if dhi < dlo:
                continue
            g = math.gcd(a, b)
            if n % g:
                continue
            bb = b // g
            if bb == 1:
                cnt += dhi - dlo + 1
                continue
            d0 = ((n // g) * pow(a // g, -1, bb)) % bb
            first = dlo + ((d0 - dlo) % bb)
            if first <= dhi:
                cnt += (dhi - first) // bb + 1
    return cnt


def ellap(curve: str, p: int) -> int:
    E = pari.ellinit(CURVES[curve])
    return int(pari.ellap(E, p))


@lru_cache(maxsize=None)
def _E(curve: str):
    return pari.ellinit(CURVES[curve])


def ap(curve: str, p: int) -> int:
    return int(pari.ellap(_E(curve), p))


def divpoly(curve: str, l: int) -> List[int]:
    return [int(c) for c in pari.Vecrev(pari.elldivpol(_E(curve), l))]


def kron(a: int, n: int) -> int:
    return int(pari.kronecker(a, n))
