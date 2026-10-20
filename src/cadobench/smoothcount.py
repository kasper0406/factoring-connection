"""Exact smoothness of NFS norms at points of real special-q lattices (diagnostic harness).

Everything that decides smoothness is exact integer arithmetic (gmpy2):
  * Bernstein's batch test: with P = prod of primes <= B, x is B-smooth iff x | (P mod x)^(2^e) for
    2^e >= log2 x. P mod x is computed for a batch of x's with a product tree / remainder tree.
  * The lim-smooth part of x is gcd(x, (P_lim mod x)^(2^e) mod x).
Floating point is used only for model values (Dickman rho, Murphy alpha, logs) and fits.

Conventions (CADO): a = i*a0 + j*a1, b = i*b0 + j*b1 (then (a, b) -> (-a, -b) if b < 0);
i in [-I/2, I/2), j in [0, J). Side 0 = rational G (linear), side 1 = algebraic F; the special-q is on
side 1, so the side-1 norm is |F(a, b)|/q.
"""
from __future__ import annotations

import math
import os
import pickle
import random
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import gmpy2
import numpy as np
from gmpy2 import mpz

CACHE = os.path.expanduser("~/cadobench-work/primeprod")


# ---------------------------------------------------------------- primes and products

def primes_upto(B: int) -> np.ndarray:
    s = np.ones(B + 1, dtype=bool)
    s[:2] = False
    s[4::2] = False
    for p in range(3, int(B ** 0.5) + 1, 2):
        if s[p]:
            s[p * p::2 * p] = False
    return np.nonzero(s)[0]


def product_tree(xs: Sequence) -> List[List]:
    level = [mpz(x) for x in xs]
    tree = [level]
    while len(level) > 1:
        nxt = [level[k] * level[k + 1] for k in range(0, len(level) - 1, 2)]
        if len(level) % 2:
            nxt.append(level[-1])
        level = nxt
        tree.append(level)
    return tree


def prime_product(B: int) -> mpz:
    """prod of primes <= B (cached on disk)."""
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, f"P_{B}.bin")
    if os.path.exists(path):
        with open(path, "rb") as fh:
            return gmpy2.from_binary(fh.read())
    ps = primes_upto(B)
    # multiply in chunks of machine-size products first to cut Python overhead
    chunk = []
    acc = []
    for p in ps.tolist():
        chunk.append(p)
        if len(chunk) == 64:
            acc.append(math.prod(chunk))
            chunk = []
    if chunk:
        acc.append(math.prod(chunk))
    P = product_tree(acc)[-1][0]
    with open(path, "wb") as fh:
        fh.write(gmpy2.to_binary(P))
    return P


def remainders(P: mpz, xs: Sequence[int]) -> List[mpz]:
    """[P mod x for x in xs] via a remainder tree (xs nonzero)."""
    if not xs:
        return []
    tree = product_tree(xs)
    rems = [P % tree[-1][0]]
    for level in reversed(tree[:-1]):
        nxt = []
        for k, node in enumerate(level):
            nxt.append(rems[k // 2] % node)
        rems = nxt
    return rems


def _batches(n: int, P_bits: int, x_bits: int) -> int:
    # choose a batch so that the batch product is ~ P/4 (keeps P mod prod cheap)
    return max(256, min(1 << 16, (P_bits // 4) // max(x_bits, 1)))


def smooth_part(xs: Sequence[int], P: mpz) -> List[mpz]:
    """For each x > 0: the largest divisor of x composed of primes dividing P (exact)."""
    out: List[mpz] = []
    if not xs:
        return out
    xbits = max(int(x).bit_length() for x in xs)
    bs = _batches(len(xs), P.bit_length(), xbits)
    for s in range(0, len(xs), bs):
        chunk = [mpz(x) for x in xs[s:s + bs]]
        rs = remainders(P, chunk)
        for x, r in zip(chunk, rs):
            e = max(1, int(x.bit_length()).bit_length())  # 2^e >= bits(x)
            y = gmpy2.powmod(r, mpz(2) ** e, x)
            out.append(gmpy2.gcd(x, y) if y else x)
    return out


def is_smooth(xs: Sequence[int], P: mpz) -> List[bool]:
    return [sp == x for sp, x in zip(smooth_part(xs, P), map(mpz, xs))]


# ---------------------------------------------------------------- polynomials

def hom_eval(coeffs: Sequence[int], a: int, b: int) -> int:
    """F(a, b) = sum c_k a^k b^(d-k) (coeffs low degree first)."""
    d = len(coeffs) - 1
    return sum(c * a ** k * b ** (d - k) for k, c in enumerate(coeffs))


def _roots_mod_pk(coeffs: Sequence[int], p: int, kmax: int) -> List[int]:
    """R_k for k=1..kmax: #{x mod p^k : F(x,1) = 0} + #{y mod p^(k-1) : F(1, p*y) = 0 mod p^k}."""
    d = len(coeffs) - 1
    f = lambda x: sum(c * x ** k for k, c in enumerate(coeffs))
    rev = lambda y: sum(c * (p * y) ** (d - k) for k, c in enumerate(coeffs))  # F(1, p y)
    out = []
    aff = [x for x in range(p) if f(x) % p == 0]
    pro_mod = [0]  # y mod p^0
    pro = [y for y in pro_mod if rev(y) % p == 0]
    out.append(len(aff) + len(pro))
    for k in range(2, kmax + 1):
        pk, pk1 = p ** k, p ** (k - 1)
        aff = [x + t * pk1 for x in aff for t in range(p) if f(x + t * pk1) % pk == 0]
        pk2 = p ** (k - 2)
        pro = [y + t * pk2 for y in pro for t in range(p) if rev(y + t * pk2) % pk == 0] if k >= 2 else pro
        out.append(len(aff) + len(pro))
        if len(aff) + len(pro) == 0:
            out.extend([0] * (kmax - k))
            break
        if len(aff) + len(pro) > 20000:
            raise ValueError("too many roots to lift")
    return out


def murphy_alpha(coeffs: Sequence[int], pmax: int = 2000, divide_q: Optional[int] = None) -> float:
    """Murphy's alpha (natural log units) for coprime (a, b): sum_p (1/(p-1) - E[v_p F]) log p.

    E[v_p] = sum_k R_k / (p^(k-1) (p+1)), R_k counting projective roots mod p^k (exact counts; the
    k-sum is truncated where p^k > 1e6 or the roots die out).
    """
    alpha = 0.0
    for p in primes_upto(pmax).tolist():
        kmax = max(1, int(math.log(1e6) / math.log(p)))
        R = _roots_mod_pk(coeffs, p, kmax)
        ev = sum(Rk / (p ** k * (p + 1)) * p for k, Rk in enumerate(R, start=1))
        alpha += (1.0 / (p - 1) - ev) * math.log(p)
    return alpha


# ---------------------------------------------------------------- Dickman rho

_RHO = None
_H = 1e-4


def _rho_table(umax: float = 30.0):
    global _RHO
    if _RHO is None:
        n = int(umax / _H) + 2
        t = np.ones(n)
        m = int(1 / _H)
        # rho' (u) = -rho(u-1)/u ; trapezoid
        for k in range(m + 1, n):
            u0, u1 = (k - 1) * _H, k * _H
            t[k] = t[k - 1] - 0.5 * _H * (t[k - 1 - m] / u0 + t[k - m] / u1)
        _RHO = t
    return _RHO


def rho(u: float) -> float:
    if u <= 1.0:
        return 1.0
    t = _rho_table()
    x = u / _H
    i = int(x)
    if i + 1 >= len(t):
        return 0.0
    f = x - i
    return float(t[i] * (1 - f) + t[i + 1] * f)


# ---------------------------------------------------------------- lattice points

@dataclass
class QLattice:
    q: int
    rho: int
    a0: int
    b0: int
    a1: int
    b1: int
    I: int
    J: int

    def ab(self, i: int, j: int) -> Tuple[int, int]:
        a = i * self.a0 + j * self.a1
        b = i * self.b0 + j * self.b1
        return (-a, -b) if b < 0 else (a, b)

    def ij(self, a: int, b: int) -> Optional[Tuple[int, int]]:
        det = self.a0 * self.b1 - self.a1 * self.b0
        i_num = a * self.b1 - b * self.a1
        j_num = -a * self.b0 + b * self.a0
        if i_num % det or j_num % det:
            return None
        i, j = i_num // det, j_num // det
        if j < 0 or (j == 0 and i < 0):
            i, j = -i, -j
        return i, j


def sample_points(lat: QLattice, n: int, rng: random.Random) -> List[Tuple[int, int, int, int]]:
    """n uniform points (i, j) of the sieve region with j >= 1 and gcd(a, b) = 1."""
    out = []
    half = lat.I // 2
    while len(out) < n:
        i = rng.randrange(-half, half)
        j = rng.randrange(1, lat.J)
        a, b = lat.ab(i, j)
        if math.gcd(a, b) == 1:
            out.append((i, j, a, b))
    return out


def window_points(lat: QLattice, i_lo: int, i_hi: int, j_lo: int, j_hi: int):
    """All coprime (a,b) points with i in [i_lo, i_hi), j in [j_lo, j_hi)."""
    out = []
    for j in range(j_lo, j_hi):
        for i in range(i_lo, i_hi):
            if j == 0 and i <= 0:
                continue
            a, b = lat.ab(i, j)
            if math.gcd(a, b) == 1:
                out.append((i, j, a, b))
    return out


# ---------------------------------------------------------------- classification

@dataclass
class SideParams:
    lim: int
    lpb: int
    mfb: int


def norms(points, f0: Sequence[int], f1: Sequence[int], q: int) -> Tuple[List[int], List[int]]:
    n0, n1 = [], []
    for (_, _, a, b) in points:
        n0.append(abs(hom_eval(f0, a, b)))
        v = abs(hom_eval(f1, a, b))
        assert v % q == 0, "special-q does not divide the side-1 norm"
        n1.append(v // q)
    return n0, n1


class Smoother:
    """Holds the prime products for one parameter set."""

    def __init__(self, side: Sequence[SideParams]):
        self.side = list(side)
        self.PB = [prime_product(1 << s.lpb) for s in self.side]
        self.Plim = [prime_product(s.lim) for s in self.side]

    def classify(self, n: Sequence[int], s: int) -> Dict[str, List]:
        """Exact flags for side s: B-smooth (B = 2^lpb) and CADO-admissible (also cofactor <= 2^mfb)."""
        idx = [k for k, x in enumerate(n) if x != 0]
        xs = [n[k] for k in idx]
        sp = smooth_part(xs, self.PB[s])
        smooth = [False] * len(n)
        cado = [False] * len(n)
        cand = []
        for k, x, y in zip(idx, xs, sp):
            if y == x:
                smooth[k] = True
                cand.append(k)
        if cand:
            lp = smooth_part([n[k] for k in cand], self.Plim[s])
            for k, l in zip(cand, lp):
                cof = mpz(n[k]) // l
                cado[k] = cof.bit_length() <= self.side[s].mfb and cof < (mpz(1) << self.side[s].mfb)
        return {"smooth": smooth, "cado": cado}
