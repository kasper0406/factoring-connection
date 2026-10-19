"""OPERATIONAL toy NFS for the RK experiments (notes/rkrank_prereg.md).

Field: monic cubic f = fhat_{N,m} + (x - m)(r1 x + r0), so f(m) = N (as in src/hypD/rank.py).
Relations: coprime (a, b), b >= 1, |a| <= A, with a - b m and F(a, b) both B-smooth (exact tests).
Pairs whose F(a, b) is divisible by a prime p at a multiple root of f mod p are skipped (toy;
index-divisor bookkeeping not done).  Monic f has no projective primes.

Column layout (bit positions):
  rational block  [0, D_rat):     0 = sign(a - b m),  1.. = primes p <= B (v_p parity)
  high block      [D_rat, ncols): degree-one ideals (p, r); real-embedding signs
                                   [a - b*rho_j < 0] for each real root rho_j; quadratic characters
                                   (q, s), q > B, f(s) = 0 mod q simple.
The high block is eliminated first (gf2.Echelon pivots on the top bit), so
W = R ∩ F_2^{rational} = span of echelon vectors with pivot < D_rat.
Exact arithmetic throughout (real roots isolated with rational intervals and exact sign tests).
"""
from __future__ import annotations

import random
from fractions import Fraction
from math import gcd
from typing import Dict, List, Optional, Tuple

import gmpy2

from ..nfsfruit.smooth import FactorBase
from ..hypD.specialq import base_m_digits, random_m


def peval(f: List[int], x, p: Optional[int] = None):
    v = 0
    for co in reversed(f):
        v = v * x + co
        if p:
            v %= p
    return v


def deriv(f: List[int]) -> List[int]:
    return [i * f[i] for i in range(1, len(f))]


def roots_mod(f: List[int], p: int) -> List[int]:
    return [r for r in range(p) if peval(f, r, p) == 0]


def sign(x) -> int:
    return (x > 0) - (x < 0)


def real_root_intervals(f: List[int]) -> List[Tuple[Fraction, Fraction]]:
    """Isolating intervals (lo, hi), f(lo) f(hi) < 0, one simple real root each (f squarefree,
    degree 3).  Exact: Cauchy bound, then bisection on the critical-point partition."""
    lead = f[-1]
    Cb = 1 + max(abs(Fraction(c, lead)) for c in f[:-1])
    # critical points of a cubic: roots of f' (quadratic); isolate them by bisection on sign
    # of f' between -Cb and Cb, then split [-Cb, Cb] at points where f is monotone.
    df = deriv(f)
    pts = [-Cb, Cb]
    # f' quadratic: real roots exist iff disc >= 0
    A2, B2, C2 = df[2], df[1], df[0]
    disc = B2 * B2 - 4 * A2 * C2
    if disc > 0:
        # bracket each critical point to width <= 2^-60 * Cb using exact bisection on f'
        vertex = Fraction(-B2, 2 * A2)
        for lo, hi in ((-Cb, vertex), (vertex, Cb)):
            if sign(peval(df, lo)) * sign(peval(df, hi)) < 0:
                for _ in range(200):
                    mid = (lo + hi) / 2
                    if sign(peval(df, mid)) == sign(peval(df, lo)):
                        lo = mid
                    else:
                        hi = mid
                pts.append(lo)
                pts.append(hi)
    pts = sorted(set(pts))
    out = []
    for lo, hi in zip(pts, pts[1:]):
        flo, fhi = peval(f, lo), peval(f, hi)
        if flo == 0 or fhi == 0:
            raise ValueError("rational root at a breakpoint")
        if sign(flo) * sign(fhi) < 0:
            out.append((lo, hi))
    return out


def below_root(a: int, b: int, f: List[int], iv: Tuple[Fraction, Fraction]) -> bool:
    """Exact: is a/b < rho for the unique root rho of f in the isolating interval iv (b > 0)?"""
    lo, hi = iv
    q = Fraction(a, b)
    if q <= lo:
        return True
    if q >= hi:
        return False
    fq = peval(f, q)
    if fq == 0:
        raise ValueError("rational root")
    # f changes sign exactly once in (lo, hi): q > rho iff sign f(q) == sign f(hi)
    return sign(fq) != sign(peval(f, hi))


class ToyField:
    def __init__(self, N: int, B: int, nchar: int, seed, H: int = 4):
        self.N, self.B = N, B
        rng = random.Random(f"rkrank-field|{N}|{seed}")
        while True:
            self.m = random_m(N, rng)
            self.r0, self.r1 = rng.randrange(-H, H + 1), rng.randrange(-H, H + 1)
            c0, c1, c2 = base_m_digits(N, self.m)
            r0, r1, m = self.r0, self.r1, self.m
            self.f = [c0 - m * r0, c1 + r0 - m * r1, c2 + r1, 1]
            # irreducible over Q <=> no rational (= integer, monic) root: test divisors of f(0)
            if self._no_integer_root():
                break
        assert peval(self.f, self.m) == N
        self.df = deriv(self.f)
        self.fb = FactorBase(B)
        self.ps = self.fb.primes
        col = {"sign": 0}
        for p in self.ps:
            col[("r", p)] = len(col)
        self.D_rat = len(col)
        self.mult = {}                      # p -> set of multiple roots mod p
        for p in self.ps:
            rs = roots_mod(self.f, p)
            for r in rs:
                if peval(self.df, r, p) == 0:
                    self.mult.setdefault(p, set()).add(r)
                else:
                    col[("a", p, r)] = len(col)
        self.D_alg = len(col) - self.D_rat
        self.real = real_root_intervals(self.f)
        for j in range(len(self.real)):
            col[("rs", j)] = len(col)
        self.chars = []
        q = B
        while len(self.chars) < nchar:
            q = int(gmpy2.next_prime(q))
            if N % q == 0:
                continue
            for s in roots_mod(self.f, q):
                if peval(self.df, s, q) != 0 and len(self.chars) < nchar:
                    self.chars.append((q, s))
        for qs in self.chars:
            col[("c",) + qs] = len(col)
        self.col = col
        self.ncols = len(col)

    def _no_integer_root(self) -> bool:
        c0 = self.f[0]
        if c0 == 0:
            return False
        # integer roots divide c0; |root| <= Cauchy bound, scan small range and the m-region
        f = self.f
        bound = 1 + max(abs(c) for c in f[:-1])
        # cheap sufficient test: f has no root mod some small prime => no integer root
        for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
            if not roots_mod(f, p):
                return True
        # fall back: exact rational root search via sympy-free integer bisection on real roots
        for lo, hi in real_root_intervals(f):
            a, b = int(lo) - 1, int(hi) + 1
            for x in range(max(a, -bound), min(b, bound) + 1):
                if peval(f, x) == 0:
                    return False
                if b - a > 10 ** 6:
                    break
        return True

    def F(self, a: int, b: int) -> int:
        f = self.f
        return ((f[3] * a + f[2] * b) * a + f[1] * b * b) * a + f[0] * b ** 3

    def row_of(self, a: int, b: int) -> Optional[int]:
        """Row bitset of a coprime pair, or None if not smooth on both sides / skipped."""
        Fv = self.F(a, b)
        if Fv == 0:
            return None
        af = self.fb.factor(Fv)
        if af is None:
            return None
        t = a - b * self.m
        rf = self.fb.factor(t)
        if rf is None:
            return None
        col = self.col
        row = 1 if t < 0 else 0
        for p, e in rf.items():
            if e & 1:
                row ^= 1 << col[("r", p)]
        for p, e in af.items():
            r = a * pow(b, -1, p) % p          # b coprime to p (monic f, gcd(a,b)=1)
            if p in self.mult and r in self.mult[p]:
                return None                    # toy: skip index-divisor primes
            if e & 1:
                row ^= 1 << col[("a", p, r)]
        for j, iv in enumerate(self.real):
            if below_root(a, b, self.f, iv) is False:
                pass
            # a - b*rho < 0  <=>  a/b < rho
            if below_root(a, b, self.f, iv):
                row ^= 1 << col[("rs", j)]
        for (q, s) in self.chars:
            if gmpy2.legendre((a - b * s) % q, q) == -1:
                row ^= 1 << col[("c", q, s)]
        return row

    def sieve(self, A: int, bmax: int, target: int):
        """Sweep b = 1.. bmax, |a| <= A; stop once `target` relations are found."""
        rows, pairs = [], []
        for b in range(1, bmax + 1):
            for a in range(-A, A + 1):
                if gcd(a, b) != 1:
                    continue
                r = self.row_of(a, b)
                if r is not None:
                    rows.append(r)
                    pairs.append((a, b))
                    if len(rows) >= target:
                        return rows, pairs, b
        return rows, pairs, bmax

    def rational_vector(self, z: int) -> Optional[int]:
        """Rational-block vector of a B-smooth integer z (sign included), else None."""
        fz = self.fb.factor(z)
        if fz is None:
            return None
        v = 1 if z < 0 else 0
        for p, e in fz.items():
            if e & 1:
                v ^= 1 << self.col[("r", p)]
        return v


def gen_functionals(T: ToyField, factors: List[int]) -> List[int]:
    """HARNESS ONLY (uses the factors): the omega quadratic-character functionals on the
    rational block, w_i(p) = [(p/p_i) = -1], w_i(sign) = [(-1/p_i) = -1]."""
    out = []
    for pi in factors:
        w = 1 if gmpy2.legendre(pi - 1, pi) == -1 else 0
        for p in T.ps:
            if gmpy2.legendre(p, pi) == -1:
                w |= 1 << T.col[("r", p)]
        out.append(w)
    return out
