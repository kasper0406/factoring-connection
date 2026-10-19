"""Toy test of the rank hypothesis RK (experiment E3 of notes/hypD_prereg.md).

OPERATIONAL: a toy cubic NFS (monic base-m cubic + (x - m)(r1 x + r0)) collects pairs (a, b)
with a - b m and F(a, b) both B-smooth.  Row of a pair over GF(2):
    [sign(a - b m)] [v_p(a - b m), p <= B] [v_P(a - b alpha), P = (p, r) degree one, p <= B]
    [quadratic characters (a - b s | q) for nchar pairs (q, s), f(s) = 0 mod q, q > B]
(index-divisor subtleties ignored: toy).  A "rational" target vector has support only on
the sign and rational-prime columns.  `in_rowspace` decides membership by elimination.
Dixon vectors come from random x with x^2 mod N B-smooth (OPERATIONAL, never factors N).
"""
from __future__ import annotations

import random
from math import gcd
from typing import Dict, List, Optional, Tuple

import gmpy2

from ..nfsfruit.smooth import FactorBase, primes_upto
from .specialq import base_m_digits, random_m


def poly_coeffs(N: int, m: int, r0: int, r1: int) -> List[int]:
    c0, c1, c2 = base_m_digits(N, m)
    # x^3 + c2 x^2 + c1 x + c0 + (x - m)(r1 x + r0)
    return [c0 - m * r0, c1 + r0 - m * r1, c2 + r1, 1]     # f[0..3]


def poly_eval(f, x, p=None):
    v = 0
    for co in reversed(f):
        v = v * x + co
        if p:
            v %= p
    return v


def roots_mod(f, p) -> List[int]:
    return [r for r in range(p) if poly_eval(f, r, p) == 0]


class ToyNFS:
    def __init__(self, N: int, B: int, nchar: int, seed: int, H: int = 4):
        self.N, self.B = N, B
        rng = random.Random(f"hypD-rank|{N}|{seed}")
        while True:
            self.m = random_m(N, rng)
            self.r0, self.r1 = rng.randrange(-H, H + 1), rng.randrange(-H, H + 1)
            self.f = poly_coeffs(N, self.m, self.r0, self.r1)
            if all(poly_eval(self.f, x) != 0 for x in range(-50, 51)):
                break                       # crude irreducibility screen (toy)
        assert poly_eval(self.f, self.m) % N == 0
        self.fb = FactorBase(B)
        self.ps = self.fb.primes
        self.col = {}
        self.col["sign"] = 0
        for p in self.ps:
            self.col[("r", p)] = len(self.col)
        for p in self.ps:
            for r in roots_mod(self.f, p):
                self.col[("a", p, r)] = len(self.col)
        self.chars = []
        q = B + 1
        while len(self.chars) < nchar:
            q = int(gmpy2.next_prime(q))
            for s in roots_mod(self.f, q) if q < 20000 else []:
                if len(self.chars) < nchar:
                    self.chars.append((q, s))
        for qs in self.chars:
            self.col[("c",) + qs] = len(self.col)
        self.ncols = len(self.col)
        self.rows: List[int] = []
        self.pairs: List[Tuple[int, int]] = []

    def F(self, a, b):
        f = self.f
        return sum(f[i] * a ** i * b ** (3 - i) for i in range(4))

    def row_of(self, a: int, b: int) -> Optional[int]:
        t = a - b * self.m
        rf = self.fb.factor(t)
        if rf is None:
            return None
        Fv = self.F(a, b)
        if Fv == 0:
            return None
        af = self.fb.factor(Fv)
        if af is None:
            return None
        row = 0
        if t < 0:
            row ^= 1 << self.col["sign"]
        for p, e in rf.items():
            if e & 1:
                row ^= 1 << self.col[("r", p)]
        for p, e in af.items():
            if e & 1:
                if b % p == 0:
                    return None          # projective prime: toy skips the pair
                r = a * pow(b, -1, p) % p
                key = ("a", p, r)
                if key not in self.col:
                    return None          # (index divisor / ramified) toy skips
                row ^= 1 << self.col[key]
        for (q, s) in self.chars:
            if gmpy2.legendre((a - b * s) % q, q) == -1:
                row ^= 1 << self.col[("c", q, s)]
        return row

    def sieve(self, A: int, Bb: int, target_rows: int):
        for b in range(1, Bb + 1):
            for a in range(-A, A + 1):
                if gcd(a, b) != 1:
                    continue
                r = self.row_of(a, b)
                if r is not None:
                    self.rows.append(r)
                    self.pairs.append((a, b))
                    if len(self.rows) >= target_rows:
                        return

    def rational_vector(self, sign: int, exps: Dict[int, int]) -> int:
        v = (1 << self.col["sign"]) if sign < 0 else 0
        for p, e in exps.items():
            if e & 1:
                v ^= 1 << self.col[("r", p)]
        return v


class RowSpace:
    def __init__(self, rows: List[int]):
        self.piv: Dict[int, int] = {}
        for r in rows:
            self.add(r)

    def reduce(self, r: int) -> int:
        while r:
            c = r.bit_length() - 1
            if c in self.piv:
                r ^= self.piv[c]
            else:
                return r
        return 0

    def add(self, r: int):
        r = self.reduce(r)
        if r:
            self.piv[r.bit_length() - 1] = r

    @property
    def rank(self):
        return len(self.piv)

    def contains(self, v: int) -> bool:
        return self.reduce(v) == 0


def dixon_vectors(N: int, fb: FactorBase, count: int, seed: int, max_trials: int):
    rng = random.Random(f"hypD-dixon|{N}|{seed}")
    out, trials = [], 0
    while len(out) < count and trials < max_trials:
        trials += 1
        x = rng.randrange(2, N - 1)
        if gcd(x, N) != 1:
            continue
        z = x * x % N
        fz = fb.factor(z)
        if fz is not None:
            out.append(fz)
    return out, trials
