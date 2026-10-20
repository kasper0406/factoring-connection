"""Small finite fields F_Q (Q = p^e) by lookup tables, and polynomials over them.

Used only for the small-characteristic POSITIVE CONTROL of the systeq experiment (the
BGJT/KW systematic-equation relation family).  Exact arithmetic; no floating point.

Elements of F_Q are the integers 0..Q-1, read as base-p digit vectors of a polynomial in w
modulo a fixed monic irreducible g(w) of degree e (the lexicographically first one).
Polynomials over F_Q are Python lists of ints, lowest degree first, without trailing zeros
(the zero polynomial is []).
"""
from __future__ import annotations

from itertools import product
from typing import List, Tuple

import numpy as np


def _poly_mod_p_irreducible(p: int, e: int) -> List[int]:
    """First monic irreducible of degree e over F_p (coefficients low->high, length e+1)."""
    if e == 1:
        return [0, 1]
    for tail in product(range(p), repeat=e):
        g = list(tail) + [1]
        if g[0] == 0:
            continue
        # irreducible iff no root-free factorisation: test by gcd with w^{p^i} - w, i <= e/2
        if _is_irreducible_fp(g, p):
            return g
    raise ValueError("no irreducible found")


def _fp_trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def _fp_divmod(a, b, p):
    a = list(a)
    inv = pow(b[-1], p - 2, p)
    q = [0] * max(0, len(a) - len(b) + 1)
    while len(a) >= len(b) and a:
        c = a[-1] * inv % p
        s = len(a) - len(b)
        q[s] = c
        for i, bi in enumerate(b):
            a[s + i] = (a[s + i] - c * bi) % p
        _fp_trim(a)
    return _fp_trim(q), a


def _fp_mul(a, b, p):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] = (r[i + j] + x * y) % p
    return _fp_trim(r)


def _fp_gcd(a, b, p):
    a, b = _fp_trim(list(a)), _fp_trim(list(b))
    while b:
        a, b = b, _fp_divmod(a, b, p)[1]
    return a


def _fp_powmod(base, k, mod, p):
    r = [1]
    base = _fp_divmod(base, mod, p)[1]
    while k:
        if k & 1:
            r = _fp_divmod(_fp_mul(r, base, p), mod, p)[1]
        base = _fp_divmod(_fp_mul(base, base, p), mod, p)[1]
        k >>= 1
    return r


def _is_irreducible_fp(g, p):
    e = len(g) - 1
    for i in range(1, e // 2 + 1):
        h = _fp_powmod([0, 1], p ** i, g, p)
        h = list(h) + [0] * max(0, 2 - len(h))
        h[1] = (h[1] - 1) % p
        if len(_fp_gcd(g, _fp_trim(h), p)) > 1:
            return False
    return True


class GF:
    """The field F_Q, Q = p^e, with full add/mul tables (Q <= ~400)."""

    def __init__(self, p: int, e: int):
        self.p, self.e, self.Q = p, e, p ** e
        Q = self.Q
        self.g = _poly_mod_p_irreducible(p, e)
        digits = [[(x // p ** i) % p for i in range(e)] for x in range(Q)]
        enc = lambda v: sum(int(c) * p ** i for i, c in enumerate(v))
        add = np.zeros((Q, Q), dtype=np.int32)
        mul = np.zeros((Q, Q), dtype=np.int32)
        for x in range(Q):
            for y in range(Q):
                add[x, y] = enc([(a + b) % p for a, b in zip(digits[x], digits[y])])
                prod_ = _fp_divmod(_fp_mul(_fp_trim(list(digits[x])), _fp_trim(list(digits[y])), p),
                                   self.g, p)[1]
                mul[x, y] = enc(prod_ + [0] * (e - len(prod_)))
        self.add, self.mul = add, mul
        self.neg = np.array([enc([(-c) % p for c in digits[x]]) for x in range(Q)], dtype=np.int32)
        self.inv = np.zeros(Q, dtype=np.int32)
        for x in range(1, Q):
            self.inv[x] = int(np.nonzero(mul[x] == 1)[0][0])
        self.elems = np.arange(Q, dtype=np.int32)

    # scalar helpers -------------------------------------------------------------------
    def a(self, x, y):
        return int(self.add[x, y])

    def m(self, x, y):
        return int(self.mul[x, y])

    def s(self, x, y):
        return int(self.add[x, self.neg[y]])

    def pw(self, x, k):
        r = 1
        while k:
            if k & 1:
                r = self.m(r, x)
            x = self.m(x, x)
            k >>= 1
        return r

    def powvec(self, k):
        """Vector x -> x^k over all of F_Q."""
        r = np.ones(self.Q, dtype=np.int32)
        b = self.elems.copy()
        while k:
            if k & 1:
                r = self.mul[r, b]
            b = self.mul[b, b]
            k >>= 1
        return r

    # polynomials ----------------------------------------------------------------------
    @staticmethod
    def trim(f):
        f = list(f)
        while f and f[-1] == 0:
            f.pop()
        return f

    def padd(self, f, g):
        n = max(len(f), len(g))
        f = list(f) + [0] * (n - len(f))
        g = list(g) + [0] * (n - len(g))
        return self.trim([self.a(x, y) for x, y in zip(f, g)])

    def psub(self, f, g):
        return self.padd(f, [int(self.neg[y]) for y in g])

    def pscale(self, c, f):
        return self.trim([self.m(c, x) for x in f])

    def pmul(self, f, g):
        if not f or not g:
            return []
        r = [0] * (len(f) + len(g) - 1)
        for i, x in enumerate(f):
            if x == 0:
                continue
            for j, y in enumerate(g):
                r[i + j] = self.a(r[i + j], self.m(x, y))
        return self.trim(r)

    def pdivmod(self, f, g):
        f = list(f)
        if not g:
            raise ZeroDivisionError
        inv = int(self.inv[g[-1]])
        q = [0] * max(0, len(f) - len(g) + 1)
        while len(f) >= len(g) and f:
            c = self.m(f[-1], inv)
            sft = len(f) - len(g)
            q[sft] = c
            for i, gi in enumerate(g):
                f[sft + i] = self.s(f[sft + i], self.m(c, gi))
            f = self.trim(f)
        return self.trim(q), f

    def pmod(self, f, g):
        return self.pdivmod(f, g)[1]

    def monic(self, f):
        return self.pscale(int(self.inv[f[-1]]), f) if f else f

    def pgcd(self, f, g):
        f, g = self.trim(f), self.trim(g)
        while g:
            f, g = g, self.pmod(f, g)
        return self.monic(f)

    def ppowmod(self, f, k, mod):
        r = [1]
        f = self.pmod(f, mod)
        while k:
            if k & 1:
                r = self.pmod(self.pmul(r, f), mod)
            f = self.pmod(self.pmul(f, f), mod)
            k >>= 1
        return r

    def peval(self, f, x):
        r = 0
        for c in reversed(f):
            r = self.a(self.m(r, x), c)
        return r

    def root_multiset(self, f) -> List[int]:
        """All roots of f in F_Q, with multiplicity (exact, by repeated division)."""
        f = self.trim(f)
        roots = []
        for x in range(self.Q):
            while len(f) > 1 and self.peval(f, x) == 0:
                roots.append(x)
                f, rem = self.pdivmod(f, [int(self.neg[x]), 1])
                assert not rem
        return roots

    def splits(self, f) -> bool:
        """True iff the nonzero polynomial f is a product of linear factors over F_Q."""
        f = self.trim(f)
        if not f:
            raise ValueError("zero polynomial")
        return len(self.root_multiset(f)) == len(f) - 1

    def irreducible_factors_ddf(self, f) -> List[Tuple[int, List[int]]]:
        """Distinct-degree factorisation of a squarefree monic f: list of (degree, product)."""
        f = self.monic(f)
        out, i, h = [], 1, [0, 1]
        while len(f) - 1 >= 2 * i:
            h = self.ppowmod(h, self.Q, f)
            g = self.pgcd(f, self.psub(h, [0, 1]))
            if len(g) > 1:
                out.append((i, g))
                f = self.pdivmod(f, g)[0]
                h = self.pmod(h, f)
            i += 1
        if len(f) > 1:
            out.append((len(f) - 1, f))
        return out


def split_count_formula(Q: int, D: int):
    """Number of split binary forms of degree D over F_Q, up to scalars, and the total.

    Split forms up to scalars <-> multisets of D points of P^1(F_Q): C(Q+D, D).
    All nonzero forms up to scalars: (Q^{D+1}-1)/(Q-1).
    """
    from math import comb
    return comb(Q + D, D), (Q ** (D + 1) - 1) // (Q - 1)
