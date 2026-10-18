"""OPERATIONAL: relation generators for toy congruence-of-squares factoring.

Every family sees only N, public parameters and a seed.  Each returns a System whose
`rows` are GF(2) parity bitmasks (one per relation) and whose `congruence(S)` turns a
kernel vector S (bitmask over relations) into (X, Y) with X^2 = Y^2 (mod N), or None if
the algebraic product is not a square (NFS only).  Fruitfulness is 1 < gcd(X - Y, N) < N.

Families
  dixon   : x uniform in [2, N-2], value x^2 mod N               (positive control)
  trivial : x uniform in [2, isqrt(N)], value x^2 (no reduction)  (negative control)
  qs      : x = isqrt(N) + 1 + i, |i| <= M, value x^2 - N         (quadratic sieve)
  lin_base: degree-1 NFS, L1 = x - m y, L2 = C1 x + C0 y, m ~ sqrt(N), n = C1 m + C0
  lin_rand: degree-1 NFS, m uniform mod N, L1, L2 random short vectors of
            {(c1,c0): c1 m + c0 = 0 (mod N)}                     (LV-style randomization)
  nfs3_base: monic cubic base-m NFS (R = 0)
  nfs3_rand: monic cubic LV-randomized NFS, f = fhat + (x - m)(r1 x + r0), |r_i| <= (Hmul*m) >> Hshift
"""
from __future__ import annotations
import random
from math import gcd, isqrt
from typing import Dict, List, Optional, Tuple

import gmpy2
from flint import fmpz_poly, nmod_poly

from .smooth import FactorBase
from .cubic import Order


class FactorFound(Exception):
    def __init__(self, g):
        super().__init__(g)
        self.g = g


class System:
    def __init__(self, N, family, params):
        self.N = N
        self.family = family
        self.params = params
        self.rows: List[int] = []
        self.ncols = 0
        self.info: Dict = {}

    def _inv_check(self, x):
        g = gcd(x % self.N, self.N)
        if 1 < g < self.N:
            raise FactorFound(g)


def _rng(family, N, seed):
    return random.Random(f"nfsfruit-op|{family}|{N}|{seed}")


# ---------------------------------------------------------------- one-sided families
class OneSided(System):
    """Relations (x_i, factorisation of value_i) with x_i^2 = value_i (mod N)."""

    def __init__(self, N, family, params, B):
        super().__init__(N, family, params)
        self.fb = FactorBase(B)
        self.col = {p: i + 1 for i, p in enumerate(self.fb.primes)}  # col 0 = sign
        self.ncols = len(self.fb.primes) + 1
        self.xs: List[int] = []
        self.facs: List[Dict[int, int]] = []
        self.signs: List[int] = []

    def add(self, x, value):
        fac = self.fb.factor(value)
        if fac is None:
            return False
        r = 1 if value < 0 else 0
        for p, e in fac.items():
            if e & 1:
                r |= 1 << self.col[p]
        self.xs.append(x)
        self.facs.append(fac)
        self.signs.append(-1 if value < 0 else 1)
        self.rows.append(r)
        return True

    def congruence(self, S: int) -> Tuple[int, int]:
        N = self.N
        X = 1
        E: Dict[int, int] = {}
        i = 0
        while S:
            if S & 1:
                X = X * self.xs[i] % N
                for p, e in self.facs[i].items():
                    E[p] = E.get(p, 0) + e
            S >>= 1
            i += 1
        Y = 1
        for p, e in E.items():
            assert e % 2 == 0
            Y = Y * pow(p, e // 2, N) % N
        return X, Y


def dixon(N, seed, B, target_excess, max_trials):
    sysm = OneSided(N, "dixon", dict(B=B), B)
    rng = _rng("dixon", N, seed)
    accidental = []
    t = 0
    while len(sysm.rows) < sysm.ncols + target_excess and t < max_trials:
        t += 1
        x = rng.randrange(2, N - 1)
        g = gcd(x, N)
        if g != 1:                      # ACCIDENTAL nonunit (prob ~ 2/sqrt(N) per trial)
            accidental.append(g)
            continue
        sysm.add(x, x * x % N)
    sysm.info = dict(trials=t, accidental=len(accidental))
    return sysm


def trivial(N, seed, B, target_excess, max_trials):
    sysm = OneSided(N, "trivial", dict(B=B), B)
    rng = _rng("trivial", N, seed)
    t = 0
    r = isqrt(N)
    seen = set()
    while len(sysm.rows) < sysm.ncols + target_excess and t < max_trials:
        t += 1
        x = rng.randrange(2, r)
        if x in seen:
            continue
        seen.add(x)
        sysm.add(x, x * x)
    sysm.info = dict(trials=t)
    return sysm


def qs(N, seed, B, M):
    sysm = OneSided(N, "qs", dict(B=B, M=M), B)
    s = isqrt(N) + 1
    for i in range(-M, M + 1):
        x = s + i
        v = x * x - N
        if v == 0:
            raise FactorFound(x)
        sysm.add(x, v)
    sysm.info = dict(trials=2 * M + 1)
    return sysm


# ---------------------------------------------------------------- degree-1 NFS
def _gauss_reduce(u, v):
    """Lagrange-Gauss reduction of a 2-dim integer lattice basis (exact)."""
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]

    def dot(a, b):
        return a[0] * b[0] + a[1] * b[1]
    if n2(u) > n2(v):
        u, v = v, u
    while True:
        # mu = round(dot/n2) exactly
        num, den = dot(u, v), n2(u)
        mu = (2 * num + den) // (2 * den)
        v = (v[0] - mu * u[0], v[1] - mu * u[1])
        if n2(v) >= n2(u):
            return u, v
        u, v = v, u


class Linear(System):
    """Two linear forms L1 = c1 x + c0 y, L2 = d1 x + d0 y with common root (m:1) mod N.
    d1 L1(a,b) = c1 L2(a,b) (mod N)."""

    def __init__(self, N, family, params, L1, L2, B):
        super().__init__(N, family, params)
        self.L1, self.L2 = L1, L2
        self.fb = FactorBase(B)
        P = len(self.fb.primes)
        self.col1 = {p: 1 + i for i, p in enumerate(self.fb.primes)}
        self.col2 = {p: 2 + P + i for i, p in enumerate(self.fb.primes)}
        # cols: 0 sign1, 1..P primes side1, P+1 sign2, P+2..2P+1 primes side2, 2P+2 parity |S|
        self.sign2 = P + 1
        self.par = 2 * P + 2
        self.ncols = 2 * P + 3
        self.f1: List[Dict[int, int]] = []
        self.f2: List[Dict[int, int]] = []
        self.ab: List[Tuple[int, int]] = []

    def try_pair(self, a, b):
        v1 = self.L1[0] * a + self.L1[1] * b
        if v1 == 0:
            return False
        fa = self.fb.factor(v1)
        if fa is None:
            return False
        v2 = self.L2[0] * a + self.L2[1] * b
        if v2 == 0:
            return False
        fb2 = self.fb.factor(v2)
        if fb2 is None:
            return False
        r = 1 << self.par
        if v1 < 0:
            r |= 1
        if v2 < 0:
            r |= 1 << self.sign2
        for p, e in fa.items():
            if e & 1:
                r |= 1 << self.col1[p]
        for p, e in fb2.items():
            if e & 1:
                r |= 1 << self.col2[p]
        self.rows.append(r)
        self.f1.append(fa)
        self.f2.append(fb2)
        self.ab.append((a, b))
        return True

    def congruence(self, S):
        N = self.N
        E1: Dict[int, int] = {}
        E2: Dict[int, int] = {}
        k = 0
        i = 0
        while S:
            if S & 1:
                k += 1
                for p, e in self.f1[i].items():
                    E1[p] = E1.get(p, 0) + e
                for p, e in self.f2[i].items():
                    E2[p] = E2.get(p, 0) + e
            S >>= 1
            i += 1
        assert k % 2 == 0
        U = 1
        for p, e in E1.items():
            U = U * pow(p, e // 2, N) % N
        W = 1
        for p, e in E2.items():
            W = W * pow(p, e // 2, N) % N
        c1, d1 = self.L1[0], self.L2[0]
        X = pow(d1 % N, k // 2, N) * U % N
        Y = pow(c1 % N, k // 2, N) * W % N
        return X, Y


def _sieve_box(sysm, A, Bb, cap=None):
    tried = 0
    for b in range(1, Bb + 1):
        for a in range(-A, A + 1):
            if gcd(a, b) != 1:
                continue
            tried += 1
            sysm.try_pair(a, b)
            if cap is not None and len(sysm.rows) >= cap(sysm):
                return tried
    return tried


def lin_base(N, seed, B, A, Bb, excess=None):
    m = isqrt(N)
    C1, C0 = divmod(N, m)
    L1, L2 = (1, -m), (C1, C0)
    sysm = Linear(N, "lin_base", dict(B=B, A=A, Bb=Bb), L1, L2, B)
    for c in (L1[0], L2[0]):
        sysm._inv_check(c)
    cap = None if excess is None else (lambda s_: s_.ncols + excess)
    sysm.info = dict(m=m, L1=L1, L2=L2, tried=_sieve_box(sysm, A, Bb, cap))
    return sysm


def lin_rand(N, seed, B, A, Bb, excess=None):
    rng = _rng("lin_rand", N, seed)
    while True:
        m = rng.randrange(2, N - 1)
        u, v = _gauss_reduce((1, -m), (0, N))
        t = rng.randrange(-2, 3)
        L1 = u
        L2 = (v[0] + t * u[0], v[1] + t * u[1])
        if L1[0] % N and L2[0] % N and L1[0] * L2[1] - L1[1] * L2[0] != 0:
            break
    sysm = Linear(N, "lin_rand", dict(B=B, A=A, Bb=Bb), L1, L2, B)
    for c in (L1[0], L2[0]):
        sysm._inv_check(c)
    cap = None if excess is None else (lambda s_: s_.ncols + excess)
    sysm.info = dict(m=m, L1=L1, L2=L2, t=t, tried=_sieve_box(sysm, A, Bb, cap))
    return sysm


# ---------------------------------------------------------------- monic cubic NFS
class Cubic(System):
    def __init__(self, N, family, params, f, m, B, nchar, rng):
        super().__init__(N, family, params)
        self.f = f                      # [f0, f1, f2, 1]
        self.m = m
        self.O = Order(f)
        self.fb = FactorBase(B)
        self.rcol = {p: 1 + i for i, p in enumerate(self.fb.primes)}   # col 0 = rational sign
        self.acol: Dict[Tuple[int, int], int] = {}
        self.base_cols = 1 + len(self.fb.primes)
        self.rel_r: List[Dict[int, int]] = []
        self.ab: List[Tuple[int, int]] = []
        self.arows: List[Dict[Tuple[int, int], int]] = []
        self.nchar = nchar
        self.rng = rng
        self.rationals: List[int] = []
        self.akeys = set()

    def F(self, a, b):
        f0, f1, f2, _ = self.f
        return a ** 3 + f2 * a * a * b + f1 * a * b * b + f0 * b ** 3

    def ncols_now(self):
        return self.base_cols + len(self.akeys) + self.nchar

    def try_pair(self, a, b):
        z = a - b * self.m
        if z == 0:
            return False
        fr = self.fb.factor(z)
        if fr is None:
            return False
        Fv = self.F(a, b)
        if Fv == 0:
            return False
        fa = self.fb.factor(Fv)
        if fa is None:
            return False
        arow = {}
        for r, e in fa.items():
            s = a * pow(b, -1, r) % r       # b coprime to r since f monic and gcd(a,b)=1
            arow[(r, s)] = e
        self.rel_r.append(fr)
        self.ab.append((a, b))
        self.arows.append(arow)
        self.akeys.update(arow)
        self.rows.append(0)          # placeholder; real rows built in finalize()
        self.rationals.append(z)
        return True

    def finalize(self, inert_start=1 << 40):
        """Assign algebraic columns, choose quadratic characters, build parity rows."""
        keys = sorted({k for ar in self.arows for k in ar})
        self.acol = {k: self.base_cols + i for i, k in enumerate(keys)}
        self.char_base = self.base_cols + len(keys)
        self.chars = self._pick_chars(self.nchar)
        self.ncols = self.char_base + len(self.chars)
        self.rows = []
        for i, (a, b) in enumerate(self.ab):
            r = 1 if self.rationals[i] < 0 else 0
            for p, e in self.rel_r[i].items():
                if e & 1:
                    r |= 1 << self.rcol[p]
            for k, e in self.arows[i].items():
                if e & 1:
                    r |= 1 << self.acol[k]
            for j, (q, s) in enumerate(self.chars):
                if gmpy2.legendre((a - b * s) % q, q) == -1:
                    r |= 1 << (self.char_base + j)
            self.rows.append(r)
        self.inert = self.O.inert_primes(inert_start, 3)

    def _pick_chars(self, k):
        out = []
        q = int(gmpy2.next_prime(1 << 26))
        disc = int(fmpz_poly(self.f).discriminant()) if hasattr(fmpz_poly(self.f), "discriminant") else 0
        while len(out) < k:
            q = int(gmpy2.next_prime(q + self.rng.randrange(1, 1 << 12)))
            if disc and disc % q == 0:
                continue
            fac = nmod_poly(self.f, q).factor()[1]
            roots = [(-int(g.coeffs()[0])) % q for g, e in fac if g.degree() == 1 and e == 1]
            if not roots:
                continue
            s = self.rng.choice(roots)
            if any((a - b * s) % q == 0 for a, b in self.ab):
                continue
            out.append((q, s))
        return out

    def congruence(self, S) -> Optional[Tuple[int, int]]:
        N, O = self.N, self.O
        elems = []
        E: Dict[int, int] = {}
        sgn = 1
        i = 0
        T = S
        while T:
            if T & 1:
                a, b = self.ab[i]
                elems.append(O.linear(a, b))
                for p, e in self.rel_r[i].items():
                    E[p] = E.get(p, 0) + e
                if self.rationals[i] < 0:
                    sgn = -sgn
            T >>= 1
            i += 1
        assert sgn == 1
        dp = O.deriv_at_alpha()
        g = O.mul(O.mul(dp, dp), O.prod(elems))
        v = None
        for l in self.inert:
            try:
                v = O.sqrt(g, l)
                break
            except ValueError:
                continue
        if v is None:
            return None
        fpm = 3 * self.m * self.m + 2 * self.f[2] * self.m + self.f[1]
        U = fpm % N
        for p, e in E.items():
            assert e % 2 == 0
            U = U * pow(p, e // 2, N) % N
        Y = O.eval_mod(v, self.m, N)
        return U, Y


def _cubic_poly(N, m, R):
    """monic base-m cubic + (x - m)(r1 x + r0); requires m^3 <= N < 2 m^3."""
    c = []
    t = N
    for _ in range(3):
        t, r = divmod(t, m)
        c.append(r)
    assert t == 1, "leading digit must be 1"
    r1, r0 = R
    # (x - m)(r1 x + r0) = r1 x^2 + (r0 - m r1) x - m r0
    f = [c[0] - m * r0, c[1] + r0 - m * r1, c[2] + r1, 1]
    return f


def nfs3(N, seed, B, A, Bb, Hmul, nchar, family, excess=None, Hshift=0):
    rng = _rng(family, N, seed)
    lo = gmpy2.iroot(gmpy2.mpz(N) // 2, 3)[0] + 1
    hi = gmpy2.iroot(gmpy2.mpz(N), 3)[0]
    tries = 0
    while True:
        tries += 1
        m = int(rng.randrange(int(lo), int(hi) + 1))
        if not (m ** 3 <= N < 2 * m ** 3):
            continue
        H = (Hmul * m) >> Hshift
        R = (rng.randrange(-H, H + 1), rng.randrange(-H, H + 1)) if H > 0 else (0, 0)
        f = _cubic_poly(N, m, R)
        fac = fmpz_poly(f).factor()[1]
        if len(fac) == 1 and fac[0][1] == 1:
            break
        for g, _e in fac:          # reducible f: try to extract a factor (record if so)
            if g.degree() == 1:
                root = -int(g.coeffs()[0]) // int(g.coeffs()[1]) if int(g.coeffs()[1]) else 0
                gg = gcd(m - root, N)
                if 1 < gg < N:
                    raise FactorFound(gg)
    sysm = Cubic(N, family, dict(B=B, A=A, Bb=Bb, Hmul=Hmul, Hshift=Hshift, nchar=nchar), f, m, B, nchar, rng)
    sysm.info = dict(m=m, f=f, R=R, poly_tries=tries)
    cap = None if excess is None else (lambda s_: s_.ncols_now() + excess)
    sysm.info["tried"] = _sieve_box(sysm, A, Bb, cap)
    sysm.finalize()
    return sysm


# ---------------------------------------------------------------- square-multiplier linear sieve
class SquareMultiplierLinear(Linear):
    """Degree-1 sieve whose multiplier lambda = c1/d1 (mod N) is r^2 for a uniformly random
    unit r.  Everything except the final Y depends on lambda only (hidden-root property).
    No |S|-parity column: odd-weight dependencies are allowed and use Y = r^{|S|} W."""

    def __init__(self, N, params, L1, L2, B, r):
        super().__init__(N, "lin_sq", params, L1, L2, B)
        self._r = r
        self.ncols = self.par            # drop the parity column (index 2P+2)

    def try_pair(self, a, b):
        ok = super().try_pair(a, b)
        if ok:
            self.rows[-1] &= ~(1 << self.par)
        return ok

    def congruence(self, S, r_override=None):
        N = self.N
        E1: Dict[int, int] = {}
        E2: Dict[int, int] = {}
        k = 0
        i = 0
        while S:
            if S & 1:
                k += 1
                for p, e in self.f1[i].items():
                    E1[p] = E1.get(p, 0) + e
                for p, e in self.f2[i].items():
                    E2[p] = E2.get(p, 0) + e
            S >>= 1
            i += 1
        U = 1
        for p, e in E1.items():
            U = U * pow(p, e // 2, N) % N
        W = 1
        for p, e in E2.items():
            W = W * pow(p, e // 2, N) % N
        r = self._r if r_override is None else r_override
        return U, pow(r, k, N) * W % N


def lin_sq(N, seed, B, A, Bb, excess=None):
    rng = _rng("lin_sq", N, seed)
    while True:
        r = rng.randrange(2, N - 1)
        g = gcd(r, N)
        if g != 1:
            raise FactorFound(g)              # ACCIDENTAL
        lam = r * r % N
        e1, e2 = _gauss_reduce((lam, 1), (N, 0))   # lattice {(c,d): c = lam d (mod N)}
        t = rng.randrange(-2, 3)
        c1, d1 = e1
        c0, d0 = e2[0] + t * e1[0], e2[1] + t * e1[1]
        if c1 % N and d1 % N and c1 * d0 - c0 * d1 != 0:
            break
    for c in (c1, d1):
        g = gcd(c, N)
        if 1 < g < N:
            raise FactorFound(g)
    L1, L2 = (c1, c0), (d1, d0)
    sysm = SquareMultiplierLinear(N, dict(B=B, A=A, Bb=Bb), L1, L2, B, r)
    cap = None if excess is None else (lambda s_: s_.ncols + excess)
    sysm.info = dict(L1=L1, L2=L2, t=t, tried=_sieve_box(sysm, A, Bb, cap))
    return sysm
