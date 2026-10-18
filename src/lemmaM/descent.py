"""OPERATIONAL: Dixon relations by descent of uniformly random squares (toy, QS special-q).

The candidate sees only N, public parameters and a seed.  It never factors N.
Smoothness tests on auxiliary integers (|u|, |v| <= sqrt(N), QS cofactors) use exact
trial division by the factor base plus Pollard rho / BPSW (gmpy2.is_prime) on the
cofactor; this is the smoothness test, documented as such (it never touches N).
Modular square roots are taken modulo factor-base-size primes ell only.

Descent relation (one per accepted r):
  r uniform in [2, N-2], R = r^2 mod N.   Only R is passed to `split_and_descend`.
  (u, v) = shortest nonzero vector of {(u, v): u = R v (mod N)}  (Gauss reduction).
  Accept iff u, v are B1-smooth and every prime ell in (B, B1] dividing uv has (N/ell) = 1
  and admits a QS special-q relation  x_ell^2 - N = ell * c_ell,  c_ell B-smooth, found by
  scanning x = x0 + t, t = t_ell (mod ell), 0 <= t < J*ell, x0 = isqrt(N) + 1.
  Then  r^2 v = u  and  ell = x_ell^2 / c_ell  (mod N)  give
      (r * W)^2 = T  (mod N),   W = prod_{ell|v} x_ell^e / prod_{ell|u} x_ell^e,
      T = (u_small * prod_{ell|v} c_ell^e) / (v_small * prod_{ell|u} c_ell^e)  in Q,
  with T a B-smooth rational number.  The relation stores (r, W, exponent vector of T).
  Everything except the factor r is a function of R (and the public seed-free search).

Plain QS relations (x, x^2 - N) can be mixed in.  A dependency S has
  X = prod X_i,  Y = prod_p p^{E_p/2} (signed, inverses mod N),  X^2 = Y^2 (mod N).
Theorem DD (notes/lemmaM_results.md): if S contains a descent relation, X/Y is uniform on
mu_2(Z/N) conditional on everything except the hidden roots, so Pr[fruitful] = 1 - 2^{1-w}.
"""
from __future__ import annotations

import random
from math import gcd, isqrt
from typing import Dict, List, Optional, Tuple

import gmpy2

from ..nfsfruit.smooth import FactorBase


class FactorFound(Exception):
    def __init__(self, g):
        super().__init__(g)
        self.g = int(g)


# ------------------------------------------------------------------ small helpers
def gauss_short(R: int, N: int) -> Tuple[int, int]:
    """Shortest nonzero vector (u, v) of the lattice {u = R v (mod N)}, v >= 0."""
    a, b = (N, 0), (R % N, 1)
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]
    if n2(a) < n2(b):
        a, b = b, a
    while True:
        # b is the shorter one
        dot = a[0] * b[0] + a[1] * b[1]
        nb = n2(b)
        k = (2 * dot + nb) // (2 * nb)   # round(dot/nb)
        a = (a[0] - k * b[0], a[1] - k * b[1])
        if n2(a) >= nb:
            break
        a, b = b, a
    u, v = b
    if v < 0 or (v == 0 and u < 0):
        u, v = -u, -v
    return u, v


def tonelli(n: int, p: int) -> int:
    """Square root of n modulo an odd prime p (n a QR mod p)."""
    n %= p
    if n == 0:
        return 0
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


def _rho(n: int, seed: int) -> int:
    n = gmpy2.mpz(n)
    c = gmpy2.mpz(seed % (n - 1) + 1)
    x = y = gmpy2.mpz(2)
    d = gmpy2.mpz(1)
    while d == 1:
        x = (x * x + c) % n
        y = (y * y + c) % n
        y = (y * y + c) % n
        d = gmpy2.gcd(abs(x - y), n)
    return int(d)


def factor_cofactor(c: int, B1: int, max_primes: int) -> Optional[Dict[int, int]]:
    """Exact factorisation of c (> 0) into primes <= B1 with at most max_primes prime
    factors (with multiplicity), else None.  Smoothness test on an auxiliary integer."""
    out: Dict[int, int] = {}
    stack = [c]
    count = 0
    while stack:
        z = stack.pop()
        if z == 1:
            continue
        if gmpy2.is_prime(z):
            if z > B1:
                return None
            out[z] = out.get(z, 0) + 1
            count += 1
            if count > max_primes:
                return None
            continue
        if z > B1 ** max_primes:
            return None
        s = 1
        d = z
        while d in (1, z):
            d = _rho(z, s)
            s += 1
        stack += [d, z // d]
    return out


def jacobi_ok(N: int, ell: int) -> bool:
    return pow(N % ell, (ell - 1) // 2, ell) == 1


# ------------------------------------------------------------------ the system
class DescentSystem:
    """Rows over columns: 0 = sign, 1.. = primes <= B."""

    def __init__(self, N: int, B: int, B1: int, J: int, max_large: int = 3):
        self.N = N
        self.B, self.B1, self.J, self.max_large = B, B1, J, max_large
        self.fb = FactorBase(B)
        self.col = {p: i + 1 for i, p in enumerate(self.fb.primes)}
        self.ncols = len(self.fb.primes) + 1
        self.x0 = isqrt(N) + 1
        self.rows: List[int] = []
        self.X: List[int] = []          # X_i mod N (includes hidden r for descent rels)
        self.E: List[Dict[int, int]] = []   # signed exponent vector of T_i (key -1 = sign)
        self.kind: List[str] = []       # 'D' (descent) or 'Q' (plain QS)
        self.hidden_r: List[Optional[int]] = []  # the algorithm's own coin r (D) / None (Q)
        self.sq_cache: Dict[int, Optional[Tuple[int, Dict[int, int]]]] = {}
        self.counters = dict(r_trials=0, split_smooth=0, split_rejected_qr=0,
                             sq_attempts=0, sq_success=0, sq_tests=0, desc_accepted=0,
                             qs_tests=0, accidental=0)

    def _check(self, z):
        g = gcd(z % self.N, self.N)
        if 1 < g < self.N:
            raise FactorFound(g)

    # -- factor-base factorisation of a signed integer; None if not B-smooth
    def _fbfac(self, z: int) -> Optional[Dict[int, int]]:
        if z == 0:
            return None
        f = self.fb.factor(abs(z))
        if f is None:
            return None
        f = dict(f)
        if z < 0:
            f[-1] = 1
        return f

    def _row(self, E: Dict[int, int]) -> int:
        r = 0
        for p, e in E.items():
            if e & 1:
                r |= 1 << (0 if p == -1 else self.col[p])
        return r

    def special_q(self, ell: int) -> Optional[Tuple[int, Dict[int, int]]]:
        """QS special-q relation x^2 - N = ell * c with c B-smooth (deterministic in ell)."""
        if ell in self.sq_cache:
            return self.sq_cache[ell]
        self.counters["sq_attempts"] += 1
        N, x0 = self.N, self.x0
        s = tonelli(N, ell)
        res = None
        for root in sorted({s, (-s) % ell}):
            t = (root - x0) % ell
            for j in range(self.J):
                x = x0 + t + j * ell
                Q = x * x - N
                self.counters["sq_tests"] += 1
                c = Q // ell
                f = self._fbfac(c)
                if f is not None:
                    res = (x, f)
                    break
            if res is not None:
                break
        if res is not None:
            self.counters["sq_success"] += 1
        self.sq_cache[ell] = res
        return res

    def split_and_descend(self, R: int) -> Optional[Tuple[int, Dict[int, int]]]:
        """Function of R only.  Returns (W mod N, signed exponent vector of T) or None."""
        N = self.N
        u, v = gauss_short(R, N)
        if v == 0 or u == 0:
            return None
        self._check(v)
        self._check(u)
        parts = []
        for z in (u, v):
            za = abs(z)
            small = {}
            for p in self.fb.primes:
                if za % p == 0:
                    e = 0
                    while za % p == 0:
                        za //= p
                        e += 1
                    small[p] = e
                if p * p > za and za > 1 and za <= self.B:
                    break
            if za > 1 and za <= self.B:      # leftover is a factor-base prime
                small[za] = small.get(za, 0) + 1
                za = 1
            large = factor_cofactor(za, self.B1, self.max_large) if za > 1 else {}
            if large is None:
                return None
            if z < 0:
                small[-1] = 1
            parts.append((small, large))
        self.counters["split_smooth"] += 1
        for _, large in parts:
            for ell in large:
                if ell <= self.B or not jacobi_ok(N, ell):
                    self.counters["split_rejected_qr"] += 1
                    return None
        E: Dict[int, int] = {}
        def addE(f, sgn):
            for p, e in f.items():
                E[p] = E.get(p, 0) + sgn * e
        W = 1
        (us, ul), (vs, vl) = parts
        addE(us, +1)
        addE(vs, -1)
        for large, side in ((ul, +1), (vl, -1)):
            for ell, e in large.items():
                sq = self.special_q(ell)
                if sq is None:
                    return None
                x, c = sq
                self._check(x)
                # ell = x^2 / c ; u-side: T gets c^{-e}, W gets x^{-e}; v-side opposite
                addE(c, -side * e)
                W = W * pow(x, -side * e, N) % N
        E = {p: e for p, e in E.items() if e != 0}
        if -1 in E:
            E[-1] %= 2
            if E[-1] == 0:
                del E[-1]
        return W, E

    def add_descent(self, rng: random.Random, max_trials: int) -> bool:
        N = self.N
        for _ in range(max_trials):
            self.counters["r_trials"] += 1
            r = rng.randrange(2, N - 1)
            g = gcd(r, N)
            if g != 1:
                self.counters["accidental"] += 1
                raise FactorFound(g)
            R = r * r % N
            out = self.split_and_descend(R)
            if out is None:
                continue
            W, E = out
            self.X.append(r * W % N)
            self.E.append(E)
            self.rows.append(self._row(E))
            self.kind.append("D")
            self.hidden_r.append(r)
            self.counters["desc_accepted"] += 1
            return True
        return False

    def add_qs(self, M: int, limit: int):
        N, x0 = self.N, self.x0
        added = 0
        for t in range(M):
            if added >= limit:
                break
            x = x0 + t
            Q = x * x - N
            self.counters["qs_tests"] += 1
            f = self._fbfac(Q)
            if f is None:
                continue
            self.X.append(x % N)
            self.E.append(f)
            self.rows.append(self._row(f))
            self.kind.append("Q")
            self.hidden_r.append(None)
            added += 1
        return added

    def congruence(self, S: int) -> Tuple[int, int]:
        N = self.N
        X = 1
        E: Dict[int, int] = {}
        i = 0
        while S:
            if S & 1:
                X = X * self.X[i] % N
                for p, e in self.E[i].items():
                    E[p] = E.get(p, 0) + e
            S >>= 1
            i += 1
        Y = 1
        for p, e in E.items():
            assert e % 2 == 0, "not a square"
            if p == -1:
                if (e // 2) % 2:
                    Y = (-Y) % N
                continue
            Y = Y * pow(p, e // 2, N) % N
        return X, Y


def build(N: int, seed: int, B: int, B1: int, J: int, n_desc: int, qs_M: int = 0,
          qs_limit: int = 0, max_trials: int = 200000, max_large: int = 3) -> DescentSystem:
    rng = random.Random(f"lemmaM-desc|{N}|{seed}|{B}|{B1}|{J}")
    sysm = DescentSystem(N, B, B1, J, max_large)
    if qs_limit:
        sysm.add_qs(qs_M, qs_limit)
    for _ in range(n_desc):
        if not sysm.add_descent(rng, max_trials):
            break
    return sysm
