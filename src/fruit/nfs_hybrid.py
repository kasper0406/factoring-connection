"""OPERATIONAL: E2 toy -- monic cubic base-m NFS relations (no hidden randomness) plus h rows
obtained by Dixon descent of uniform squares INTO THE SAME FIELD AND FACTOR BASE.

Descent row j (all quantities except r_j are functions of R_j = r_j^2 mod N and public data):
  (u, v) = shortest vector of {u = R v (mod N)};  u = s_u u_s prod_{l in L_u} l^e,
  v = s_v v_s prod_{l in L_v} l^e  (u_s, v_s B-smooth, B < l <= B1);
  special-q in the cubic field: a_l - b_l m = l * c_l, c_l and F(a_l, b_l) B-smooth.
  With P_u = prod_{L_u} (a_l - b_l alpha)^e, C_u = prod_{L_u} c_l^e (same for v):
     X_j := r_j v_s C_u phi(P_v),  Z_j := s_u s_v u_s v_s C_u C_v,  Pi_j := P_u P_v,
     X_j^2 = Z_j phi(Pi_j)  (mod N).
  Its parity row is (rational parities of Z_j | ideal parities of Pi_j | characters of Pi_j).
For S = I (base pairs) + J (descent rows) with prod_I z_i prod_J Z_j = U0^2 and
f'(alpha)^2 prod_I (a_i - b_i alpha) prod_J Pi_j = V^2 in Z[alpha]:
     X = f'(m) U0 prod_J X_j,   Y = phi(V) prod_J Z_j,   X^2 = Y^2 (mod N).
No inversion mod N is needed.  Uses nfsfruit's Cubic/Order (imported, not edited).
"""
from __future__ import annotations

import random
from math import gcd
from typing import Dict, List, Optional, Tuple

import gmpy2

from ..nfsfruit import families as fam
from ..nfsfruit.smooth import FactorBase
from ..lemmaM.descent import gauss_short, factor_cofactor

FactorFound = fam.FactorFound


def _gauss2(e1, e2):
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]
    a, b = e1, e2
    if n2(a) < n2(b):
        a, b = b, a
    while True:
        dot = a[0] * b[0] + a[1] * b[1]
        nb = n2(b)
        k = (2 * dot + nb) // (2 * nb)
        a = (a[0] - k * b[0], a[1] - k * b[1])
        if n2(a) >= nb:
            return b, a
        a, b = b, a


class CubicHybrid(fam.Cubic):
    def __init__(self, N, params, f, m, B, nchar, rng, B1, J, kmax):
        super().__init__(N, "nfs3_hybrid", params, f, m, B, nchar, rng)
        self.B1, self.J, self.kmax = B1, J, kmax
        self.desc: List[Dict] = []          # descent rows
        self.sq_cache: Dict[int, Optional[Tuple]] = {}
        self.counters = dict(r_trials=0, split_smooth=0, sq_attempts=0, sq_success=0,
                             sq_points=0, desc_accepted=0)

    # ------------------------------------------------------------ special-q in the cubic field
    def special_q(self, ell: int):
        if ell in self.sq_cache:
            return self.sq_cache[ell]
        self.counters["sq_attempts"] += 1
        e1, e2 = _gauss2((ell, 0), (self.m % ell, 1))
        pts = []
        J = self.J
        for i in range(-J, J + 1):
            for j in range(-J, J + 1):
                a = i * e1[0] + j * e2[0]
                b = i * e1[1] + j * e2[1]
                if b <= 0:
                    continue
                pts.append((i * i + j * j, a, b))
        pts.sort()
        res = None
        for _, a, b in pts:
            self.counters["sq_points"] += 1
            if gcd(a, b) != 1:
                continue
            z = a - b * self.m
            if z == 0 or z % ell:
                continue
            c = z // ell
            fc = self.fb.factor(c)
            if fc is None:
                continue
            Fv = self.F(a, b)
            if Fv == 0:
                continue
            fa = self.fb.factor(Fv)
            if fa is None:
                continue
            arow = {}
            for rr, e in fa.items():
                arow[(rr, a * pow(b, -1, rr) % rr)] = e
            res = (a, b, c, fc, arow)
            break
        if res is not None:
            self.counters["sq_success"] += 1
        self.sq_cache[ell] = res
        return res

    def _split(self, z: int):
        """(sign, small B-smooth factorisation, large primes dict) or None."""
        za = abs(z)
        g = gmpy2.gcd(za, self.fb.P)
        small = {}
        r = gmpy2.mpz(za)
        while g > 1:
            r //= g
            g = gmpy2.gcd(r, g)
        sm = za // int(r)
        if sm > 1:
            small = self.fb.factor(sm)
        r = int(r)
        large = factor_cofactor(r, self.B1, self.kmax) if r > 1 else {}
        if large is None:
            return None
        return (-1 if z < 0 else 1), small, large

    def split_and_descend(self, R: int):
        """A function of R and public data only."""
        N = self.N
        u, v = gauss_short(R, N)
        if u == 0 or v == 0:
            return None
        for t in (u, v):
            g = gcd(t, N)
            if 1 < g < N:
                raise FactorFound(g)
        su = self._split(u)
        sv = self._split(v)
        if su is None or sv is None:
            return None
        self.counters["split_smooth"] += 1
        kids = []
        for side, (_s, _sm, large) in (("u", su), ("v", sv)):
            for ell, e in large.items():
                sq = self.special_q(ell)
                if sq is None:
                    return None
                kids.append((side, ell, e, sq))
        return dict(R=R, su=su, sv=sv, kids=kids)

    def add_descent(self, rng: random.Random, max_trials: int) -> bool:
        N = self.N
        for _ in range(max_trials):
            self.counters["r_trials"] += 1
            r = rng.randrange(1, N)
            g = gcd(r, N)
            if g != 1:
                raise FactorFound(g)
            d = self.split_and_descend(r * r % N)
            if d is None:
                continue
            d["r"] = r                      # hidden coin: used only in congruence()
            self.desc.append(d)
            self.counters["desc_accepted"] += 1
            return True
        return False

    # ------------------------------------------------------------ rows
    def _desc_parts(self, d):
        """Z_j factorisation (signed), Pi_j as list of (a, b, e), C_u, v_s, P_v pieces."""
        (s_u, us, _), (s_v, vs, _) = d["su"], d["sv"]
        sign = s_u * s_v
        E: Dict[int, int] = {}
        for fdict in (us, vs):
            for p, e in fdict.items():
                E[p] = E.get(p, 0) + e
        pis = []
        for side, ell, e, (a, b, c, fc, arow) in d["kids"]:
            if c < 0 and e % 2:
                sign = -sign
            for p, ee in fc.items():
                E[p] = E.get(p, 0) + ee * e
            pis.append((side, a, b, e, c, arow))       # c signed
        return sign, E, pis

    def finalize_hybrid(self):
        """Columns/characters over base pairs AND descent children; build all rows."""
        kid_ab = [(a, b) for d in self.desc for (_, _, _, (a, b, _c, _f, _r)) in d["kids"]]
        kid_ar = [ar for d in self.desc for (_, _, _, (_a, _b, _c, _f, ar)) in d["kids"]]
        keys = sorted({k for ar in list(self.arows) + kid_ar for k in ar})
        self.acol = {k: self.base_cols + i for i, k in enumerate(keys)}
        self.char_base = self.base_cols + len(keys)
        saved = list(self.ab)
        self.ab = saved + kid_ab            # so _pick_chars avoids every pair used
        self.chars = self._pick_chars(self.nchar)
        self.ab = saved
        self.ncols = self.char_base + len(self.chars)
        self.rows = []
        for i, (a, b) in enumerate(self.ab):
            self.rows.append(self._pair_row(self.rel_r[i], self.rationals[i] < 0, self.arows[i], [(a, b, 1)]))
        self.drows = []
        for d in self.desc:
            sign, E, pis = self._desc_parts(d)
            arow: Dict = {}
            for (_side, a, b, e, _c, ar) in pis:
                for k, ee in ar.items():
                    arow[k] = arow.get(k, 0) + ee * e
            self.drows.append(self._pair_row(E, sign < 0, arow, [(a, b, e) for (_s, a, b, e, _c, _a) in pis]))
        self.inert = self.O.inert_primes(1 << 40, 3)

    def _pair_row(self, rfac, neg, arow, lin):
        row = 1 if neg else 0
        for p, e in rfac.items():
            if e & 1:
                row |= 1 << self.rcol[p]
        for k, e in arow.items():
            if e & 1:
                row |= 1 << self.acol[k]
        for j, (q, s) in enumerate(self.chars):
            par = 0
            for (a, b, e) in lin:
                if e & 1 and gmpy2.legendre((a - b * s) % q, q) == -1:
                    par ^= 1
            if par:
                row |= 1 << (self.char_base + j)
        return row

    # ------------------------------------------------------------ congruence
    def hybrid_congruence(self, base_idx: List[int], desc_idx: List[int]):
        """Returns (X, Y) or None if the algebraic product is not a square."""
        N, O, m = self.N, self.O, self.m
        elems = []
        E: Dict[int, int] = {}
        sgn = 1
        for i in base_idx:
            a, b = self.ab[i]
            elems.append(O.linear(a, b))
            for p, e in self.rel_r[i].items():
                E[p] = E.get(p, 0) + e
            if self.rationals[i] < 0:
                sgn = -sgn
        Xd = 1
        Zprod = 1
        for j in desc_idx:
            d = self.desc[j]
            sign, Ej, pis = self._desc_parts(d)
            sgn *= sign
            for p, e in Ej.items():
                E[p] = E.get(p, 0) + e
            Zj = sign
            for p, e in Ej.items():
                Zj *= p ** e
            Zprod = Zprod * Zj % N
            # X_j = r v_s C_u phi(P_v)
            (_su, _us, _), (_sv, vs, _) = d["su"], d["sv"]
            vs_int = 1
            for p, e in vs.items():
                vs_int *= p ** e
            Xj = d["r"] * vs_int % N
            for (side, a, b, e, c, _ar) in pis:
                for _ in range(e):
                    elems.append(O.linear(a, b))
                if side == "u":
                    Xj = Xj * pow(c % N, e, N) % N
                else:
                    Xj = Xj * pow((a - b * m) % N, e, N) % N
            Xd = Xd * Xj % N
        if sgn != 1:
            return None
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
        fpm = 3 * m * m + 2 * self.f[2] * m + self.f[1]
        U = fpm % N
        for p, e in E.items():
            assert e % 2 == 0
            U = U * pow(p, e // 2, N) % N
        X = U * Xd % N
        Y = O.eval_mod(v, m, N) * Zprod % N
        return X, Y


def build_base(N: int, seed: int, B: int, A: int, Bb: int, nchar: int, excess: int,
               B1: int, J: int, kmax: int) -> CubicHybrid:
    """Monic base-m cubic (R = 0), as nfsfruit nfs3_base, with descent support."""
    rng = random.Random(f"fruit-e2|{N}|{seed}")
    lo = gmpy2.iroot(gmpy2.mpz(N) // 2, 3)[0] + 1
    hi = gmpy2.iroot(gmpy2.mpz(N), 3)[0]
    while True:
        m = int(rng.randrange(int(lo), int(hi) + 1))
        if not (m ** 3 <= N < 2 * m ** 3):
            continue
        f = fam._cubic_poly(N, m, (0, 0))
        from flint import fmpz_poly
        fac = fmpz_poly(f).factor()[1]
        if len(fac) == 1 and fac[0][1] == 1:
            break
        for g_, _e in fac:
            if g_.degree() == 1 and int(g_.coeffs()[1]):
                root = -int(g_.coeffs()[0]) // int(g_.coeffs()[1])
                gg = gcd(m - root, N)
                if 1 < gg < N:
                    raise FactorFound(gg)
    sysm = CubicHybrid(N, dict(B=B, A=A, Bb=Bb, nchar=nchar), f, m, B, nchar, rng, B1, J, kmax)
    sysm.info = dict(m=m, f=f)
    cap = lambda s_: s_.ncols_now() + excess
    sysm.info["tried"] = fam._sieve_box(sysm, A, Bb, cap)
    return sysm
