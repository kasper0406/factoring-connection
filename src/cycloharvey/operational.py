"""OPERATIONAL code: Harvey-style one-sided collision search (Harvey 2021
Alg. 4.1/4.2, randomised base instead of Hittmeir's Prop. 2.7), and the
cyclotomic-handle variant (guess r = p mod m, handles via sigma_r in
R = (Z/N)[X]/(Phi_m)).

Receives only N, public parameters and a seed.  Never uses p, q or group
orders.  Every inversion mod N is gcd-checked; every gcd hit is returned.
flint usage: nmod_poly multiplication and remainder by monic polynomials
only (src/common/README.md audit); no gcd/xgcd/factor on composite moduli.
"""
from __future__ import annotations

import random
from math import gcd, isqrt
from typing import Dict, List, Optional, Sequence, Tuple

import flint

from src.common.metrics import Metrics, NullMetrics
from .calibration import ceil_sqrt_ratio


class FactorFound(Exception):
    def __init__(self, d: int, how: str):
        self.d, self.how = d, how


def _check(x: int, N: int, met, how: str):
    met.inc("gcd")
    g = gcd(x % N, N)
    if 1 < g < N:
        raise FactorFound(g, how)
    return g


def inv_mod(a: int, N: int, met) -> int:
    g = _check(a, N, met, "inverse")
    if g != 1:
        raise ZeroDivisionError("non-invertible (gcd = N)")
    met.inc("inv")
    return pow(a, -1, N)


# --------------------------------------------------------------------------
# Lemma 3.1 (recover p, q from u = a q + b p): roots of y^2 - u y + abN
# --------------------------------------------------------------------------
def recover_from_u(N: int, a: int, b: int, u: int) -> Optional[int]:
    disc = u * u - 4 * a * b * N
    if disc < 0:
        return None
    s = isqrt(disc)
    if s * s != disc:
        return None
    for y in ((u + s) // 2, (u - s) // 2):
        if (u + s) % 2:
            return None
        for d in (gcd(y, N),):
            if 1 < d < N:
                return d
    return None


# --------------------------------------------------------------------------
# Algorithm 4.1: collisions modulo an unknown factor
# --------------------------------------------------------------------------
def product_tree(vals: List[int], N: int, met) -> flint.nmod_poly:
    polys = [flint.nmod_poly([(-v) % N, 1], N) for v in vals]
    while len(polys) > 1:
        nxt = []
        for i in range(0, len(polys) - 1, 2):
            met.inc("poly_mul")
            met.inc("poly_mul_dimsum", polys[i].degree() + polys[i + 1].degree())
            nxt.append(polys[i] * polys[i + 1])
        if len(polys) % 2:
            nxt.append(polys[-1])
        polys = nxt
    return polys[0]


def bluestein_eval(f: flint.nmod_poly, alpha: int, kappa: int, N: int, met) -> List[int]:
    """f(alpha^i) for i = 0..kappa-1 via ij = C(i,2) + C(-j,2) - C(i-j,2)."""
    n = f.degree()
    fc = [int(c) for c in f.coeffs()] + [0] * (n + 1 - len(f.coeffs()))
    ainv = inv_mod(alpha, N, met)
    # h_i = alpha^{C(i,2)}, i >= 0 ; f'_j = alpha^{C(j+1,2)} f_j ; g_k = alpha^{-C(k,2)}
    def powers_binom(base: int, ks: range) -> Dict[int, int]:
        out = {}
        for k in ks:
            e = k * (k - 1) // 2
            out[k] = pow(base, e, N)
            met.inc("rn_pow")
        return out
    h = powers_binom(alpha, range(0, kappa))
    fpj = [(pow(alpha, j * (j + 1) // 2, N) * fc[j]) % N for j in range(n + 1)]
    met.inc("rn_pow", n + 1)
    gk = powers_binom(ainv, range(-n, kappa))
    F = flint.nmod_poly(fpj, N)
    G = flint.nmod_poly([gk[k] for k in range(-n, kappa)], N)
    met.inc("poly_mul")
    met.inc("poly_mul_dimsum", n + kappa + n)
    P = F * G
    pc = [int(c) for c in P.coeffs()]
    out = []
    for i in range(kappa):
        idx = i + n
        s = pc[idx] if idx < len(pc) else 0
        out.append((h[i] * s) % N)
    met.inc("rn_mul", 2 * kappa)
    return out


def find_collisions(N: int, alpha: int, kappa: int, vals: List[int], met) -> Optional[int]:
    """Harvey Alg. 4.1.  Precondition: no v equals alpha^i mod N."""
    if not vals:
        return None
    f = product_tree(vals, N, met)
    ev = bluestein_eval(f, alpha, kappa, N, met)
    ai = 1
    for i in range(kappa):
        met.inc("gcd")
        g = gcd(ev[i], N)
        if g > 1:
            # locate the colliding giant (Harvey Step 3c; also used to
            # classify the hit as structural/accidental in the harness)
            for h, v in enumerate(vals):
                met.inc("gcd_locate")
                g2 = gcd((v - ai) % N, N)
                if 1 < g2 < N:
                    met.set("hit", [i, h])
                    return g2
        ai = ai * alpha % N
    return None


# --------------------------------------------------------------------------
# Harvey Alg. 4.2 over (Z/N)^*  (F_p^* handles)
# --------------------------------------------------------------------------
def harvey_targets(N: int, r: int, m: int):
    """(a, b, shift, W) with shift = ceil(sqrt(4abN)), W = #giant steps."""
    for a in range(1, r + 1):
        for b in range(1, r // a + 1):
            shift = ceil_sqrt_ratio(4 * a * b * N, 1)
            J = max(1, ceil_sqrt_ratio(N, 16 * r * r * m * m * a * b))
            yield a, b, shift, J


def harvey_search(N: int, r: int, m: int, seed: int = 0, met=None,
                  early_exit: bool = True) -> Dict:
    """early_exit=False builds the complete giant list (cost = model count)
    before collision detection; exact matches are still verified."""
    met = met or Metrics()
    found = None
    rng = random.Random(f"cycloharvey-harvey|{N}|{seed}")
    try:
        alpha = rng.randrange(2, N - 1)
        _check(alpha, N, met, "base")
        baby, x = {}, 1
        for i in range(m):
            if i:
                _check(x - 1, N, met, "baby_order")
            if x in baby:
                return {"factor": None, "status": "small_order", "met": met.to_record()}
            baby[x] = i
            x = x * alpha % N
            met.inc("rn_mul")
        am_inv = inv_mod(pow(alpha, m, N), N, met)
        giants: List[int] = []
        meta: List[Tuple[int, int, int, int]] = []
        for a, b, shift, J in harvey_targets(N, r, m):
            met.inc("handles")
            e = a * N + b - shift
            t = pow(alpha, e, N)
            met.inc("rn_pow")
            # square-and-multiply count, comparable with Cyclo.pow's R_mul
            met.inc("rn_mul_pow", e.bit_length() + bin(e).count("1") - 2)
            v = t
            for j in range(J):
                if v in baby:                      # exact match mod N
                    met.inc("exact_match")
                    u = baby[v] + j * m + shift
                    d = recover_from_u(N, a, b, u)
                    if d and early_exit:
                        return {"factor": d, "status": "exact_match", "met": met.to_record()}
                    if d and found is None:
                        found = d
                else:
                    giants.append(v)
                    meta.append((a, b, j, shift))
                v = v * am_inv % N
                met.inc("rn_mul")
        met.set("baby", m)
        met.set("giant", len(giants) + met.counts.get("exact_match", 0))
        if found is not None:
            return {"factor": found, "status": "exact_match", "met": met.to_record()}
        d = find_collisions(N, alpha, m, giants, met)
        hit = met.fields.get("hit") if hasattr(met, "fields") else None
        info = None
        if d and hit:
            i, h = hit
            info = {"baby": i, "target": list(meta[h])}
        return {"factor": d, "status": "collision" if d else "none", "hit": info,
                "met": met.to_record()}
    except FactorFound as e:
        return {"factor": e.d, "status": "gcd_" + e.how, "met": met.to_record()}


# --------------------------------------------------------------------------
# cyclotomic ring R = (Z/N)[X]/(Phi_m) and Frobenius handles
# --------------------------------------------------------------------------
class Cyclo:
    def __init__(self, N: int, m: int, met=None):
        self.N, self.m = N, m
        self.phi_poly = flint.fmpz_poly.cyclotomic(m)
        self.D = self.phi_poly.degree()
        self.mod = flint.nmod_poly([int(c) for c in self.phi_poly.coeffs()], N)
        self.met = met or NullMetrics()
        self.units = [s for s in range(1, m) if gcd(s, m) == 1]

    def elt(self, coeffs: Sequence[int]) -> flint.nmod_poly:
        return flint.nmod_poly([c % self.N for c in coeffs], self.N) % self.mod

    def mul(self, a, b):
        self.met.inc("R_mul")
        self.met.inc("R_mul_D", self.D)
        return (a * b) % self.mod

    def pow(self, a, e: int):
        assert e >= 0
        res = self.elt([1])
        base = a
        while e:
            if e & 1:
                res = self.mul(res, base)
            e >>= 1
            if e:
                base = self.mul(base, base)
        return res

    def sigma(self, a, s: int):
        """X -> X^s (s a unit mod m); exact via X^m = 1."""
        self.met.inc("R_sigma")
        c = [int(x) for x in a.coeffs()]
        out = [0] * self.m
        for i, ci in enumerate(c):
            out[(i * s) % self.m] = (out[(i * s) % self.m] + ci) % self.N
        return flint.nmod_poly(out, self.N) % self.mod

    def norm(self, a) -> int:
        prod = self.elt([1])
        for s in self.units:
            prod = self.mul(prod, self.sigma(a, s))
        c = [int(x) for x in prod.coeffs()]
        assert all(x == 0 for x in c[1:]), "norm not in base ring"
        return c[0] if c else 0

    def inv(self, a):
        """a^{-1} = prod_{s != 1} sigma_s(a) / Norm(a), Norm gcd-checked."""
        nm = self.norm(a)
        ninv = inv_mod(nm, self.N, self.met)
        prod = self.elt([ninv])
        for s in self.units:
            if s != 1:
                prod = self.mul(prod, self.sigma(a, s))
        return prod

    def coeffs(self, a) -> List[int]:
        c = [int(x) for x in a.coeffs()]
        return c + [0] * (self.D - len(c))

    def check_unit_coeffs(self, a, how: str):
        """gcd of every coefficient with N (zero test mod a hidden factor)."""
        for c in self.coeffs(a):
            _check(c, self.N, self.met, how)


def handle(R: Cyclo, alpha, r: int, F: Dict[Tuple[int, int], int], alpha_inv=None):
    """T_{r,F} = prod sigma_r^{i-j}(alpha)^{c_ij N^j}  (== alpha^{F(p,q)} mod p
    when r == p mod m).  sigma_r^k for negative k uses the inverse unit r^{-1}."""
    m, N = R.m, R.N
    rinv = pow(r, -1, m)
    out = R.elt([1])
    for (i, j), c in F.items():
        k = i - j
        s = pow(r, k, m) if k >= 0 else pow(rinv, -k, m)
        base = R.sigma(alpha, s)
        e = c * N ** j
        if e < 0:
            base = R.inv(base)
            e = -e
        out = R.mul(out, R.pow(base, e))
    return out


# --------------------------------------------------------------------------
# general multipoint evaluation (remainder tree; monic divisors only)
# --------------------------------------------------------------------------
def multipoint_eval(f: flint.nmod_poly, pts: List[int], N: int, met) -> List[int]:
    if not pts:
        return []
    levels = [[flint.nmod_poly([(-x) % N, 1], N) for x in pts]]
    while len(levels[-1]) > 1:
        cur, nxt = levels[-1], []
        for i in range(0, len(cur) - 1, 2):
            met.inc("poly_mul")
            nxt.append(cur[i] * cur[i + 1])
        if len(cur) % 2:
            nxt.append(cur[-1])
        levels.append(nxt)
    rems = [f % levels[-1][0]]
    for lev in range(len(levels) - 2, -1, -1):
        cur, new = levels[lev], []
        for i, node in enumerate(cur):
            met.inc("poly_rem")
            new.append(rems[i // 2] % node)
        rems = new
    return [int(r.coeffs()[0]) if r.coeffs() else 0 for r in rems]


def cyclo_harvey_search(N: int, m: int, r: int, mb: int, seed: int = 0, met=None,
                        guesses: Sequence[int] = None) -> Dict:
    """Harvey-style search with cyclotomic Frobenius handles.  For every guess
    s in (Z/m)^* (or the given subset), targets u = a q + b p - shift with
    handle T = sigma_s^{-1}(alpha)^{aN} sigma_s(alpha)^{b} alpha^{-shift};
    giant steps T alpha^{-j mb}; collision detection on a random Z/N-linear
    projection lambda (alpha^i == v mod p  =>  lambda equal mod p).
    Cost is reported in R-multiplications (each = D-dimensional)."""
    met = met or Metrics()
    rng = random.Random(f"cycloharvey-cyclo|{N}|{m}|{seed}")
    try:
        if gcd(m, N) != 1:
            _check(m, N, met, "conductor")
        R = Cyclo(N, m, met)
        alpha = R.elt([rng.randrange(N) for _ in range(R.D)])
        ainv = R.inv(alpha)                       # gcd-checked via the norm
        lam = [rng.randrange(1, N) for _ in range(R.D)]
        proj = lambda x: sum(c * l for c, l in zip(R.coeffs(x), lam)) % N
        baby_full, baby_proj, x = {}, [], R.elt([1])
        for i in range(mb):
            key = tuple(R.coeffs(x))
            baby_full[key] = i
            baby_proj.append(proj(x))
            x = R.mul(x, alpha)
        amb_inv = R.pow(ainv, mb)
        giants, gmeta = [], []
        guesses = list(guesses) if guesses is not None else R.units
        for s in guesses:
            sp = R.sigma(alpha, s)
            sq = R.sigma(alpha, pow(s, -1, m))
            for a, b, shift, J in harvey_targets(N, r, mb):
                met.inc("handles")
                T = R.mul(R.mul(R.pow(sq, a * N), R.pow(sp, b)), R.pow(ainv, shift))
                v = T
                for j in range(J):
                    key = tuple(R.coeffs(v))
                    if key in baby_full:
                        met.inc("exact_match")
                        u = baby_full[key] + j * mb + shift
                        d = recover_from_u(N, a, b, u)
                        if d:
                            return {"factor": d, "status": "exact_match", "met": met.to_record()}
                    else:
                        giants.append(proj(v))
                        gmeta.append((s, a, b, j, shift))
                    v = R.mul(v, amb_inv)
        met.set("giant", len(giants))
        met.set("baby", mb)
        met.set("D", R.D)
        f = product_tree(giants, N, met)
        vals = multipoint_eval(f, baby_proj, N, met)
        for i, val in enumerate(vals):
            met.inc("gcd")
            g = gcd(val, N)
            if g > 1:
                for h, gv in enumerate(giants):
                    met.inc("gcd_locate")
                    g2 = gcd((gv - baby_proj[i]) % N, N)
                    if 1 < g2 < N:
                        return {"factor": g2, "status": "collision",
                                "hit": {"baby": i, "target": list(gmeta[h])},
                                "met": met.to_record()}
        return {"factor": None, "status": "none", "met": met.to_record()}
    except FactorFound as e:
        return {"factor": e.d, "status": "gcd_" + e.how, "met": met.to_record()}
