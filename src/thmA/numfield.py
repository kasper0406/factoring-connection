"""Design (a'''_k; mu) fields and exact arithmetic in Z[beta] (agent thmA, independent implementation).
OPERATIONAL: receives N and public parameters only.

Conventions (notes/thmA_proof.md §1):
  F_hat = Lemma GB applied to T = N*w at the point (g0, g1);  R_c(x, y) = sum_{i<d} c_i x^(d-1-i) y^i;
  F_A = F_hat + G * R_c with G(x, y) = g1 x - g0 y;  f(x) = F_A(x, 1), f[i] = coefficient of x^i;
  F(a, b) = sum_i f[i] a^i b^(d-i);  t = f_d;  beta = t*alpha with monic g(y) = t^(d-1) f(y/t);
  z_P = t*a - b*beta in Z[beta];  psi: beta -> t*g0/g1 (mod N);  psi(z_P) = t G(P)/g1.
Polynomial factorisation over F_p (python-flint nmod_poly.factor / roots) is a standard subroutine
(Berlekamp / Cantor-Zassenhaus); it is never applied to N.
"""
from __future__ import annotations

import random
from math import gcd, isqrt
from typing import Dict, List, Optional, Sequence, Tuple

import flint

from .arith import is_prime


# ------------------------------------------------------------------ Lemma GB_w
def _centered(x: int, m: int) -> int:
    r = x % m
    return r - m if 2 * r > m else r


def gbe(T: int, d: int, g0: int, g1: int) -> List[int]:
    """a[0..d] with sum_i a_i g0^i g1^(d-i) = T, |a_i| <= g1/2 (i >= 1). Needs gcd(g0, g1) = 1, g1 >= 2."""
    assert gcd(g0, g1) == 1 and g1 >= 2
    a = [0] * (d + 1)
    for i in range(d, 0, -1):
        a[i] = _centered(T * pow(g0, -i, g1), g1)
        T -= a[i] * g0 ** i
        assert T % g1 == 0
        T //= g1
    a[0] = T
    return a


def hom_eval(f: Sequence[int], a: int, b: int) -> int:
    d = len(f) - 1
    return sum(f[i] * a ** i * b ** (d - i) for i in range(d + 1))


def make_f(N: int, w: int, d: int, g0: int, g1: int, c: Sequence[int]) -> List[int]:
    """f = F_hat_{g,w}(x, 1) + (g1 x - g0) * sum_i c_i x^(d-1-i)."""
    f = gbe(N * w, d, g0, g1)
    for i, ci in enumerate(c):
        j = d - 1 - i                    # R contributes c_i x^j
        f[j + 1] += g1 * ci
        f[j] -= g0 * ci
    assert hom_eval(f, g0, g1) == N * w
    return f


def disc(f: Sequence[int]) -> int:
    """Discriminant of f (coefficients low->high), exact (Sylvester resultant via flint)."""
    d = len(f) - 1
    if d == 3:
        e, c, b, a = f
        return b * b * c * c - 4 * a * c ** 3 - 4 * b ** 3 * e - 27 * a * a * e * e + 18 * a * b * c * e
    P = flint.fmpz_poly(list(f))
    dP = P.derivative()
    res = int(P.resultant(dP))
    sign = -1 if (d * (d - 1) // 2) % 2 else 1
    assert res % f[d] == 0
    return sign * res // f[d]


def roots_mod(f: Sequence[int], r: int) -> List[int]:
    """Simple roots of f mod prime r (Las Vegas root finding)."""
    P = flint.nmod_poly([x % r for x in f], r)
    if P.degree() <= 0:
        return []
    return sorted(int(z) for z, m in P.roots() if m == 1)


def irreducible_mod(f: Sequence[int], r: int) -> bool:
    """f mod r irreducible of full degree d (certifies irreducibility over Q when r does not divide f_d)."""
    d = len(f) - 1
    if f[d] % r == 0:
        return False
    P = flint.nmod_poly([x % r for x in f], r)
    _, facs = P.factor()
    return len(facs) == 1 and facs[0][1] == 1 and facs[0][0].degree() == d


# ------------------------------------------------------------------ Z[beta]
def monic_g(f: Sequence[int]) -> List[int]:
    d = len(f) - 1
    t = f[d]
    return [f[i] * t ** (d - 1 - i) for i in range(d)] + [1]


def reduce_g(c: List[int], g: Sequence[int], m: int | None = None) -> List[int]:
    d = len(g) - 1
    c = list(c)
    for i in range(len(c) - 1, d - 1, -1):
        ti = c[i]
        if ti:
            for j in range(d):
                c[i - d + j] -= ti * g[j]
        c[i] = 0
    c = (c + [0] * d)[:d]
    if m is not None:
        c = [x % m for x in c]
    return c


def zmul(x: Sequence[int], y: Sequence[int], g: Sequence[int], m: int | None = None) -> List[int]:
    d = len(g) - 1
    pr = [0] * (2 * d - 1)
    for i, xi in enumerate(x):
        if xi:
            for j, yj in enumerate(y):
                if yj:
                    pr[i + j] += xi * yj
    if m is not None:
        pr = [p % m for p in pr]
    return reduce_g(pr, g, m)


def zprod(elems: List[List[int]], g: Sequence[int]) -> List[int]:
    d = len(g) - 1
    level = [list(e) for e in elems] or [[1] + [0] * (d - 1)]
    while len(level) > 1:
        nxt = [zmul(level[i], level[i + 1], g) for i in range(0, len(level) - 1, 2)]
        if len(level) % 2:
            nxt.append(level[-1])
        level = nxt
    return level[0]


def z_atom(a: int, b: int, t: int, d: int) -> List[int]:
    e = [0] * d
    e[0] = t * a
    e[1] = -b
    return e


def gprime(g: Sequence[int]) -> List[int]:
    d = len(g) - 1
    return [(i + 1) * g[i + 1] for i in range(d)]


def psi(elem: Sequence[int], beta_img: int, N: int) -> int:
    r = 0
    for c in reversed(elem):
        r = (r * beta_img + c) % N
    return r


def sqrt_bound(g: Sequence[int], atom_l1: List[Tuple[int, int]]) -> int:
    """Lemma SQ(ii) bound A on ||eta||_inf for eta^2 = g'(beta)^2 * prod z^{e} (atom_l1 = [(||z||_1, e)]).
    Integer arithmetic only; log2 A is computed through bit lengths to avoid huge powers."""
    d = len(g) - 1
    Mh = isqrt(sum(c * c for c in g)) + 1
    gl1 = sum(abs(c) for c in g)
    gd = d * gl1 * Mh ** (d - 1)
    U = gd * gd
    for n1, e in atom_l1:
        U *= (n1 * Mh ** (d - 1)) ** e
    return d * (isqrt(U) + 1) * 2 ** (d - 1) * Mh * gd ** (d - 1)


# ------------------------------------------------------------------ F_p[y]/phi arithmetic (toy sizes)
def _pmod(c: List[int], phi: List[int], p: int) -> List[int]:
    """c mod phi over F_p; phi monic (low->high)."""
    c = [x % p for x in c]
    k = len(phi) - 1
    for i in range(len(c) - 1, k - 1, -1):
        ti = c[i]
        if ti:
            for j in range(k):
                c[i - k + j] = (c[i - k + j] - ti * phi[j]) % p
        c[i] = 0
    return (c + [0] * k)[:k]


def _pmul(x, y, phi, p):
    pr = [0] * (len(x) + len(y) - 1)
    for i, xi in enumerate(x):
        if xi:
            for j, yj in enumerate(y):
                pr[i + j] += xi * yj
    return _pmod(pr, phi, p)


def _ppow(x, e, phi, p):
    k = len(phi) - 1
    r = [1] + [0] * (k - 1)
    while e:
        if e & 1:
            r = _pmul(r, x, phi, p)
        x = _pmul(x, x, phi, p)
        e >>= 1
    return r


def fq_sqrt(c: List[int], phi: List[int], p: int, rng: random.Random) -> Optional[List[int]]:
    """Square root in F_q = F_p[y]/(phi), q = p^k, by Tonelli-Shanks; None for a non-square."""
    k = len(phi) - 1
    q = p ** k
    one = [1] + [0] * (k - 1)
    c = _pmod(c, phi, p)
    if not any(c):
        return c
    if _ppow(c, (q - 1) // 2, phi, p) != one:
        return None
    s, t = 0, q - 1
    while t % 2 == 0:
        s += 1
        t //= 2
    minus1 = [p - 1] + [0] * (k - 1)
    while True:
        z = [rng.randrange(p) for _ in range(k)]
        if any(z) and _ppow(z, (q - 1) // 2, phi, p) == minus1:
            break
    M, cc, tt, R = s, _ppow(z, t, phi, p), _ppow(c, t, phi, p), _ppow(c, (t + 1) // 2, phi, p)
    while tt != one:
        i, t2 = 0, tt
        while t2 != one:
            t2 = _pmul(t2, t2, phi, p)
            i += 1
        b = cc
        for _ in range(M - i - 1):
            b = _pmul(b, b, phi, p)
        M, cc = i, _pmul(b, b, phi, p)
        tt, R = _pmul(tt, cc, phi, p), _pmul(R, b, phi, p)
    return R


def hensel_prime(g: Sequence[int], D: int, lo: int) -> int:
    """Least odd prime p > lo with p not dividing D and g mod p squarefree."""
    p = lo + 1
    while True:
        if p > 2 and is_prime(p) and D % p != 0:
            P = flint.nmod_poly([x % p for x in g], p)
            if P.gcd(P.derivative()).degree() == 0:
                return p
        p += 1


def zsqrt(Gamma: List[int], g: Sequence[int], A: int, p: int, rng: random.Random,
          stats: Dict | None = None) -> Optional[List[int]]:
    """Algorithm SQ (notes/lv6_results.md §3.5), independent implementation: eta in Z[beta] with
    eta^2 == Gamma exactly, or None.  Hensel lifting of the inverse square root at ONE prime p (g mod p
    squarefree, Gamma a unit mod p), 2^(tau-1) sign patterns, exact verification."""
    d = len(g) - 1
    P = flint.nmod_poly([x % p for x in g], p)
    _, facs = P.factor()
    phis = [[int(c) for c in phi.coeffs()] for phi, m in facs]
    assert all(m == 1 for _, m in facs)
    gb = [x % p for x in g]
    roots = []
    for phi in phis:
        cr = _pmod(list(Gamma), phi, p)
        if not any(cr):
            raise ValueError("Gamma not a unit at p")
        rr = fq_sqrt(cr, phi, p, rng)
        if rr is None:
            if stats is not None:
                stats["residue_nonsquare"] = stats.get("residue_nonsquare", 0) + 1
            return None
        roots.append(rr)
    # CRT idempotents in F_p[y]/(g)
    gP = flint.nmod_poly(gb, p)
    idem = []
    for phi in phis:
        ph = flint.nmod_poly(phi, p)
        Mj = gP // ph
        inv = flint.nmod_poly([int(c) for c in Mj.coeffs()], p)
        # inverse of Mj modulo phi via xgcd
        gg, s_, _t = inv.xgcd(ph)
        assert gg.degree() == 0
        s_ = s_ * pow(int(gg.coeffs()[0]), -1, p)
        e = (Mj * s_) % gP
        idem.append([int(c) for c in e.coeffs()])
    e_prec = 1
    while p ** e_prec <= 2 * A:
        e_prec *= 2
    for mask in range(2 ** len(phis)):
        if mask & 1:
            continue
        rho = [0] * d
        for j in range(len(phis)):
            sgn = -1 if (mask >> j) & 1 else 1
            term = _poly_mul_mod(idem[j], roots[j], gb, p)
            rho = [(x + sgn * y) % p for x, y in zip(rho, term)]
        rP = flint.nmod_poly(rho, p)
        gg, s_, _t = rP.xgcd(gP)
        if gg.degree() != 0:
            continue
        s_ = s_ * pow(int(gg.coeffs()[0]), -1, p)
        z = ([int(c) for c in (s_ % gP).coeffs()] + [0] * d)[:d]
        prec = 1
        while prec < e_prec:
            prec = min(2 * prec, e_prec)
            m = p ** prec
            Gm = [c % m for c in Gamma]
            Gz2 = zmul(zmul(z, z, g, m), Gm, g, m)
            corr = [(-c) % m for c in Gz2]
            corr[0] = (corr[0] + 1) % m
            inv2 = pow(2, -1, m)
            dz = zmul(z, corr, g, m)
            z = [(zi + inv2 * di) % m for zi, di in zip(z, dz)]
        m = p ** e_prec
        eta = zmul([c % m for c in Gamma], z, g, m)
        eta = [c - m if 2 * c > m else c for c in eta]
        if zmul(eta, eta, g) == list(Gamma):
            if stats is not None:
                stats["tau"] = len(phis)
            return eta
    return None


def _mul_raw(x, y):
    pr = [0] * (len(x) + len(y) - 1)
    for i, xi in enumerate(x):
        for j, yj in enumerate(y):
            pr[i + j] += xi * yj
    return pr


def _poly_mul_mod(x, y, gmod, p):
    """x*y mod (gmod, p) where gmod is monic of degree d; inputs low->high (y may be shorter)."""
    return _pmod(_mul_raw(x, y), list(gmod), p)
