"""Algorithm SQ (notes/lv6_results.md §3): decide whether Gamma in Z[beta] is a square in K and, if so,
return its square root eta in Z[beta], by Hensel lifting at ONE unramified prime p with an enumeration of
the 2^t sign choices (t = number of prime factors of g mod p, t <= d; d is fixed).

Certification: a returned eta always satisfies eta^2 == Gamma exactly (checked).  If Gamma = eta^2 with
eta in Z[beta] and ||eta||_inf <= A, the algorithm finds it (Lemma SQ).  Operational: public data only.
"""
from __future__ import annotations

import random
from typing import List, Optional, Sequence, Tuple

import flint
import sympy

from . import numfield as nf


def _fq_pow(x, e, phi):
    return x.pow_mod(e, phi)


def _fq_sqrt(c, phi, p: int, rng: random.Random):
    """Square root of c in F_p[y]/(phi) (phi irreducible of degree k), or None if c is a non-square.
    Tonelli-Shanks in F_q, q = p^k (randomized non-residue search: Las Vegas)."""
    ctx = phi.context()
    k = phi.degree()
    q = p ** k
    c = c % phi
    if c.is_zero():
        return c
    if _fq_pow(c, (q - 1) // 2, phi) != ctx([1]):
        return None
    s, t = 0, q - 1
    while t % 2 == 0:
        s += 1
        t //= 2
    while True:
        z = ctx([rng.randrange(p) for _ in range(k)]) % phi
        if not z.is_zero() and _fq_pow(z, (q - 1) // 2, phi) == ctx([p - 1]):
            break
    M = s
    cc = _fq_pow(z, t, phi)
    tt = _fq_pow(c, t, phi)
    R = _fq_pow(c, (t + 1) // 2, phi)
    one = ctx([1])
    while tt != one:
        i, t2 = 0, tt
        while t2 != one:
            t2 = (t2 * t2) % phi
            i += 1
        bpow = cc
        for _ in range(M - i - 1):
            bpow = (bpow * bpow) % phi
        M = i
        cc = (bpow * bpow) % phi
        tt = (tt * cc) % phi
        R = (R * bpow) % phi
    assert (R * R - c) % phi == ctx([0]) or ((R * R) % phi) == c
    return R


def choose_prime(g: Sequence[int], lower: int, avoid: int) -> int:
    """Least prime p > lower with p odd, p not dividing avoid (e.g. disc(g) * f_d), g mod p squarefree."""
    p = max(3, lower + 1)
    while True:
        p = int(sympy.nextprime(p - 1)) if not sympy.isprime(p) else p
        if p % 2 == 1 and avoid % p != 0:
            ctx = flint.fmpz_mod_poly_ctx(p)
            if ctx([c % p for c in g]).is_squarefree():
                return p
        p += 1


def hensel_sqrt(Gamma: Sequence[int], g: Sequence[int], A: int, p: int,
                rng: random.Random | None = None, stats: dict | None = None) -> Optional[List[int]]:
    """Return eta in Z[beta] with eta^2 = Gamma, or None (then Gamma is not the square of any eta in Z[beta]
    with ||eta||_inf <= A).  Requires p odd, g mod p squarefree, Gamma a unit mod p."""
    rng = rng or random.Random(0)
    d = len(g) - 1
    ctx = flint.fmpz_mod_poly_ctx(p)
    gbar = ctx([c % p for c in g])
    Gbar = ctx([c % p for c in Gamma])
    _, facs = gbar.factor()
    phis = [phi for phi, mult in facs]
    assert all(m == 1 for _, m in facs)
    roots = []
    for phi in phis:
        if (Gbar % phi).is_zero():
            raise ValueError("Gamma is not a unit at p")
        r = _fq_sqrt(Gbar, phi, p, rng)
        if r is None:
            if stats is not None:
                stats["residue_nonsquare"] = stats.get("residue_nonsquare", 0) + 1
            return None
        roots.append(r)
    # CRT idempotents e_j with e_j = 1 mod phi_j, 0 mod phi_i (i != j)
    idem = []
    for j, phi in enumerate(phis):
        Mj = gbar.exact_division(phi)
        inv = Mj.inverse_mod(phi)
        idem.append((Mj * inv) % gbar)
    # precision: p^e > 2A
    e = 1
    while p ** e <= 2 * A:
        e *= 2
    tries = 0
    for mask in range(2 ** len(phis)):
        if mask & 1:      # eta and -eta give the same pair: fix the sign on the first factor
            continue
        tries += 1
        rho = ctx([0])
        for j in range(len(phis)):
            s = -1 if (mask >> j) & 1 else 1
            rho = (rho + idem[j] * roots[j] * s) % gbar
        rho_inv = rho.inverse_mod(gbar)
        z = [int(c) for c in rho_inv.coeffs()]
        z = (z + [0] * d)[:d]
        # Newton for the inverse square root: z <- z + z (1 - Gamma z^2) / 2  (mod p^prec)
        prec = 1
        while prec < e:
            prec = min(2 * prec, e)
            m = p ** prec
            Gz2 = nf.mul(nf.mul(z, z, g, m), [c % m for c in Gamma], g, m)
            corr = [(-c) % m for c in Gz2]
            corr[0] = (corr[0] + 1) % m
            inv2 = pow(2, -1, m)
            dz = nf.mul(z, corr, g, m)
            z = [(zi + inv2 * di) % m for zi, di in zip(z, dz)]
        m = p ** e
        eta = nf.mul([c % m for c in Gamma], z, g, m)
        eta = [c - m if 2 * c > m else c for c in eta]
        if nf.mul(eta, eta, g) == list(Gamma):
            if stats is not None:
                stats["sign_tries"] = stats.get("sign_tries", 0) + tries
                stats["t"] = len(phis)
                stats["prec_bits"] = (p ** e).bit_length()
            return eta
    if stats is not None:
        stats["sign_tries"] = stats.get("sign_tries", 0) + tries
        stats["t"] = len(phis)
    return None
