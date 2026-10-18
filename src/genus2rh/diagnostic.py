"""DIAGNOSTIC (factor-aware) helpers for genus2rh. Never imported by operational code.

- jac_order(p, F): #J(F_p) for y^2 = F(X), F monic quintic, from exact point counts over
  F_p and F_{p^2} (L-polynomial 1 + c1 T + c2 T^2 + p c1 T^3 + p^2 T^4, #J = L(1)).
- squarefree_mod(p, F): F mod p has nonzero discriminant.
- order_of(p, F, D, n): exact order of D in J(F_p) given a multiple n of it (sympy factorint).
"""
import numpy as np
from sympy import factorint

from .cantor import Jac
from .hecm import resultant_F_dF

__all__ = ["jac_order", "squarefree_mod", "order_of", "largest_prime"]


def _legendre_table(p):
    t = np.full(p, -1, dtype=np.int64)
    t[0] = 0
    sq = (np.arange(1, p, dtype=np.int64) ** 2) % p
    t[sq] = 1
    return t


def _nonresidue(p, chi):
    for d in range(2, p):
        if chi[d] == -1:
            return d
    raise ValueError


def point_counts(p, F):
    """(#C(F_p), #C(F_{p^2})) for the smooth model with one point at infinity."""
    F = [c % p for c in F]
    chi = _legendre_table(p)
    x = np.arange(p, dtype=np.int64)
    val = np.zeros(p, dtype=np.int64)
    for c in reversed(F):
        val = (val * x + c) % p
    N1 = 1 + p + int(chi[val].sum())
    d = _nonresidue(p, chi)
    # elements re + im*s with s^2 = d
    re0 = np.repeat(np.arange(p, dtype=np.int64), p)
    im0 = np.tile(np.arange(p, dtype=np.int64), p)
    vr = np.zeros(p * p, dtype=np.int64)
    vi = np.zeros(p * p, dtype=np.int64)
    for c in reversed(F):
        nr = (vr * re0 + d * ((vi * im0) % p) + c) % p
        ni = (vr * im0 + vi * re0) % p
        vr, vi = nr, ni
    norm = (vr * vr - d * ((vi * vi) % p)) % p
    N2 = 1 + p * p + int(chi[norm].sum())
    return N1, N2


def jac_order(p, F):
    N1, N2 = point_counts(p, F)
    c1 = N1 - p - 1
    two_c2 = N2 - p * p - 1 + c1 * c1
    assert two_c2 % 2 == 0
    c2 = two_c2 // 2
    return 1 + c1 + c2 + p * c1 + p * p


def squarefree_mod(p, F):
    return resultant_F_dF([c % p for c in F]) % p != 0


def order_of(p, F, D, n):
    J = Jac(p, F)
    assert Jac.is_identity(J.smul(n, D))
    o = n
    for ell, e in factorint(n).items():
        for _ in range(e):
            if Jac.is_identity(J.smul(o // ell, D)):
                o //= ell
            else:
                break
    return o


def largest_prime(n):
    return max(factorint(n)) if n > 1 else 1
