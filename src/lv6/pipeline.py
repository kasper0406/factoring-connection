"""Toy end-to-end LV6 pipeline (notes/lv6_results.md; pre-registration notes/lv6_prereg.md).

OPERATIONAL part (receives N, public config, seed): toy design-(a''_k) field, ordinary two-sided relations
(a, b) with G = g1 a - g0 b and F(a, b) B-smooth, "rows" = disjoint signed groups of 1..gmax relations
(mimicking Dixon-descent rows u/v), the LV6 column set (sign, rational primes, ONE parity column,
first-degree algebraic primes (r, s) with r not dividing D = disc(f) f_d, quadratic characters at large
first-degree primes), the left kernel over F2, Algorithm SQ (src/lv6/sqrt.py) and the congruence X^2 = Y^2.

No factor of N is used anywhere.  Smoothness is decided exactly by repeated gcd with the primorial
(Bernstein-style smooth-part computation, one number at a time).  Every gcd with N is checked; a nonunit
is recorded as an ACCIDENTAL factor hit (random gcd).
"""
from __future__ import annotations

import random
from math import gcd
from typing import Dict, List, Sequence, Tuple

import sympy
from sympy import integer_nthroot

from src.pcest.geometry import field as pc_field, hom_eval

from . import chars, numfield as nf
from .sqrt import choose_prime, hensel_sqrt


# ----------------------------------------------------------------------------- smoothness (exact)
def smooth_factor(x: int, primorial: int) -> Dict[int, int] | None:
    """Factorisation of |x| if it is B-smooth (primorial = prod of primes <= B), else None."""
    x = abs(x)
    if x == 0:
        return None
    rad = gcd(x, primorial)
    y = x
    while True:
        gg = gcd(y, primorial)
        if gg == 1:
            break
        y //= gg
    if y != 1:
        return None
    fac: Dict[int, int] = {}
    for p in sympy.factorint(rad):           # rad is squarefree, all prime factors <= B
        e = 0
        while x % p == 0:
            x //= p
            e += 1
        fac[p] = e
    return fac


# ----------------------------------------------------------------------------- field
def toy_field(N: int, d: int, hbits: int, rng: random.Random, events: List[str]):
    M = integer_nthroot(N, d + 1)[0]
    H = 2 ** hbits
    while True:
        g1 = rng.randrange(max(2, M // 2), M + 1)
        g0 = rng.randrange(max(1, M // 2), M + 1)
        if gcd(g0, g1) != 1:
            continue
        R = [rng.randrange(-H, H + 1) for _ in range(d)]
        f = pc_field(N, d, g0, g1, R)
        fd = f[d]
        if fd == 0:
            events.append("degree_drop")
            continue
        df = nf.disc_poly(f)
        if df == 0:
            events.append("disc_zero")
            continue
        gg = gcd(fd * g1 * df, N)
        if gg != 1:
            events.append("gcd_fd_g1_disc_nonunit" + ("_proper" if gg != N else "_N"))
            continue
        x = sympy.Symbol("x")
        if not sympy.Poly(list(reversed(f)), x).is_irreducible:
            events.append("reducible")
            continue
        return {"g0": g0, "g1": g1, "R": R, "f": f, "disc_f": df}


# ----------------------------------------------------------------------------- relations
def collect_relations(N: int, fld: dict, B: int, A: int, max_atoms: int, events: List[str]):
    f, g0, g1 = fld["f"], fld["g0"], fld["g1"]
    primorial = 1
    for p in sympy.primerange(2, B + 1):
        primorial *= p
    atoms = []
    for b in range(1, A + 1):
        for a in range(-A, A + 1):
            if gcd(a, b) != 1:
                continue
            F = hom_eval(f, a, b)
            if F == 0:
                continue
            fF = smooth_factor(F, primorial)
            if fF is None:
                continue
            G = g1 * a - g0 * b
            if G == 0:
                continue
            fG = smooth_factor(G, primorial)
            if fG is None:
                continue
            gg = gcd(F * G, N)
            if gg != 1:
                events.append("ACCIDENTAL_gcd_FG")
                continue
            atoms.append({"a": a, "b": b, "F": F, "G": G, "fF": fF, "fG": fG})
            if len(atoms) >= max_atoms:
                return atoms
    return atoms


def atom_columns(atom: dict, D: int, fd: int) -> Dict[tuple, int]:
    """Integer exponent vector of one relation over the LV6 columns (before reduction mod 2)."""
    cols: Dict[tuple, int] = {("par",): 1}
    if atom["G"] < 0:
        cols[("sgn",)] = 1
    for p, e in atom["fG"].items():
        cols[("q", p)] = e
    a, b = atom["a"], atom["b"]
    for r, e in atom["fF"].items():
        if D % r == 0:
            continue                       # bad prime: no column (absorbed in H, Lemma H)
        assert b % r != 0                  # r | b would force r | f_d
        s = a * pow(b, -1, r) % r
        cols[("a", r, s)] = e
    return cols


# ----------------------------------------------------------------------------- F2 linear algebra
def left_kernel(rows: List[int], ncols: int) -> List[int]:
    """Basis of {S subset rows : sum of rows in S = 0 over F2}, as bitsets over row indices
    (Gaussian elimination with history bits).  Deterministic given the row order."""
    n = len(rows)
    work = [(rows[i], 1 << i) for i in range(n)]
    pivots: Dict[int, Tuple[int, int]] = {}
    kernel = []
    for vec, hist in work:
        while vec:
            top = vec.bit_length() - 1
            if top in pivots:
                pv, ph = pivots[top]
                vec ^= pv
                hist ^= ph
            else:
                pivots[top] = (vec, hist)
                break
        if vec == 0:
            kernel.append(hist)
    return kernel
