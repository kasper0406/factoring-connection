"""Norm-compression (Lemma NC / ptecons LEDGER) check on exact presentations (DIAGNOSTIC).

Presentation: monic F in Z[X] of degree d, small height; N := |F(m)| for a chosen m (SNFS-type, so that
N | F(m) exactly).  Members: multisets A of size M from [-B, B] with ||P_A mod F||_inf <= h.
For each member we record the number of non-unit pieces (|F(a)| >= 2) and check
    prod_j |F(a_j)| = |N(r(theta))| <= ||r||_1^d * ||F||_2^(d-1)   (F irreducible).
"""
from __future__ import annotations

import math
from itertools import combinations_with_replacement

from flint import fmpz_poly


def poly_from_roots(roots) -> fmpz_poly:
    p = fmpz_poly([1])
    for a in roots:
        p *= fmpz_poly([-a, 1])
    return p


def members(F: fmpz_poly, M: int, B: int, h: int):
    out = []
    for A in combinations_with_replacement(range(-B, B + 1), M):
        r = poly_from_roots(A) % F
        cs = [int(c) for c in r.coeffs()]
        if cs and max(abs(c) for c in cs) <= h and any(cs):
            out.append((A, cs))
    return out


def nc_row(F: fmpz_poly, M: int, B: int, h: int) -> dict:
    d = F.degree()
    fl2 = math.sqrt(sum(int(c) ** 2 for c in F.coeffs()))
    irreducible = len(F.factor()[1]) == 1 and F.factor()[1][0][1] == 1  # factor of F itself, not of N
    mem = members(F, M, B, h)
    viol = 0
    max_nonunit = 0
    nonunit_parts = set()
    rs = set()
    many_units = 0
    for A, cs in mem:
        normprod = math.prod(abs(int(F(a))) for a in A)
        r1 = sum(abs(c) for c in cs)
        bound = (r1 ** d) * fl2 ** (d - 1)
        res = abs(int(fmpz_poly(cs).resultant(F))) if hasattr(fmpz_poly, "resultant") else None
        if irreducible and (normprod > bound * (1 + 1e-12) or (res is not None and res != normprod)):
            viol += 1
        nu = tuple(sorted(a for a in A if abs(int(F(a))) >= 2))
        max_nonunit = max(max_nonunit, len(nu))
        nonunit_parts.add((tuple(cs), nu))
        rs.add(tuple(cs))
        if M - len(nu) >= M - 2:
            many_units += 1
    units = sorted(a for a in range(-B, B + 1) if abs(int(F(a))) == 1)
    return dict(F=[int(c) for c in F.coeffs()], d=d, M=M, B=B, h=h, irreducible=irreducible,
                n_members=len(mem), n_distinct_r=len(rs), n_informative=len(nonunit_parts),
                max_nonunit=max_nonunit, frac_le2_nonunit=(many_units / len(mem) if mem else None), unit_pieces=units, violations=viol)
