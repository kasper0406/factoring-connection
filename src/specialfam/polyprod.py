"""Resultant packaging: block products of polynomial values by one monic remainder tree.

OPERATIONAL. Every function here receives only N and public data (polynomials and
point sets built from public parameters). Nothing uses p, q, phi(N) or a factoring
routine. The only flint operations on (Z/N)[Y] are multiplication and ``%`` by a
MONIC polynomial (src/common/README.md flint audit). No element of Z/N is ever
inverted here, except leading coefficients handled by ``monic_from`` (gcd-checked).

Identity used (Theorem G of notes/specialfam_results.md):
  for a finite S and monic chi of degree d over any commutative ring,
      prod_{s in S} chi(s) = Res(F_S, chi),   F_S(Y) = prod_{s in S} (Y - s),
      Res(F_S, chi) = (-1)^{|S| d} Res(chi, F_S) = (-1)^{|S| d} Nm_{A/(Z/N)}(F_S mod chi),
  where A = (Z/N)[Y]/(chi) and Nm is the determinant of multiplication (division-free).
So the K block products prod_{s in S} chi_k(s) cost one product tree for F_S, one
product tree of the chi_k, one remainder tree, and K small determinants.
"""
from __future__ import annotations

from math import gcd
from typing import List, Optional, Sequence, Tuple

import flint

from src.common.metrics import Metrics, NullMetrics
from src.common.ring import NonUnit


def poly_ctx(N: int):
    return flint.fmpz_mod_poly_ctx(N)


def _coeffs(P, d: int) -> List[int]:
    cs = [int(c) for c in P.coeffs()]
    return cs + [0] * (d - len(cs))


def product_tree(ctx, leaves: Sequence, mt: Metrics) -> List[List]:
    """Levels of a product tree, levels[0] = leaves, levels[-1] = [root]."""
    levels = [list(leaves)]
    while len(levels[-1]) > 1:
        cur = levels[-1]
        nxt = []
        for i in range(0, len(cur) - 1, 2):
            nxt.append(cur[i] * cur[i + 1])
            mt.inc("poly_mul")
            mt.inc("poly_mul_deg", cur[i].degree() + cur[i + 1].degree())
        if len(cur) % 2:
            nxt.append(cur[-1])
        levels.append(nxt)
    return levels


def remainder_tree(F, levels: List[List], mt: Metrics) -> List:
    """F mod every leaf of a product tree of MONIC polynomials."""
    rems = [F % levels[-1][0]]
    mt.inc("poly_mod")
    mt.inc("poly_mod_deg", max(F.degree(), 0))
    for lev in range(len(levels) - 2, -1, -1):
        cur = levels[lev]
        nxt = []
        for i, node in enumerate(cur):
            parent = rems[i // 2]
            if parent.degree() < node.degree():
                nxt.append(parent)  # already reduced (also covers odd carry)
            else:
                nxt.append(parent % node)
                mt.inc("poly_mod")
                mt.inc("poly_mod_deg", parent.degree())
        rems = nxt
    return rems


def det_mod(M: List[List[int]], N: int) -> int:
    """Division-free determinant (Laplace) mod N; used for d <= 6 only."""
    n = len(M)
    if n == 1:
        return M[0][0] % N
    if n == 2:
        return (M[0][0] * M[1][1] - M[0][1] * M[1][0]) % N
    tot = 0
    for j in range(n):
        if M[0][j] % N == 0:
            continue
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        term = M[0][j] * det_mod(minor, N)
        tot += -term if j % 2 else term
    return tot % N


def norm_mod(chi: List[int], r: List[int], N: int) -> int:
    """Nm of r(Y) in (Z/N)[Y]/(chi), chi monic (low->high, len d+1), deg r < d.

    Equals Res(chi, r) = prod over roots theta of chi of r(theta). Division-free.
    """
    d = len(chi) - 1
    if d == 1:
        return r[0] % N
    if d == 2:
        c0, c1 = chi[0], chi[1]
        r0, r1 = r[0], r[1]
        return (r0 * r0 - c1 * r0 * r1 + c0 * r1 * r1) % N
    # multiplication matrix: column j = Y^j * r mod chi
    cols = []
    cur = list(r) + [0] * (d - len(r))
    for _ in range(d):
        cols.append(cur)
        # cur *= Y mod chi
        top = cur[-1]
        nxt = [0] + cur[:-1]
        nxt = [(nxt[i] - top * chi[i]) % N for i in range(d)]
        cur = nxt
    M = [[cols[j][i] for j in range(d)] for i in range(d)]
    return det_mod(M, N)


def block_products(N: int, S: Sequence[int], chis: Sequence[Sequence[int]],
                   metrics: Optional[Metrics] = None) -> List[int]:
    """Return [prod_{s in S} chi_k(s) mod N for each k]; chi_k monic int lists low->high."""
    mt = metrics if metrics is not None else NullMetrics()
    ctx = poly_ctx(N)
    # F_S(Y) = prod (Y - s)
    lin = [ctx([(-s) % N, 1]) for s in S]
    FS = product_tree(ctx, lin, mt)[-1][0]
    mt.hi("max_degree", FS.degree())
    leaves = []
    for c in chis:
        if c[-1] % N != 1:
            raise ValueError("chi must be monic")
        leaves.append(ctx([x % N for x in c]))
    levels = product_tree(ctx, leaves, mt)
    mt.hi("max_degree", levels[-1][0].degree())
    rems = remainder_tree(FS, levels, mt)
    out = []
    nS = len(S)
    for c, rm in zip(chis, rems):
        d = len(c) - 1
        v = norm_mod([x % N for x in c], _coeffs(rm, d), N)
        mt.inc("small_norm")
        if (nS * d) % 2:
            v = (-v) % N
        out.append(v)
    mt.inc("stored_coeffs", sum(len(c) for c in chis) + nS)
    return out


def monic_from(N: int, coeffs: Sequence[int]) -> Tuple[List[int], int]:
    """Scale an integer polynomial to be monic mod N. Returns (monic, lc).

    Raises NonUnit (rule 3: the caller outputs a proper divisor) if gcd(lc, N) > 1.
    """
    lc = int(coeffs[-1]) % N
    if lc == 1:
        return [int(c) % N for c in coeffs], 1
    g = gcd(lc, N)
    if g != 1:
        raise NonUnit(g, N, lc)
    inv = pow(lc, -1, N)
    return [(int(c) * inv) % N for c in coeffs], lc
