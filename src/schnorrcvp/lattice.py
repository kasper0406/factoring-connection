"""Schnorr/Adleman prime-number lattices as exact integer bases, and CVP calls through fpylll.

Agent `schnorrcvp`.  Operational inputs: N, n, public lattice parameters, a seed.  No factor of N
is used.  Logarithms enter only through correctly rounded integers (mpmath at 80 digits), so the
lattice is an exact integer lattice; all residual / relation checks are exact integer arithmetic.

Lattice (Schnorr 1991 / 2021 shape, CVP version):
    rows b_i = (d_i * e_i, c_i),  i = 1..n,   c_i = round(S * C * ln p_i)
    target t = (0, ..., 0, round(S * C * ln N))
  d_i = round(S * sqrt(weight_i)) with weight_i one of
    "lnp"  : ln p_i                       (Schnorr 1991, Adleman 1995: ||e||^2 ~ ln(uv) for squarefree)
    "perm" : f(i), f a random permutation of 1..n  (Schnorr 2021: diagonal f(i) [no square root])
    "sqrtperm": sqrt(f(i))                (Schnorr 2021 variant with sqrt; Yan et al.)
    "lnpperm": ln p_{f(i)}               (Aboumrad-Widdows-Kaushik variant: diag log p_{sigma(j)})
  C = N^c (given as the exponent c, or Schnorr's N' = N^{1/(n+1)}).
"""
from __future__ import annotations

import random
from typing import Dict, List, Optional, Sequence, Tuple

import mpmath as mp
from fpylll import BKZ, CVP, FPLLL, GSO, LLL, Enumeration, EnumerationError, IntegerMatrix
from fpylll.fplll.bkz_param import BKZParam

from .core import factor_over, first_primes, is_smooth, uvw

mp.mp.dps = 80


def _round(x) -> int:
    return int(mp.nint(x))


def build_lattice(N: int, n: int, *, c_exp=None, weight: str = "lnp", S: int = 1000, seed: int = 0,
                  schnorr_scale: bool = False):
    """Return (B, target, diag, P, meta) with B an IntegerMatrix (n rows, n+1 cols)."""
    P = first_primes(n)
    rng = random.Random(f"schnorrcvp|lattice|{N}|{n}|{weight}|{seed}")
    perm = list(range(1, n + 1))
    rng.shuffle(perm)
    if weight == "lnp":
        w = [mp.log(p) for p in P]
        d = [_round(S * mp.sqrt(x)) for x in w]
    elif weight == "perm":
        d = [S * perm[i] for i in range(n)]
    elif weight == "sqrtperm":
        d = [_round(S * mp.sqrt(perm[i])) for i in range(n)]
    elif weight == "lnpperm":
        d = [_round(S * mp.log(P[perm[i] - 1])) for i in range(n)]
    elif weight == "sqrtlnpperm":
        d = [_round(S * mp.sqrt(mp.log(P[perm[i] - 1]))) for i in range(n)]
    else:
        raise ValueError(weight)
    lnN = mp.log(N)
    if schnorr_scale:
        C = mp.mpf(N) ** (mp.mpf(1) / (n + 1))
    else:
        C = mp.mpf(N) ** mp.mpf(c_exp)
    col = [_round(S * C * mp.log(p)) for p in P]
    tgt_last = _round(S * C * lnN)
    B = IntegerMatrix(n, n + 1)
    for i in range(n):
        B[i, i] = d[i]
        B[i, n] = col[i]
    target = [0] * n + [tgt_last]
    meta = {"N": N, "n": n, "weight": weight, "S": S, "c_exp": None if schnorr_scale else str(c_exp),
            "schnorr_scale": schnorr_scale, "seed": seed}
    return B, target, d, P, meta


def coords_from_vector(vec: Sequence[int], d: Sequence[int]) -> List[int]:
    """Exponent vector e of a lattice vector (d_i e_i, *) -- exact division."""
    e = []
    for x, di in zip(vec[: len(d)], d):
        q, r = divmod(int(x), di)
        assert r == 0
        e.append(q)
    return e


def reduce_basis(B: IntegerMatrix, beta: int = 0):
    """LLL (beta = 0) or LLL followed by BKZ-beta.  In place."""
    LLL.reduction(B)
    if beta and beta > 2:
        BKZ.reduction(B, BKZParam(block_size=beta, max_loops=8,
                                  flags=BKZ.AUTO_ABORT | BKZ.MAX_LOOPS))
        LLL.reduction(B)
    return B


def exact_cvp(B: IntegerMatrix, target: Sequence[int]):
    """Exact closest vector (fplll enumeration, method='proved'); B must be LLL-reduced."""
    return tuple(int(x) for x in CVP.closest_vector(B, list(target), method="proved"))


def babai(B: IntegerMatrix, target: Sequence[int]):
    """Babai nearest plane on the given (reduced) basis; returns the lattice vector."""
    M = GSO.Mat(B)
    M.update_gso()
    coeffs = M.babai(list(target))
    v = B.multiply_left(coeffs)
    return tuple(int(x) for x in v)


def list_cvp(B: IntegerMatrix, target: Sequence[int], k: int, radius_factor: float = 1.01,
             precision: int = 240):
    """Up to k closest lattice vectors to target (fplll enumeration on an mpfr GSO), exactly re-sorted.

    Radius = radius_factor * (squared Babai distance), so the exact closest vector is always inside
    and is returned first (all candidates are re-checked and re-sorted with exact integers).  If fewer
    than k vectors lie in that radius, fewer are returned (the radius is enlarged up to 3 times).
    fplll takes the target in Gram-Schmidt coordinates (M.from_canonical)."""
    FPLLL.set_precision(precision)
    M = GSO.Mat(B, float_type="mpfr")
    M.update_gso()
    bv = babai(B, target)
    d0 = sum((a - b) ** 2 for a, b in zip(bv, target))
    rad = max(d0, 1) * radius_factor
    n = B.nrows
    tt = M.from_canonical(tuple(float(x) for x in target))
    sols = []
    for attempt in range(4):
        enum = Enumeration(M, nr_solutions=k)
        try:
            sols = enum.enumerate(0, n, float(rad), 0, target=tt)
        except EnumerationError:
            sols = []
        if len(sols) >= k or attempt == 3:
            break
        rad *= 2
    out = {bv}
    for dist, coeffs in sols:
        v = B.multiply_left([int(round(x)) for x in coeffs])
        out.add(tuple(int(x) for x in v))
    out = sorted(out, key=lambda v: sum((a - b) ** 2 for a, b in zip(v, target)))
    return out[:k]


def analyse_vector(vec, d, P, N, target) -> Dict:
    """Exact (u, v, w), smoothness of w, and the exact squared distance."""
    e = coords_from_vector(vec, d)
    u, v, w = uvw(e, P, N)
    dist2 = sum((a - b) ** 2 for a, b in zip(vec, target))
    rel = (w != 0) and is_smooth(w, P)
    return {"e": e, "u": u, "v": v, "w": w, "rel": rel, "dist2": dist2,
            "nnz": sum(1 for x in e if x), "trivial": all(x == 0 for x in e)}
