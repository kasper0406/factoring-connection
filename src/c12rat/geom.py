"""Exact lattice geometry for the child boxes (notes/c12rat_results.md §3: Lemmas LS, MULT, TF).  DIAGNOSTIC.

label (q, r): Lambda = {(a, b) : a = r b mod q}, Gauss-reduced basis e1, e2 (|e1| <= |e2|).
Box(q, r): primitive P = i e1 + j e2 with gcd(i, j) = 1, max(|i|, |j|) <= S, P primitive in Z^2, modulo +-1
(the conservative geometry of uprime §4.2 / pc3 §1).  We reuse pc3's toy code by import (read-only).
"""
from __future__ import annotations

from math import gcd, isqrt, pi, sqrt
from typing import Dict, List, Sequence, Tuple

import numpy as np

from src.pc3 import toy as pc3toy

label_basis = pc3toy.label_basis
child_box = pc3toy.child_box


def shape(q: int, r: int) -> Tuple[float, float]:
    """(sigma, t2): sigma = lambda1^2 / q, t2 = lambda1 / lambda2 (<= 1)."""
    e1, e2 = label_basis(q, r)
    l1 = sqrt(e1[0] ** 2 + e1[1] ** 2)
    l2 = sqrt(e2[0] ** 2 + e2[1] ** 2)
    if l1 > l2:
        l1, l2 = l2, l1
    return l1 * l1 / q, l1 / l2


def shape_cdf_bound(sigma: float, q: int) -> float:
    """Lemma LS (proved): Pr_r[lambda1^2 <= sigma q] <= pi (sqrt(sigma q) + 1)^2 / (2 q)."""
    return pi * (sqrt(sigma * q) + 1.0) ** 2 / (2.0 * q)


def box_arrays(q: int, r: int, S: int):
    """Box points as arrays: coords (i, j) and vectors (a, b)."""
    e1, e2 = label_basis(q, r)
    box = child_box(e1, e2, S)
    ij = np.array([(b[0], b[1]) for b in box], dtype=np.int64)
    ab = np.array([(b[2], b[3]) for b in box], dtype=np.int64)
    return ij, ab


def pair_stats(ij: np.ndarray, ab: np.ndarray):
    """For all unordered pairs (k < l): Delta = |i j' - j i'| (coordinate determinant), |det| = q Delta (checked by
    caller), kappa = |P'||P''| / |det(P', P'')|."""
    n = len(ij)
    iu, ju = np.triu_indices(n, k=1)
    Delta = np.abs(ij[iu, 0] * ij[ju, 1] - ij[iu, 1] * ij[ju, 0])
    det = np.abs(ab[iu, 0] * ab[ju, 1] - ab[iu, 1] * ab[ju, 0])
    norm = np.sqrt((ab[:, 0].astype(float) ** 2 + ab[:, 1].astype(float) ** 2))
    kappa = norm[iu] * norm[ju] / det
    return iu, ju, Delta, det, kappa


def thin_fraction(kappa: np.ndarray, Ks: Sequence[float]) -> Dict[float, float]:
    return {float(K): float((kappa > K).mean()) for K in Ks}


def thin_bound(S: int, n: int, t2: float, K: float) -> float:
    """Lemma TF (proved): fraction of unordered pairs with kappa > K is
    <= min(1, ceil(2 sqrt3 S^2/(t2 K)) * (10 S^2 + 2 S) / (n (n - 1)))   (t2 = lambda1/lambda2)."""
    X = int(np.ceil(2.0 * sqrt(3.0) * S * S / (t2 * K)))
    return min(1.0, X * (10.0 * S * S + 2.0 * S) / (n * (n - 1.0)))


def mult_max(Delta: np.ndarray) -> int:
    """max over Delta >= 1 of #unordered pairs with that coordinate determinant (Lemma MULT: <= 5 S^2 + S)."""
    d = Delta[Delta > 0]
    if len(d) == 0:
        return 0
    return int(np.bincount(d).max())
