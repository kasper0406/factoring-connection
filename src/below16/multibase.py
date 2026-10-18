"""HARNESS: prime-coordinate reduction of multi-base searches (Prop. MB,
notes/below16_results.md Section 3).

Bases alpha_1..alpha_k are positive rationals (here: positive integers).  The
identity-provable relations among them are exactly the relations in Q^*, so in
prime-exponent coordinates (one coordinate per prime dividing some base) an
identity-provable collision of prod alpha_i^{u_i} with 1 means E^T u = 0, where
E is the k x r prime-exponent matrix.  A multi-base target with integer
exponent forms F_i(p, q) = B_i p + A_i q, offsets c_i and box
u in prod [0, W_i) covers (p, q) iff  E^T (F(p,q) - c) in E^T(box).
Each coordinate j gives a 1-D target L_j(F) = sum_i E_ij F_i with the finite
window pi_j(E^T box) (a sum of APs); coverage of the multi-base target implies
coverage of every coordinate target.
"""
from __future__ import annotations

import itertools
from typing import List, Sequence, Set, Tuple

import numpy as np
import sympy


def prime_exponents(bases: Sequence[int]) -> Tuple[List[int], np.ndarray]:
    primes = sorted({l for a in bases for l in sympy.factorint(a)})
    E = np.zeros((len(bases), len(primes)), dtype=np.int64)
    for i, a in enumerate(bases):
        for l, e in sympy.factorint(a).items():
            E[i, primes.index(l)] = e
    return primes, E


def image_of_box(E: np.ndarray, W: Sequence[int]) -> Set[Tuple[int, ...]]:
    out = set()
    for u in itertools.product(*[range(w) for w in W]):
        out.add(tuple(int(v) for v in np.array(u, dtype=np.int64) @ E))
    return out


def covered_multibase(p, q, forms, c, W, E) -> np.ndarray:
    T = image_of_box(E, W)
    F = np.stack([b * p + a * q - ci for (b, a), ci in zip(forms, c)], axis=1)
    V = F @ E
    return np.array([tuple(int(t) for t in row) in T for row in V], dtype=bool)


def covered_coordinate(p, q, forms, c, W, E, j) -> Tuple[np.ndarray, int]:
    T = image_of_box(E, W)
    Sj = {t[j] for t in T}
    F = np.stack([b * p + a * q - ci for (b, a), ci in zip(forms, c)], axis=1)
    Lj = F @ E[:, j]
    return np.isin(Lj, np.array(sorted(Sj), dtype=np.int64)), len(Sj)
