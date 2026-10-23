"""MODEL computations (no primes, no instances): Haar/Chebotarev constants on GL2(F_l), S3, S4.

n(t, d) = #{g in GL2(F_l): tr g = t, det g = d} = l (l + chi(t^2 - 4d)), chi the Legendre symbol, chi(0) = 0.
"""
from __future__ import annotations

from itertools import product
from math import log2
from typing import Dict, List, Tuple

import numpy as np


def legendre(a: int, l: int) -> int:
    a %= l
    if a == 0:
        return 0
    return 1 if pow(a, (l - 1) // 2, l) == 1 else -1


def class_counts(l: int) -> np.ndarray:
    """n[t, d] for t in F_l, d in F_l (d = 0 row is zero)."""
    n = np.zeros((l, l), dtype=np.int64)
    for t in range(l):
        for d in range(1, l):
            n[t, d] = l * (l + legendre(t * t - 4 * d, l))
    return n


def class_counts_bruteforce(l: int) -> np.ndarray:
    n = np.zeros((l, l), dtype=np.int64)
    for a, b, c, d in product(range(l), repeat=4):
        det = (a * d - b * c) % l
        if det:
            n[(a + d) % l, det] += 1
    return n


def _mi(joint: np.ndarray) -> float:
    """Mutual information (bits) between the two axes of a nonnegative 2-D array (normalised here)."""
    P = joint / joint.sum()
    px = P.sum(1, keepdims=True)
    py = P.sum(0, keepdims=True)
    m = P > 0
    return float((P[m] * np.log2(P[m] / (px @ py)[m])).sum())


def info_ordered(l: int) -> float:
    """I(det g_p ; (t_p, t_q) | det g_p det g_q) in bits, g_p, g_q iid Haar on GL2(F_l)."""
    n = class_counts(l).astype(float)
    tot = 0.0
    I = 0.0
    for c in range(1, l):
        # joint over d (= det g_p) and (t_p, t_q): n[t_p, d] n[t_q, c/d]
        J = np.zeros((l - 1, l * l))
        for i, d in enumerate(range(1, l)):
            dq = c * pow(d, -1, l) % l
            J[i] = np.outer(n[:, d], n[:, dq]).ravel()
        w = J.sum()
        I += w * _mi(J)
        tot += w
    return I / tot


def prod_likelihood(l: int) -> np.ndarray:
    """L[d, dq, u] = sum_{t t' = u} n[t, d] n[t', dq]  (unnormalised)."""
    n = class_counts(l).astype(float)
    L = np.zeros((l, l, l))
    for t in range(l):
        for t2 in range(l):
            u = t * t2 % l
            L[:, :, u] += np.outer(n[t, :], n[t2, :])
    return L


def info_product(l: int) -> float:
    """I(det g_p ; t_p t_q | det g_p det g_q) in bits (= information about s mod l, since symmetric)."""
    L = prod_likelihood(l)
    tot = 0.0
    I = 0.0
    for c in range(1, l):
        J = np.zeros((l - 1, l))
        for i, d in enumerate(range(1, l)):
            J[i] = L[d, c * pow(d, -1, l) % l, :]
        w = J.sum()
        I += w * _mi(J)
        tot += w
    return I / tot


# ------------------------------------------------------------------ fixed-l count distributions (model)

def gl2_elements(l: int):
    for a, b, c, d in product(range(l), repeat=4):
        if (a * d - b * c) % l:
            yield (a, b, c, d)


def dist_M_mod_l(l: int) -> Dict[int, float]:
    n = class_counts(l).astype(float)
    tot = n.sum()
    one: Dict[int, float] = {}
    for t in range(l):
        for d in range(1, l):
            v = (1 - t + d) % l
            one[v] = one.get(v, 0.0) + n[t, d] / tot
    return _prod_dist(one, l)


def dist_aN_mod_l(l: int) -> Dict[int, float]:
    n = class_counts(l).astype(float)
    tot = n.sum()
    one: Dict[int, float] = {}
    for t in range(l):
        one[t] = n[t].sum() / tot
    return _prod_dist(one, l)


def _prod_dist(one: Dict[int, float], mod: int = 0) -> Dict[int, float]:
    out: Dict[int, float] = {}
    for a, pa in one.items():
        for b, pb in one.items():
            v = (a * b) % mod if mod else a * b
            out[v] = out.get(v, 0.0) + pa * pb
    return out


def eigenlines_F3(g) -> int:
    """Number of the 4 lines of F_3^2 stable under g (= roots of psi_3 mod r when g = Frob_r on E[3])."""
    a, b, c, d = g
    lines = [(1, 0), (0, 1), (1, 1), (1, 2)]
    cnt = 0
    for (x, y) in lines:
        gx, gy = (a * x + b * y) % 3, (c * x + d * y) % 3
        if (gx * y - gy * x) % 3 == 0:
            cnt += 1
    return cnt


def psi3_table() -> Dict[Tuple[int, int], int]:
    """Counts over GL2(F_3) of (det g, #stable lines)."""
    T: Dict[Tuple[int, int], int] = {}
    for g in gl2_elements(3):
        a, b, c, d = g
        key = ((a * d - b * c) % 3, eigenlines_F3(g))
        T[key] = T.get(key, 0) + 1
    return T


def dist_r_psi3() -> Dict[int, float]:
    T = psi3_table()
    tot = sum(T.values())
    one: Dict[int, float] = {}
    for (det, r), c in T.items():
        one[r] = one.get(r, 0.0) + c / tot
    return _prod_dist(one)


def dist_r_S3() -> Dict[int, float]:
    one = {3: 1 / 6, 1: 1 / 2, 0: 1 / 3}
    return _prod_dist(one)


def pt_floor(dist: Dict[int, float]) -> float:
    """1 - max_b Pr[F = b]; Theorem PT: a pure tree computing F exactly splits with probability >= half of this."""
    return 1.0 - max(dist.values())


def adv_psi3_beta() -> float:
    """Bayes advantage for (p/3) on N = 1 mod 3 from r_N(psi_3) (model)."""
    T = psi3_table()
    tot = sum(T.values())
    joint: Dict[int, List[float]] = {}
    mass = 0.0
    for (d1, r1), c1 in T.items():
        for (d2, r2), c2 in T.items():
            if d1 != d2:   # N = d1 d2 = 1 mod 3 iff d1 = d2 (F_3^* = {1, 2})
                continue
            w = c1 * c2 / tot ** 2
            joint.setdefault(r1 * r2, [0.0, 0.0])[0 if d1 == 1 else 1] += w
            mass += w
    return sum(max(v) for v in joint.values()) / mass - 0.5


def adv_aN_beta(l: int) -> float:
    """Bayes advantage for (p/l) on (N/l) = 1 from (N mod l, a_N mod l) (model)."""
    n = class_counts(l).astype(float)
    tot = n.sum()
    joint: Dict[Tuple[int, int], List[float]] = {}
    mass = 0.0
    for d1 in range(1, l):
        for d2 in range(1, l):
            s1, s2 = legendre(d1, l), legendre(d2, l)
            if s1 != s2:
                continue
            for t1 in range(l):
                for t2 in range(l):
                    w = n[t1, d1] * n[t2, d2] / tot ** 2
                    key = (d1 * d2 % l, t1 * t2 % l)
                    joint.setdefault(key, [0.0, 0.0])[0 if s1 == 1 else 1] += w
                    mass += w
    return sum(max(v) for v in joint.values()) / mass - 0.5


# ------------------------------------------------------------------ Haar simulation of the decoders

def sample_traces(rng: np.random.Generator, l: int, d: int, K: int) -> np.ndarray:
    n = class_counts(l)[:, d].astype(float)
    return rng.choice(l, size=K, p=n / n.sum())


# ------------------------------------------------------------------ identifiability (post hoc, PH1/PH2)

def _law_ordered(l: int, d: int, c: int) -> np.ndarray:
    n = class_counts(l).astype(float)
    J = np.outer(n[:, d], n[:, c * pow(d, -1, l) % l])
    return J / J.sum()


def law(l: int, kind: str, d: int, c: int) -> np.ndarray:
    """Law of the per-curve observation given p = d, q = c/d (mod l), N = c.  kind in
    {'ordered' (t_p, t_q), 'unordered' {t_p, t_q}, 'product' t_p t_q, 'M' det(1-g_p)det(1-g_q), 'single' t_p}."""
    J = _law_ordered(l, d, c)
    if kind == "ordered":
        return J.ravel()
    if kind == "single":
        return J.sum(1)
    if kind == "unordered":
        U = J + J.T - np.diag(np.diag(J))   # {t, u}: J[t,u] + J[u,t] for t < u, J[t,t] on the diagonal (cx2 fix)
        return U[np.triu_indices(l)]
    if kind == "product":
        out = np.zeros(l)
        for t in range(l):
            for t2 in range(l):
                out[t * t2 % l] += J[t, t2]
        return out
    if kind == "M":
        n = class_counts(l).astype(float)
        dq = c * pow(d, -1, l) % l
        a = np.zeros(l)
        b = np.zeros(l)
        for t in range(l):
            a[(1 - t + d) % l] += n[t, d]
            b[(1 - t + dq) % l] += n[t, dq]
        a /= a.sum()
        b /= b.sum()
        out = np.zeros(l)
        for x in range(l):
            for y in range(l):
                out[x * y % l] += a[x] * b[y]
        return out
    raise ValueError(kind)


def identifiable_classes(l: int, kind: str, c: int, tol: float = 1e-12) -> List[List[int]]:
    """Partition of d in F_l^* into classes with identical per-curve laws (what infinitely many independent curves
    can identify about p mod l, given N mod l = c)."""
    classes: List[List[int]] = []
    reps: List[np.ndarray] = []
    for d in range(1, l):
        v = law(l, kind, d, c)
        for i, r in enumerate(reps):
            if np.abs(v - r).max() < tol:
                classes[i].append(d)
                break
        else:
            reps.append(v)
            classes.append([d])
    return classes


def info_kind(l: int, kind: str) -> float:
    """I(p mod l ; one curve's observation | N mod l) in bits, Haar model (post hoc PH2)."""
    I = 0.0
    for c in range(1, l):
        J = np.array([law(l, kind, d, c) for d in range(1, l)])
        I += _mi(J)
    return I / (l - 1)
