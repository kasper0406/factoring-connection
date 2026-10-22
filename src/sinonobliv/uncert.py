"""Uncertified outputs of a Howgrave-Graham / Coppersmith divisor lattice (sihard's 'not excluded' item 5).

Call (x, V, P, ell, s): find r in [0, V] with large d = gcd(x + r, P). Basis (f(r) = x + r):
  P^{ell-j} f^j  (j = 0..ell),   r^i f^ell  (i = 1..s);   coefficients of r^k scaled by V^k.
LLL (python-flint). Each reduced row b gives h_b(r) = sum_k (b_k / V^k) r^k in Z[r].
Certificate: if d^ell > sqrt(n) * ||b_1|| then h_{b_1}(r0) = 0 for every r0 in [0, V] with d | f(r0) (Howgrave-Graham).
Roots of h_{b_1..b_3} in [0, V] with gcd < T_cert := (sqrt(n)||b_1||)^{1/ell} are 'uncertified outputs'.
Operational: sees only x, V, P (public), y.
"""
from __future__ import annotations

import math
from typing import Dict, List

import flint
import numpy as np

from .core import is_smooth_int


def _poly_mul(a: List[int], b: List[int]) -> List[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        if u:
            for j, v in enumerate(b):
                out[i + j] += u * v
    return out


def hg_basis(x: int, V: int, P: int, ell: int, s: int) -> List[List[int]]:
    n = ell + s + 1
    f = [x, 1]
    fp = [[1]]
    for _ in range(ell):
        fp.append(_poly_mul(fp[-1], f))
    rows = []
    for j in range(ell + 1):
        c = [P ** (ell - j) * v for v in fp[j]]
        rows.append(c)
    for i in range(1, s + 1):
        rows.append([0] * i + fp[ell])
    B = []
    for c in rows:
        c = c + [0] * (n - len(c))
        B.append([c[k] * V ** k for k in range(n)])
    return B


def hg_call(x: int, V: int, P: int, ell: int, s: int, nvec: int = 3) -> Dict:
    B = hg_basis(x, V, P, ell, s)
    n = len(B)
    R = flint.fmpz_mat(B).lll()
    rows = [[int(R[i, k]) for k in range(n)] for i in range(min(nvec, n))]
    norm1 = math.isqrt(sum(v * v for v in rows[0])) + 1
    T_cert_log = (math.log(math.sqrt(n)) + math.log(norm1)) / ell
    roots = set()
    roots_first = set()
    for idx, b in enumerate(rows):
        coeffs = []
        for k in range(n):
            q, rem = divmod(b[k], V ** k)
            assert rem == 0
            coeffs.append(q)
        if not any(coeffs):
            continue
        p = flint.fmpz_poly(coeffs)
        for (rt, _mult) in p.roots():
            rt = int(rt)
            if 0 <= rt <= V:
                roots.add(rt)
                if idx == 0:
                    roots_first.add(rt)
    return {"n": n, "T_cert_log": T_cert_log, "roots": sorted(roots), "roots_first": sorted(roots_first),
            "h1_coeffs_first": rows[0]}


def classify_call(x: int, V: int, P: int, ell: int, s: int, y: int) -> Dict:
    res = hg_call(x, V, P, ell, s)
    Tlog = res["T_cert_log"]
    out = {"cert_roots": 0, "uncert_roots": 0, "uncert": [], "cert_violations": 0, "n_cert_points": 0}
    # sanity: every point with gcd >= T_cert must be a root of h_1
    r = np.arange(V + 1, dtype=object)
    for rr in range(V + 1):
        g = math.gcd(x + rr, P)
        if g > 1 and math.log(g) >= Tlog:
            out["n_cert_points"] += 1
            if rr not in res["roots_first"]:
                out["cert_violations"] += 1
    for rt in res["roots"]:
        z = x + rt
        g = math.gcd(z, P)
        rec = {"r": rt, "log_gcd": math.log(g), "smooth": is_smooth_int(z, y), "cof_log2": math.log2(z / g)}
        if g > 1 and math.log(g) >= Tlog:
            out["cert_roots"] += 1
        else:
            out["uncert_roots"] += 1
            out["uncert"].append(rec)
    out["T_cert_log"] = Tlog
    return out
