"""sieve2d Theorem 1 (two-variable fourth moment): the explicit error bound, an
exhaustive checker for small boxes, and the theorem-scale evaluation.

HARNESS / DIAGNOSTIC code (it enumerates integer points; it is not an algorithm).
Notation as in notes/sieve2d_results.md §2:
  box   B = (0, X]^2 minus the lines L = U_v R v
  z_l(x) = #{v in V : l | B a - A b},  W(x) = #{l in band : z_l(x) > 0}
  s = sum 1/l, lam = |V| s,  E(Z)_j = j! e_j(|V|/l)  (independent Bernoulli model)

Theorem 1:  #{x in B \\ L : W = 0} <= [X^2 (3 lam^2 + lam) + 2 H X lam^4
                                   + sum_j kappa_j Err_j] / lam^4
with kappa = (15 lam^3, 25 lam^2, 10 lam, 1) for j = 1..4 and Err_j = err_terms(...)[j].
"""
from __future__ import annotations

from math import comb, gcd, sqrt
from typing import Dict, List, Sequence, Tuple

import numpy as np


def err_terms(H: int, X: int, Y: float, band: Sequence[int], V: Sequence[Tuple[int, int]],
              j: int) -> Dict[str, float]:
    """Explicit upper bounds (float) for the pieces of Err_j, j = 1..4 (notes §2.4)."""
    K = len(band)
    nV = len(V)
    s = sum(1.0 / l for l in band)
    lam = nV * s
    up = -(-j // 2)          # ceil(j/2)
    dn = j // 2
    t = {}
    t["A_plus1"] = float(K * nV) ** j
    t["B_boundary"] = 28.3 * X * (2 * Y) ** (j / 2) * (2 * nV) ** j
    t["C_otherlines"] = 2 * H * X * sum(comb(j, r) * lam ** r * K ** (j - r)
                                        for r in range(up, j + 1))
    t["D_pure4"] = X * X * nV * s ** 4 if j == 4 else 0.0
    rr = range(max(0, j - 3), dn + 1)
    t["E_lines1"] = 6 * H * X * sum(comb(j, r) * lam ** r * K ** (j - r) for r in rr)
    t["E_lines2"] = 2 * H * X * nV * sum(comb(j, r) * s ** (j - r) * (K * nV) ** r for r in rr)
    t["F_D"] = (3 * nV) ** j * X * X * s / Y
    t["total"] = sum(t.values())
    return t


def theorem1_bound(H, X, Y, band, V) -> Dict[str, float]:
    s = sum(1.0 / l for l in band)
    lam = len(V) * s
    assert lam >= 1, "Theorem 1 is stated for lam >= 1"
    kap = {1: 15 * lam ** 3, 2: 25 * lam ** 2, 3: 10 * lam, 4: 1.0}
    errs = {j: err_terms(H, X, Y, band, V, j)["total"] for j in (1, 2, 3, 4)}
    main = X * X * (3 * lam ** 2 + lam)
    err = 2 * H * X * lam ** 4 + sum(kap[j] * errs[j] for j in errs)
    return {"lam": lam, "main": main / lam ** 4, "err": err / lam ** 4,
            "bound": (main + err) / lam ** 4, "errs": errs}


def hypotheses_ok(H, X, Y, band) -> bool:
    return Y >= 16 * H * H and H * X < Y ** 4 and all(Y < l <= 2 * Y for l in band)


# ---------------------------------------------------------------------------
# exhaustive check on small boxes
# ---------------------------------------------------------------------------
def model_factorial_moments(band, nV, jmax=4) -> List[float]:
    """E(Z)_j = j! e_j(p_l), p_l = nV / l (float)."""
    e = [1.0] + [0.0] * jmax
    for l in band:
        p = nV / l
        for j in range(jmax, 0, -1):
            e[j] += e[j - 1] * p
    out = []
    f = 1
    for j in range(jmax + 1):
        if j:
            f *= j
        out.append(e[j] * f)
    return out


def exhaustive_stats(H, X, band, V):
    """Exact sums over x in (0,X]^2 \\ L of (W)_j, j = 0..4, uncovered count, and the
    sum of (W - lam)^4."""
    a = np.arange(1, X + 1, dtype=np.int64)
    A2, B2 = np.meshgrid(a, a, indexing="ij")
    xa, xb = A2.ravel(), B2.ravel()
    online = np.zeros(xa.shape, dtype=bool)
    for (A, B) in V:
        online |= (B * xa == A * xb)
    xa, xb = xa[~online], xb[~online]
    W = np.zeros(xa.shape, dtype=np.int64)
    for l in band:
        z = np.zeros(xa.shape, dtype=bool)
        for (A, B) in V:
            z |= ((B * xa - A * xb) % l) == 0
        W += z
    sums = []
    for j in range(5):
        prod = np.ones(W.shape, dtype=np.int64)
        for i in range(j):
            prod = prod * (W - i)
        sums.append(int(prod.sum()))
    s = sum(1.0 / l for l in band)
    lam = len(V) * s
    m4 = float((((W.astype(float)) - lam) ** 4).sum())
    return {"npts": int(W.size), "uncovered": int((W == 0).sum()), "fact_sums": sums,
            "sum_dev4": m4, "lam": lam}


def check_box(H, X, Y, band, V) -> Dict:
    st = exhaustive_stats(H, X, band, V)
    mod = model_factorial_moments(band, len(V))
    actual = [abs(st["fact_sums"][j] - X * X * mod[j]) for j in range(5)]
    bnd = [2 * H * X] + [err_terms(H, X, Y, band, V, j)["total"] for j in (1, 2, 3, 4)]
    out = {"H": H, "X": X, "Y": Y, "K": len(band), "lam": st["lam"],
           "actual_err": actual, "bound_err": bnd,
           "uncovered": st["uncovered"], "npts": st["npts"]}
    if st["lam"] >= 1:
        out["thm1_bound"] = theorem1_bound(H, X, Y, band, V)["bound"]
        out["dev4_bound"] = out["thm1_bound"] * st["lam"] ** 4
        out["sum_dev4"] = st["sum_dev4"]
    return out


# ---------------------------------------------------------------------------
# theorem-scale evaluation (PNT approximations; illustration only)
# ---------------------------------------------------------------------------
def scale_terms(n: int, H: int, k: int = 0) -> Dict[str, float]:
    """Normalised pieces of Theorem 1 at the D' scale, using |P| ~ Y/(2 ln Y)
    (primes 3 mod 4 in (Y, 2Y]) and s ~ ln2/(2 ln Y).  Each entry is the term divided
    by X^2 (so the uncovered integer fraction is <= the sum)."""
    from math import log
    P = 2.0 ** (n / 2)
    Y = 2 ** k * P ** (1 / 3) / H ** (1 / 3)
    X = sqrt(2) * P
    K = Y / (2 * log(Y))
    s = log(2) / (2 * log(Y))
    nV = 0.6079 * H * H
    lam = nV * s
    kap = {1: 15 * lam ** 3, 2: 25 * lam ** 2, 3: 10 * lam, 4: 1.0}
    out = {"lam": lam, "main": (3 * lam ** 2 + lam) / lam ** 4, "j0": 2 * H / X}
    for j in (1, 2, 3, 4):
        # synthetic band with the right K and s for err_terms
        t = {}
        up, dn = -(-j // 2), j // 2
        t["A"] = (K * nV) ** j
        t["B"] = 28.3 * X * (2 * Y) ** (j / 2) * (2 * nV) ** j
        t["C"] = 2 * H * X * sum(comb(j, r) * lam ** r * K ** (j - r) for r in range(up, j + 1))
        t["D"] = X * X * nV * s ** 4 if j == 4 else 0.0
        rr = range(max(0, j - 3), dn + 1)
        t["E1"] = 6 * H * X * sum(comb(j, r) * lam ** r * K ** (j - r) for r in rr)
        t["E2"] = 2 * H * X * nV * sum(comb(j, r) * s ** (j - r) * (K * nV) ** r for r in rr)
        t["F"] = (3 * nV) ** j * X * X * s / Y
        for key, v in t.items():
            out[f"j{j}_{key}"] = kap[j] * v / (X * X * lam ** 4)
    out["err_total"] = sum(v for kk, v in out.items() if kk.startswith("j"))
    out["hyp_Y16H2"] = Y >= 16 * H * H
    out["hyp_HX_lt_Y4"] = H * X < Y ** 4
    return out
