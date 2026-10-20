"""E1 (notes/almostall13_prereg.md): exact parameter count of design (a''') under N-averaging (Prop NB).

DIAGNOSTIC code: it builds fields from N exactly as the design does (src/thmA/numfield.make_f) and never factors N.
All linear algebra is exact over Z (python-flint fmpz_mat).

Usage: .venv/bin/python -m src.almostall13.rank experiments/raw/almostall13_e1.json
"""
from __future__ import annotations

import json
import random
import sys
from math import gcd
from typing import Dict, List, Sequence, Tuple

import flint

from src.thmA import arith
from src.thmA.numfield import hom_eval, make_f

N_BITS = 40
G_BITS = 10
C_BITS = 12


def affine_dim(vectors: Sequence[Sequence[int]]) -> int:
    """Dimension of the affine hull of the given integer vectors (exact)."""
    if len(vectors) <= 1:
        return 0
    v0 = vectors[0]
    rows = [[x - y for x, y in zip(v, v0)] for v in vectors[1:]]
    return int(flint.fmpz_mat(rows).rank())


def eval_matrix(points: Sequence[Tuple[int, int]], d: int) -> List[List[int]]:
    """E[i][j] = a_i^j b_i^(d-j), so that (F(P_i))_i = E f for f = (f_0, ..., f_d) (low -> high, hom_eval)."""
    return [[a ** j * b ** (d - j) for j in range(d + 1)] for (a, b) in points]


def left_kernel(points: Sequence[Tuple[int, int]], d: int) -> List[List[int]]:
    """Integer vectors lambda with lambda^T E = 0 (a basis of the rational left kernel, scaled to Z)."""
    E = eval_matrix(points, d)
    Et = flint.fmpz_mat([[E[i][j] for i in range(len(points))] for j in range(d + 1)])
    X, nul = Et.nullspace()
    out = []
    for k in range(nul):
        out.append([int(X[i, k]) for i in range(len(points))])
    return out


def draw_points(d: int, s: int) -> List[Tuple[int, int]]:
    rng = random.Random(f"almostall13|E1|pts|{d}|{s}")
    pts: List[Tuple[int, int]] = []
    while len(pts) < s:
        a = rng.randrange(-2 ** 10, 2 ** 10 + 1)
        b = rng.randrange(1, 2 ** 10 + 1)
        if gcd(a, b) != 1:
            continue
        if any(a * b2 - a2 * b == 0 for (a2, b2) in pts):
            continue
        pts.append((a, b))
    return pts


def _rand_N(rng: random.Random) -> int:
    return rng.randrange(2 ** (N_BITS - 1), 2 ** N_BITS)


def _rand_g(rng: random.Random) -> Tuple[int, int]:
    while True:
        g0, g1 = rng.randrange(2 ** (G_BITS - 1), 2 ** G_BITS), rng.randrange(2 ** (G_BITS - 1), 2 ** G_BITS)
        if gcd(g0, g1) == 1 and g1 >= 2:
            return g0, g1


def _rand_w(rng: random.Random) -> int:
    while True:
        w = rng.randrange(2 ** 5 + 1, 2 ** 6)
        if arith.is_prime(w):
            return w


def _rand_c(rng: random.Random, d: int) -> List[int]:
    return [rng.randrange(-2 ** C_BITS, 2 ** C_BITS + 1) for _ in range(d)]


def sample_regime(regime: str, d: int, K: int, seed: str) -> List[Dict]:
    """K samples (N, g0, g1, w, f) of design (a''') with the regime's coordinates held fixed."""
    rng = random.Random(seed)
    N0, (g00, g10), w0 = _rand_N(rng), _rand_g(rng), _rand_w(rng)
    out = []
    for _ in range(K):
        N = N0 if regime in ("A", "B") else _rand_N(rng)
        g0, g1 = (g00, g10) if regime in ("A", "B", "C") else _rand_g(rng)
        w = w0 if regime == "A" else _rand_w(rng)
        c = _rand_c(rng, d)
        f = make_f(N, w, d, g0, g1, c)
        out.append({"N": N, "g0": g0, "g1": g1, "w": w, "f": f})
    return out


def run_case(d: int, s: int) -> Dict:
    pts = draw_points(d, s)
    K = 3 * (d + 3) + 10
    ker = left_kernel(pts, d)
    res: Dict = {"d": d, "s": s, "K": K, "kernel_dim": len(ker)}
    viol = 0
    for reg in ("A", "B", "C", "D"):
        samp = sample_regime(reg, d, K, f"almostall13|E1|{reg}|{d}|{s}")
        vf = [[hom_eval(x["f"], a, b) for (a, b) in pts] for x in samp]
        res[f"aF_{reg}"] = affine_dim(vf)
        for v in vf:
            for lam in ker:
                if sum(l * y for l, y in zip(lam, v)) != 0:
                    viol += 1
        if reg == "D":
            vfg = [[hom_eval(x["f"], a, b) for (a, b) in pts] + [x["g1"] * a - x["g0"] * b for (a, b) in pts]
                   for x in samp]
            res["aFG_D"] = affine_dim(vfg)
    res["kernel_violations"] = viol
    res["R1"] = (res["aF_A"] == min(s, d) and all(res[f"aF_{r}"] == min(s, d + 1) for r in "BCD"))
    res["R2"] = res["aF_C"] == res["aF_B"]
    res["R3"] = res["aFG_D"] == min(s, d + 1) + min(s, 2)
    res["R4"] = (s <= d + 1) or viol == 0
    return res


def main():
    out = []
    for d in (3, 5):
        for s in (d - 1, d, d + 1, d + 2, 2 * d + 2):
            r = run_case(d, s)
            out.append(r)
            print(json.dumps(r), flush=True)
    summary = {"cases": out, "verdict": {k: all(r[k] for r in out) for k in ("R1", "R2", "R3", "R4")}}
    json.dump(summary, open(sys.argv[1], "w"), indent=1)
    print(json.dumps(summary["verdict"]))


if __name__ == "__main__":
    main()
