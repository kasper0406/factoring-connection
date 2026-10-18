"""E2 driver: structured GAP families at larger n, prime mode (notes/divcover_prereg.md §2).

(a) CRT/LLL designs: assign primes > L_1 to box points, LLL-reduce the (r+1)-dim lattice
    of (b, c) solving all congruences, verify the cover exactly.
(b) Polynomial-value (cyclotomic) GAPs {f(t) : deg f <= r, |coeffs| <= (K-1)/2}.

    .venv/bin/python -m src.divcover.run_structured
"""
from __future__ import annotations

import itertools
import json
import math
import os
import time

import numpy as np
import sympy
from fpylll import IntegerMatrix, LLL

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CFG = os.path.join(ROOT, "experiments/configs/divcover_structured.json")


def balanced_shape(r: int, L: int):
    side = max(2, round(L ** (1.0 / r)))
    shape = [side] * r
    # adjust the first side so that the product is close to L
    rest = math.prod(shape[1:])
    shape[0] = max(2, round(L / rest))
    return tuple(sorted(shape, reverse=True))


def box_points(shape):
    return np.array(list(itertools.product(*[range(l) for l in shape])), dtype=np.int64)


def gap_cover_stats(b: int, c: tuple, shape: tuple, primes):
    """Exact: which primes divide some NONZERO element of the GAP."""
    X = box_points(shape)
    # exact zero elements (python ints)
    zeros = 0
    for x in X:
        if b + sum(int(xi) * ci for xi, ci in zip(x, c)) == 0:
            zeros += 1
    covered = []
    for p in primes:
        cm = np.array([ci % p for ci in c], dtype=np.int64)
        v = (b % p + X @ cm) % p
        if int((v == 0).sum()) > zeros:
            covered.append(p)
    return covered, zeros


def gap_height(b: int, c: tuple, shape: tuple) -> int:
    hi = b + sum(max(0, (l - 1) * ci) for l, ci in zip(shape, c))
    lo = b + sum(min(0, (l - 1) * ci) for l, ci in zip(shape, c))
    return max(abs(hi), abs(lo))


def crt_design(n: int, r: int, L_target: int):
    shape = balanced_shape(r, L_target)
    L = math.prod(shape)
    X = box_points(shape)
    primes_all = list(sympy.primerange(2, n + 1))
    assigned = [p for p in primes_all if p > shape[0]]
    extra = []
    for attempt in range(6):
        targets = sorted(set(assigned) | set(extra), reverse=True)
        load = np.zeros(len(X))
        assign = {}
        for p in targets:
            j = int(np.argmin(load))
            assign[p] = j
            load[j] += math.log(p)
        # lattice: b + sum_k Xk ck == 0 mod M, Xk = CRT of x_{p,k}
        M = math.prod(targets)
        Xk = []
        for k in range(r):
            residues = [int(X[assign[p], k]) % p for p in targets]
            Xk.append(int(sympy.ntheory.modular.crt(targets, residues)[0]) % M)
        # scale columns: c_k weighted by (L_k - 1) to approximate the height norm
        w = [1] + [max(1, l - 1) for l in shape]
        rows = [[M * w[0]] + [0] * r]
        for k in range(r):
            row = [(-Xk[k]) % M * w[0]] + [0] * r
            row[1 + k] = w[1 + k]
            rows.append(row)
        A = IntegerMatrix.from_matrix(rows)
        LLL.reduction(A)
        best = None
        for i in range(A.nrows):
            v = [A[i, j] for j in range(r + 1)]
            bb = v[0] // w[0]
            cc = tuple(v[1 + k] // w[1 + k] for k in range(r))
            if all(ci == 0 for ci in cc):
                continue
            h = gap_height(bb, cc, shape)
            if best is None or h < best[0]:
                best = (h, bb, cc)
        h, bb, cc = best
        covered, zeros = gap_cover_stats(bb, cc, shape, primes_all)
        missing = sorted(set(primes_all) - set(covered))
        if not missing:
            break
        extra = sorted(set(extra) | set(missing))
    mass = sum(math.log(p) for p in targets)
    logH = math.log(h)
    return {
        "family": "crt_lll", "n": n, "rank": r, "shape": list(shape), "L": L,
        "logH": logH, "alpha_eff": math.log(logH) / math.log(n),
        "beta_eff": math.log(L) / (2 * math.log(n)),
        "mass_assigned": mass, "minkowski_logH_est": mass / (r + 1),
        "covered_all": not missing, "missing": missing[:20], "n_missing": len(missing),
        "zero_elements": zeros, "attempts": attempt + 1,
    }


def poly_family(n: int, t: int, r: int, K: int):
    shape = tuple([K] * (r + 1))
    c = tuple(t ** i for i in range(r + 1))
    b = -((K - 1) // 2) * sum(c)
    box = K ** (r + 1)
    primes = [p for p in sympy.primerange(2, n + 1) if p > box]
    rec = {"family": "poly_values", "n": n, "t": t, "rank": r + 1, "K": K, "box": box,
           "n_primes_above_box": len(primes)}
    if not primes:
        rec["vacuous"] = True
        return rec
    X = box_points(shape)
    vals = [b + sum(int(xi) * ci for xi, ci in zip(x, c)) for x in X]
    distinct_nonzero = len({v for v in vals if v != 0})
    covered, zeros = gap_cover_stats(b, c, shape, primes)
    expect = sum(1 - (1 - 1 / p) ** distinct_nonzero for p in primes)
    rec.update({
        "vacuous": False, "distinct_nonzero": distinct_nonzero, "zero_elements": zeros,
        "covered": len(covered), "expected_random": expect,
        "ratio_covered_over_expected": len(covered) / expect if expect else None,
        "logH": math.log(gap_height(b, c, shape)),
    })
    return rec


def main():
    cfg = json.load(open(CFG))
    out = os.path.join(ROOT, cfg["raw"])
    with open(out, "w") as fh:
        for n in cfg["crt"]["n_grid"]:
            for r in cfg["crt"]["ranks"]:
                for e in cfg["crt"]["L_exponents"]:
                    t0 = time.time()
                    rec = crt_design(n, r, round(n ** e))
                    rec["L_exponent_target"] = e
                    rec["wall_s"] = round(time.time() - t0, 2)
                    fh.write(json.dumps(rec) + "\n"); fh.flush()
                    print(json.dumps(rec), flush=True)
        for n in cfg["poly"]["n_grid"]:
            for t in cfg["poly"]["t_values"]:
                for r, K in cfg["poly"]["rK"]:
                    rec = poly_family(n, t, r, K)
                    fh.write(json.dumps(rec) + "\n"); fh.flush()
                    print(json.dumps(rec), flush=True)


if __name__ == "__main__":
    main()
