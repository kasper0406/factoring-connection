"""Reproductions of the r3 audit's Track B/D items (agent `repair`).  DIAGNOSTIC / HARNESS.

Imports `src.trackBD.linalg` read-only where the claim is about that code's behaviour.
Run:  .venv/bin/python -m src.repair.audit_trackBD --out experiments/raw/repair_audit_trackBD.jsonl
"""
from __future__ import annotations

import argparse
import json
import random
from fractions import Fraction
from math import comb, gcd
from typing import Dict, List

import flint

from src.common import instances as I
from src.repair.c2r_search import is_prime
from src.trackBD.linalg import bm_mod_N, elim_profile


def rref_pivots(M: List[List[int]], p: int) -> List[int]:
    A = [[x % p for x in row] for row in M]
    piv, row = [], 0
    for j in range(len(A[0])):
        k = next((i for i in range(row, len(A)) if A[i][j]), None)
        if k is None:
            continue
        A[row], A[k] = A[k], A[row]
        inv = pow(A[row][j], -1, p)
        A[row] = [x * inv % p for x in A[row]]
        for i in range(len(A)):
            if i != row and A[i][j]:
                t = A[i][j]
                A[i] = [(a - t * b) % p for a, b in zip(A[i], A[row])]
        piv.append(j)
        row += 1
    return piv


def rank_vs_pivot() -> Dict:
    M = [[3, 5]]
    res = elim_profile(M, 15)
    # Krylov claim for multiplication-by-S in the power basis: pivots = {0..rank-1}
    rng = random.Random(0)
    bad = n = 0
    for _ in range(3000):
        p = rng.choice([3, 5, 7, 11])
        r = rng.randint(2, 6)
        f = [rng.randrange(p) for _ in range(r)] + [1]
        S = [rng.randrange(p) for _ in range(r)]  # small p: rank deficiency is frequent
        cols = []
        cur = S[:]
        for _j in range(r):
            cols.append(cur[:])
            # multiply cur by X mod f
            top = cur[-1]
            cur = [0] + cur[:-1]
            cur = [(cur[i] - top * f[i]) % p for i in range(r)]
        mat = [[cols[j][i] for j in range(r)] for i in range(r)]
        pv = rref_pivots(mat, p)
        n += 1
        bad += pv != list(range(len(pv)))
    return {"item": "BD-rank-vs-pivot", "matrix": M, "N": 15, "elim_result": res,
            "pivots_mod_3": rref_pivots(M, 3), "pivots_mod_5": rref_pivots(M, 5),
            "krylov_mult_matrices_checked": n, "krylov_violations": bad}


def roots_Fe(p: int, e: int) -> List[int]:
    return [x for x in range(p) if (pow(1 + x, e, p) - 1 - pow(x, e, p)) % p == 0]


def is_frob_d1(x: int, p: int, e: int) -> bool:
    # FROB (d = 1): e = 1 mod L_x, i.e. x^(e-1) = 1 = (1+x)^(e-1) (x not in {0, -1})
    return pow(x, e - 1, p) == 1 and pow(1 + x, e - 1, p) == 1


def sporadic_six(pmax: int = 400) -> Dict:
    p, q = 7, 11
    e = q % (p - 1)
    rts = roots_Fe(p, e)
    lit = [x for x in rts if x not in (0, p - 1) and not is_frob_d1(x, p, e)]
    omega = [x for x in range(p) if (x * x + x + 1) % p == 0]
    lit_bad = rep_bad = n = 0
    for pp in range(5, pmax):
        if not is_prime(pp):
            continue
        om = {x for x in range(pp) if (x * x + x + 1) % pp == 0}
        for ee in range(1, pp - 1):
            if gcd(ee, pp - 1) != 1:
                continue
            r = roots_Fe(pp, ee)
            literal = [x for x in r if x not in (0, pp - 1) and not is_frob_d1(x, pp, ee)]
            repaired = [x for x in literal if x not in om]
            n += 1
            lit_bad += len(literal) % 6 != 0
            rep_bad += len(repaired) % 6 != 0
    return {"item": "BD-six", "p": p, "q": q, "e": e, "roots": rts, "cube_roots_of_unity": omega,
            "literal_SPORADIC": lit, "literal_count": len(lit),
            "census_(p,e)": n, "literal_def_not_0_mod_6": lit_bad, "repaired_def_not_0_mod_6": rep_bad}


def _side(l: int, m: int, r: int):
    """Own rank distribution of multiplication-by-sketch over admissible c in F_l."""
    X = flint.nmod_poly([0, 1], l)
    one = flint.nmod_poly([1], l)
    dist = {}
    adm = 0
    for c in range(l):
        if c == 0 or (pow(-1, r) - c) % l == 0:
            continue
        if r % 3 == 0:
            if (1 - c) % l == 0:
                continue
        elif (c * c + c + 1) % l == 0:
            continue
        adm += 1
        f = flint.nmod_poly([(-c) % l] + [0] * (r - 1) + [1], l)
        F = ((X + one).pow_mod(m, f) - one - X.pow_mod(m, f)) % f
        dg = r if F.is_zero() else f.gcd(F).degree()
        dist[r - dg] = dist.get(r - dg, 0) + 1
    return adm, dist


def pgap(N: int, p: int, q: int, r: int) -> Fraction:
    ap, dp = _side(p, q, r)
    aq, dq = _side(q, p, r)
    same = sum(dp[k] * dq.get(k, 0) for k in dp)
    return 1 - Fraction(same, ap * aq)


def n209() -> Dict:
    N, p, q = 209, 11, 19
    out = {"item": "BD-B6B7", "N": N}
    for l, m in ((p, q), (q, p)):
        X = flint.nmod_poly([0, 1], l)
        one = flint.nmod_poly([1], l)
        f = flint.nmod_poly([1, 0, 0, 0, 0, 0, 1], l)
        F = ((X + one).pow_mod(m, f) - one - X.pow_mod(m, f)) % f
        out[f"gcd(X^6+1, F_{m}) mod {l}"] = str(f.gcd(F))
    out["P_gap_r2"] = str(pgap(N, p, q, 2))
    out["P_gap_r6"] = str(pgap(N, p, q, 6))
    try:
        from src.trackBD.exp_rankdist import side as their_side
        for r in (2, 6):
            sp, sq = their_side(p, q, r), their_side(q, p, r)
            num = sum(v * sq["dist"].get(k, 0) for k, v in sp["dist"].items())
            out[f"their_P_gap_r{r}"] = str(1 - Fraction(num, sp["adm"] * sq["adm"]))
    except Exception as exc:  # pragma: no cover
        out["their_side_error"] = repr(exc)
    return out


def n247() -> Dict:
    N = 247
    out = {"item": "BD-M3b", "N": N, "f": "X^3+5"}
    for l in (13, 19):
        rts = [x for x in range(l) if (x ** 3 + 5) % l == 0]
        f = flint.nmod_poly([5, 0, 0, 1], l)
        X = flint.nmod_poly([0, 1], l)
        cols = [X.pow_mod(i * N, f) for i in range(3)]
        Q = [[int(cols[j][i]) for j in range(3)] for i in range(3)]
        QmI = [[(Q[i][j] - (i == j)) % l for j in range(3)] for i in range(3)]
        out[f"mod{l}"] = {"roots": rts, "N-th powers": [pow(x, N, l) for x in rts], "Q_N": Q,
                          "rank(Q_N - I)": len(rref_pivots(QmI, l))}
    return out


def elim_bound() -> Dict:
    rows = {}
    for N in (15, 105, 3 * 5 * 7 * 11 * 13, 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23):
        ls = [l for l in range(3, N + 1) if N % l == 0 and is_prime(l)]
        phi = 1
        for l in ls:
            phi *= l - 1
        fail1 = Fraction(N - phi, N)
        rows[N] = {"phi/N": float(Fraction(phi, N)), "fail_64_trials_lower_bound": float(fail1 ** 64),
                   "exceeds_2^-64": fail1 ** 64 > Fraction(1, 2 ** 64)}
    # Monte Carlo of the per-trial unit probability for the column [N/l]
    N = 3 * 5 * 7 * 11 * 13
    ls = [3, 5, 7, 11, 13]
    col = [N // l for l in ls]
    rng = random.Random(1)
    T = 200000
    units = sum(gcd(sum(rng.randrange(N) * x for x in col) % N, N) == 1 for _ in range(T))
    res = elim_profile([[x] for x in col], N, seed=3)
    return {"item": "BD-B3-2^-64", "rows": rows, "column_[N/l]_N": N, "unit_rate_MC": units / T,
            "phi/N": 5760 / 15015, "elim_profile_on_column": {k: res[k] for k in ("factor", "profile")}}


def stirling_window(N: int = 391) -> Dict:
    out = []
    for K in range(1, 20):
        v = sum((-1) ** (K - i) * comb(K, i) * ((pow(1 + i, N, N) - pow(i, N, N)) % N) for i in range(K + 1)) % N
        out.append((K, gcd(v, N)))
    return {"item": "BD-Stirling-391", "N": N, "e_p": 23 % 16, "gcds_K=1..19": out}


def model_M_77() -> Dict:
    return {"item": "BD-modelM", "N": 77, "constant": 2 ** 6 - 1, "elim_result": elim_profile([[2 ** 6 - 1]], 77)}


def bm_step_bound(bits=(16, 20), seeds: int = 30) -> Dict:
    """Prop. (repair): e_p >= 3  =>  BM over Z/N on s_c = F_N(c) meets a nonunit discrepancy at
    index <= p + e_p - 1; any nonunit at index >= 2 e_p reveals p."""
    rows = []
    viol = 0
    for b in bits:
        for s in range(seeds):
            inst = I.balanced_semiprime(b, s)
            N = inst.N
            p, q = inst.labels["factors"]
            e = q % (p - 1)
            if e < 3:
                continue
            T = p + e + 5
            seq = [(pow(1 + c, N, N) - 1 - pow(c, N, N)) % N for c in range(T)]
            r = bm_mod_N(seq, N)
            step = r.get("step")
            ok = r.get("factor") is not None and step is not None and step <= p + e - 1
            if step is not None and step >= 2 * e:
                ok &= r["factor"] == p
            viol += not ok
            rows.append((b, s, p, e, step, r.get("factor") == p, step == 2 * e))
    return {"item": "BD-BM-bound", "runs": len(rows), "violations": viol,
            "step==2e": sum(x[6] for x in rows), "revealed_p": sum(x[5] for x in rows),
            "max_step_over_(p+e-1)": max(x[4] / (x[2] + x[3] - 1) for x in rows)}


ALL = [rank_vs_pivot, sporadic_six, n209, n247, elim_bound, stirling_window, model_M_77, bm_step_bound]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="experiments/raw/repair_audit_trackBD.jsonl")
    a = ap.parse_args(argv)
    with open(a.out, "w") as fh:
        for f in ALL:
            r = f()
            fh.write(json.dumps(r, default=str) + "\n")
            print(json.dumps(r, default=str))


if __name__ == "__main__":
    main()
