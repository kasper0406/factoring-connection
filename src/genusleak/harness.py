"""genusleak HARNESS: ground truth (p, q, class numbers via PARI) and the experiment drivers E1-E5.

Never imported by operational code.  PARI's class-group routines factor D internally: harness only.
"""
from __future__ import annotations

import json
import math
import os
import random
import sys
import time
from fractions import Fraction
from math import gcd, isqrt
from typing import Dict, List

import gmpy2

from src.common.instances import balanced_semiprime

from . import cayley, cf, forms, sqrtred

RAW = os.path.join(os.path.dirname(__file__), "..", "..", "experiments", "raw")


def _pari():
    import cypari2
    return cypari2.Pari()


def instances(n: int, count: int, mod4: int = None):
    """balanced_semiprime(n, seed) for seeds 0, 1, ...; skip those failing N % 4 == mod4 (recorded)."""
    out, skipped, seed = [], [], 0
    while len(out) < count:
        inst = balanced_semiprime(n, seed)
        if mod4 is None or inst.N % 4 == mod4:
            out.append(inst)
        else:
            skipped.append(seed)
        seed += 1
    return out, skipped


def class_number(D: int) -> int:
    return int(_pari().qfbclassno(D))


def _write(name: str, recs: List[Dict]):
    os.makedirs(RAW, exist_ok=True)
    with open(os.path.join(RAW, name), "w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")


# ---------------------------------------------------------------- E1
def e1_instance(inst, tmax: int = 40) -> Dict:
    N = inst.N
    p, q = inst.factors
    D, ells, gens, elems, idx, table = cayley.setup(N, 1, 6)
    H = len(elems)
    W = cayley.closed_walks(table, tmax)
    Wf = cayley.free_closed_walks(len(ells), tmax)
    g0 = cayley.girth_bound(N, 1, ells)
    first_excess = next((t for t in range(tmax + 1) if W[t] != Wf[t]), None)
    m = 2 * len(ells)
    sp = 2 * sum(gmpy2.legendre(l, p) for l in ells)
    sq = 2 * sum(gmpy2.legendre(l, q) for l in ells)
    # characters as functions on the generators (each value counted for f and f^{-1})
    vp = tuple(int(gmpy2.legendre(l, p)) for l in ells)
    vq = tuple(int(gmpy2.legendre(l, q)) for l in ells)
    v4 = tuple(1 if l % 4 == 1 else -1 for l in ells)
    triv = tuple(1 for _ in ells)
    # deviation D1 (disclosed): a hidden genus character that coincides on S (hence on H) with a public character
    # (trivial or chi_-4) or with the other hidden one is not a separate character of H; count it once / not at all.
    hidden_chars = []
    for v in (vp, vq):
        if v not in (triv, v4) and v not in hidden_chars:
            hidden_chars.append(v)
    sums_reg = [sp, sq]
    sums_cor = [2 * sum(v) for v in hidden_chars]

    def mins(sums):
        Mret, Mcoll = [], []
        for t in range(2, tmax + 1, 2):
            P = Fraction(W[t], m ** t)
            G = Fraction(sum(x ** t for x in sums), H * m ** t)
            Mret.append(float(P / (G * G)) if G else math.inf)
        for t in range(1, tmax // 2 + 1):
            P2 = Fraction(W[2 * t], m ** (2 * t))
            G2 = Fraction(sum(x ** (2 * t) for x in sums), H * m ** (2 * t))
            Mcoll.append(math.sqrt(2 * float(P2)) / float(G2) if G2 else math.inf)
        return min(Mret) / H, min(Mcoll) / math.sqrt(H)

    reg = mins(sums_reg)
    cor = mins(sums_cor) if sums_cor else (math.inf, math.inf)
    h = class_number(D)
    hidden = [f for f in elems if forms.is_ambiguous(f) and forms.ambiguous_split(N, f)]
    return {"n": N.bit_length(), "seed": inst.seed, "N": N, "ells": ells, "H": H, "h": h,
            "H_contains_hidden_ambiguous": bool(hidden), "g0": g0, "first_excess": first_excess,
            "P1_ok": all(W[t] == Wf[t] for t in range(min(g0, tmax + 1))),
            "lam_p": sp / m, "lam_q": sq / m,
            "n_hidden_chars_on_H": len(hidden_chars),
            "min_Mret_over_H": reg[0], "min_Mcoll_over_sqrtH": reg[1],
            "min_Mret_over_H_D1": cor[0], "min_Mcoll_over_sqrtH_D1": cor[1]}


def run_e1():
    recs = []
    for n in (20, 24, 28):
        insts, skipped = instances(n, 20, mod4=1)
        for inst in insts:
            t0 = time.time()
            r = e1_instance(inst)
            r["skipped_seeds_before"] = [s for s in skipped if s < inst.seed]
            r["secs"] = round(time.time() - t0, 2)
            recs.append(r)
            print("E1", n, inst.seed, r["H"], r["h"], r["g0"], r["first_excess"], r["P1_ok"],
                  round(r["min_Mret_over_H"], 3), round(r["min_Mcoll_over_sqrtH"], 3), flush=True)
    _write("genusleak_E1.jsonl", recs)
    return recs


# ---------------------------------------------------------------- E2
def run_e2():
    recs = []
    for n in (28, 32, 36):
        insts, _ = instances(n, 30)
        for inst in insts:
            t0 = time.time()
            s = cf.squfof_scan(inst.N, 1)
            L, pr = s["L"], s["proper"]
            rho = len(pr) / L
            rec = {"n": n, "seed": inst.seed, "N": inst.N, "L": L, "n_proper": len(pr), "n_improper": len(s["improper"]),
                   "first": s["first"], "rev_steps_first": s["rev_steps_first"], "rho": rho,
                   "first_times_rho": (s["first"] * rho) if s["first"] else None,
                   "secs": round(time.time() - t0, 2)}
            p, q = inst.factors
            recs.append(rec)
            print("E2", n, inst.seed, L, len(pr), len(s["improper"]), s["first"], flush=True)
    _write("genusleak_E2.jsonl", recs)
    return recs


# ---------------------------------------------------------------- E3
def run_e3():
    recs = []
    for n in (24, 32):
        insts, skipped = instances(n, 10, mod4=1)
        for inst in insts:
            N = inst.N
            D = -4 * N
            ells = forms.split_primes(D, 10)
            gens = [forms.prime_form(D, l) for l in ells]
            t0 = time.time()
            elems, idx, table = cayley.enumerate_subgroup(D, gens)
            roots: Dict = {}
            for z in elems:
                y = forms.compose(z, z)
                if y not in roots or z < roots[y]:
                    roots[y] = z  # canonical: lexicographically least root (deterministic, x-blind)
            oracle = lambda y: roots[y]
            h = class_number(D)
            rng = random.Random(f"genusleak-E3|{n}|{inst.seed}")
            succ, correct = 0, 0
            for _ in range(40):
                out = sqrtred.sqrt_reduction_query(N, D, gens, oracle, rng)
                if out["split"]:
                    succ += 1
                    correct += int(out["split"] in inst.factors)
            rec = {"n": n, "seed": inst.seed, "N": N, "H": len(elems), "h": h, "H_eq_Cl": len(elems) == h,
                   "queries": 40, "success": succ, "correct": correct, "secs": round(time.time() - t0, 2)}
            recs.append(rec)
            print("E3", rec, flush=True)
    _write("genusleak_E3.jsonl", recs)
    return recs


# ---------------------------------------------------------------- E4
MULTS = (1, 2, 3, 5, 6, 7, 10, 11)


def run_e4():
    recs = []
    for n, cnt in ((28, 100), (32, 50)):
        insts, _ = instances(n, cnt)
        for inst in insts:
            row = {"n": n, "seed": inst.seed, "N": inst.N, "by_k": {}}
            for k in MULTS:
                r = cf.midpoint_split(inst.N, k)
                ok = r["split"] is not None
                if ok:
                    assert r["split"] in inst.factors
                row["by_k"][str(k)] = {"L": r["L"], "parity": r["parity"], "split": ok}
            row["k1"] = row["by_k"]["1"]["split"]
            row["any"] = any(v["split"] for v in row["by_k"].values())
            recs.append(row)
            print("E4", n, inst.seed, row["k1"], row["any"], flush=True)
    _write("genusleak_E4.jsonl", recs)
    return recs


# ---------------------------------------------------------------- E5
def run_e5():
    import mpmath as mp
    mp.mp.dps = 60
    xs = [2 ** j for j in range(8, 21)]
    primes = forms.small_primes(xs[-1])
    recs = []
    for n in (32, 40, 48):
        insts, _ = instances(n, 20)
        for inst in insts:
            D = -4 * inst.N
            h = class_number(D)
            prod = mp.mpf(1)
            out, pi_count, j = [], 0, 0
            for x in xs:
                while j < len(primes) and primes[j] <= x:
                    l = primes[j]
                    chi = 0 if l == 2 else int(gmpy2.jacobi(D % l, l))
                    prod *= 1 / (1 - mp.mpf(chi) / l)
                    j += 1
                hhat = mp.sqrt(-D) / mp.pi * prod
                delta = abs(h - hhat)
                out.append({"x": x, "pi": j, "delta": float(delta), "cost": float(j + 2 * mp.sqrt(delta))})
            lx = [math.log(o["x"]) for o in out]
            ld = [math.log(max(o["delta"], 1e-30)) for o in out]
            mx, md = sum(lx) / len(lx), sum(ld) / len(ld)
            slope = sum((a - mx) * (b - md) for a, b in zip(lx, ld)) / sum((a - mx) ** 2 for a in lx)
            cmin = min(o["cost"] for o in out)
            rec = {"n": n, "seed": inst.seed, "N": inst.N, "h": h, "grid": out, "slope": slope,
                   "logC_over_logD": math.log(cmin) / math.log(-D)}
            recs.append(rec)
            print("E5", n, inst.seed, h, round(slope, 3), round(rec["logC_over_logD"], 3), flush=True)
    _write("genusleak_E5.jsonl", recs)
    return recs


if __name__ == "__main__":
    which = sys.argv[1]
    {"E1": run_e1, "E2": run_e2, "E3": run_e3, "E4": run_e4, "E5": run_e5}[which]()
