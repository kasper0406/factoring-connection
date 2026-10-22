"""dsr harness (DIAGNOSTIC). E1 uses the factorisation (discrete-log tables) and float64 FFTs; E2's random controls
use seeded `random`. Never imported by src/dsr/core.py.

Usage:
  PYTHONPATH=. .venv/bin/python -m src.dsr.harness e1 experiments/raw/dsr_e1.jsonl
  PYTHONPATH=. .venv/bin/python -m src.dsr.harness e2 experiments/raw/dsr_e2.jsonl
"""
from __future__ import annotations

import json
import math
import random
import resource
import sys
import time
from typing import Dict, List

import numpy as np
from sympy import primitive_root

from src.common.instances import balanced_semiprime, prime_instance
from src.dsr import core

A_CAN, A_BAL = 0.8, 0.5


def _rss_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20  # macOS: bytes


# ---------------------------------------------------------------- E1 (Fourier-algebra norm)

def _dlog_table(p: int) -> np.ndarray:
    g = primitive_root(p)
    ind = np.zeros(p, dtype=np.int64)
    x = 1
    for k in range(p - 1):
        ind[x] = k
        x = x * g % p
    return ind


def group_array(mask: np.ndarray, factors: List[int]) -> np.ndarray:
    """Indicator of mask (over [0, M)) on G = prod (Z/p)^x, indexed by discrete logs (squarefree M only)."""
    M = len(mask)
    xs = np.nonzero(mask)[0].astype(np.int64)
    if len(factors) == 1:
        p = factors[0]
        arr = np.zeros(p - 1)
        arr[_dlog_table(p)[xs % p]] = 1.0
        return arr
    p, q = factors
    arr = np.zeros((p - 1, q - 1))
    arr[_dlog_table(p)[xs % p], _dlog_table(q)[xs % q]] = 1.0
    return arr


def fourier_stats(arr: np.ndarray) -> Dict:
    F = np.abs(np.fft.fftn(arr))
    phi = arr.size
    s = float(arr.sum())
    pars = float((F ** 2).sum()) / (phi * s)
    F.flat[0] = 0.0
    return {"A1": float(F.sum()) / phi, "beta": float(F.max()), "parseval_rel_err": abs(pars - 1.0)}


def run_e1(out) -> None:
    mods = [("prime", b, s) for b in (16, 18, 20, 22) for s in range(3)] + \
           [("semi", b, s) for b in (20, 22) for s in range(3)]
    for fam, bits, seed in mods:
        inst = prime_instance(bits, seed) if fam == "prime" else balanced_semiprime(bits, seed)
        M, fac = inst.N, inst.factors
        G = core.units_mask(M)
        phi = int(G.sum())
        for name, a in (("S_can", A_CAN), ("S_bal", A_BAL)):
            t0 = time.time()
            B, P = core.factor_base(M, a)
            S = (core.s_can(M, P) if name == "S_can" else core.s_bal(M, P)) & G
            size = int(S.sum())
            st = fourier_stats(group_array(S, fac))
            ctrl = []
            for k in range(3):
                rng = np.random.default_rng(1000 * bits + seed + k)
                units = np.nonzero(G)[0]
                R = np.zeros(M, dtype=bool)
                R[rng.choice(units, size=size, replace=False)] = True
                ctrl.append(fourier_stats(group_array(R, fac)))
            row = {"exp": "E1", "fam": fam, "bits": bits, "seed": seed, "M": M, "phi": phi, "set": name, "a": a,
                   "B": B, "nP": len(P), "size": size, "delta": size / phi, **st,
                   "ctrl_A1": [c["A1"] for c in ctrl], "ctrl_beta": [c["beta"] for c in ctrl],
                   "ctrl_parseval": max(c["parseval_rel_err"] for c in ctrl),
                   "secs": round(time.time() - t0, 2), "rss_mb": round(_rss_mb(), 1)}
            out.write(json.dumps(row) + "\n")
            out.flush()
            print(fam, bits, seed, name, "delta=%.3f rA=%.3f rinf=%.3f" % (
                row["delta"], st["A1"] / np.mean(row["ctrl_A1"]), st["beta"] / np.mean(row["ctrl_beta"])),
                flush=True)


# ---------------------------------------------------------------- E2 (worst-case coverage)

def r_g(N: int, rng: random.Random):
    while True:
        x = rng.randrange(1, N)
        if math.gcd(x, N) == 1:
            yield x


def r_walk(P, t: int, rng: random.Random):
    while True:
        w = 1
        for _ in range(t):
            w *= rng.choice(P)
        yield w


def run_e2(out) -> None:
    for bits in (20, 22, 24):
        for seed in range(10):
            inst = balanced_semiprime(bits, seed)
            N = inst.N
            G = core.units_mask(N)
            phi = int(G.sum())
            H = 8 * math.ceil(math.log2(N))
            for name, a, lists in (("S_can", A_CAN, ("R_G", "R_walk", "D_box", "D_slide")),
                                   ("S_bal", A_BAL, ("R_G", "D_box"))):
                B, P = core.factor_base(N, a)
                S = (core.s_can(N, P) if name == "S_can" else core.s_bal(N, P)) & G
                delta = int(S.sum()) / phi
                lam = -math.log1p(-delta)
                cap = math.ceil(50 * math.log(phi) / delta)
                m0 = math.ceil(math.log(phi) / delta)
                t = core.walk_length(N, B)
                targets = [pow(q, H, N) for q in P]

                def make(lst: str, tag: str):
                    rng = random.Random(f"dsr|E2|{lst}|{tag}|{bits}|{seed}|{name}")
                    if lst == "R_G":
                        return r_g(N, rng)
                    if lst == "R_walk":
                        return r_walk(P, t, rng)
                    if lst == "D_box":
                        return core.d_box(P)
                    return core.d_slide(P, t)

                for lst in lists:
                    t0 = time.time()
                    cov = core.coverage(N, S, G, make(lst, "cov"), cap, checkpoints=(m0, 2 * m0, 4 * m0))
                    fh = core.first_hits(N, S, targets, make(lst, "fh"), cap)
                    vals = [v for v in fh.values()]
                    row = {"exp": "E2", "bits": bits, "seed": seed, "N": N, "phi": phi, "set": name, "a": a, "B": B,
                           "nP": len(P), "t": t, "delta": delta, "cap": cap, "list": lst,
                           "m_star": cov["m_star"],
                           "nu": (cov["m_star"] * lam / math.log(phi)) if cov["m_star"] else None,
                           "uncovered_at": {str(k): v for k, v in cov["uncovered_at"].items()},
                           "uncovered_final": cov["uncovered_final"],
                           "T_star": (max(vals) if all(v is not None for v in vals) else None),
                           "target_hits_mean": float(np.mean([v for v in vals if v is not None])) if any(
                               v is not None for v in vals) else None,
                           "targets_unhit": sum(v is None for v in vals),
                           "secs": round(time.time() - t0, 2), "rss_mb": round(_rss_mb(), 1)}
                    out.write(json.dumps(row) + "\n")
                    out.flush()
                    print(bits, seed, name, lst, "delta=%.4f m*=%s nu=%s T*=%s" % (
                        delta, row["m_star"], None if row["nu"] is None else round(row["nu"], 3), row["T_star"]),
                        flush=True)
                if name == "S_can":
                    ob = core.orbit_first_hits(N, S, P, H, cap)
                    ov = list(ob.values())
                    row = {"exp": "E2", "bits": bits, "seed": seed, "N": N, "phi": phi, "set": name, "a": a, "B": B,
                           "nP": len(P), "delta": delta, "cap": cap, "list": "D_orbit",
                           "T_star": (max(ov) + 1 if all(v is not None for v in ov) else None),
                           "target_hits_mean": float(np.mean([v + 1 for v in ov if v is not None])),
                           "targets_unhit": sum(v is None for v in ov)}
                    out.write(json.dumps(row) + "\n")
                    out.flush()


def main() -> None:
    exp, path = sys.argv[1], sys.argv[2]
    with open(path, "w") as out:
        {"e1": run_e1, "e2": run_e2}[exp](out)


if __name__ == "__main__":
    main()
