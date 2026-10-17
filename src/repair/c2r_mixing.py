"""C-2R extension X1: additive mixing of two public power branches (agent `repair`).

DIAGNOSTIC code (p known).  Definition frozen in notes/repair_C2R_prereg.md section 6
BEFORE this was run on any instance prime.

    psi(y) = y^E1 + beta * y^E2,   2 <= E1 < E2 <= 64, both even, beta in F_p^*,
    phi(x) = (psi(x + m) + v)^2,   Q(Y) = prod_{w in phi(Z_L)} (Y - w),  d = |phi(Z_L)|.

psi is even, so phi(Z_L) = {(w + v)^2 : w in W}, W = psi([0, m]).  With c = -2v,
d = |W| - #{unordered {w, w'} subset W : w != w', w + w' = c}  (exact count, any c).
Stage 1 computes min_c d for every (E1, E2, beta) exhaustively.  Stage 2 checks exactness
({x in F_p : phi(x) in phi(Z_L)} = Z_L) by evaluating phi on all of F_p for EVERY
(E1, E2, beta, c) with 8d <= 3L.  Success per instance: an exact fold with 8d <= 3L.

Run:  .venv/bin/python -m src.repair.c2r_mixing --out experiments/raw/repair_C2R_mixing.jsonl
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from typing import Dict, List

import numpy as np

from src.repair.c2r_search import instance_primes, params

EXPS = list(range(2, 65, 2))
CHUNK = 512


def _maxrun_rows(S: np.ndarray) -> np.ndarray:
    """Row-wise maximum multiplicity of values in an integer matrix."""
    S = np.sort(S, axis=1)
    n = S.shape[1]
    run = np.ones(S.shape[0], dtype=np.int64)
    k = 2
    alive = np.ones(S.shape[0], dtype=bool)
    while k <= n:
        hit = (S[:, k - 1:] == S[:, : n - k + 1]).any(axis=1) & alive
        if not hit.any():
            break
        run[hit] = k
        alive = hit
        k += 1
    return run


def exact_d_rowwise(W_row: np.ndarray, p: int) -> List[tuple]:
    """Exact (d, c) list minimising d for one row (handles repeated values in W)."""
    Wd = np.unique(W_row)
    best = {}
    for a, b in combinations(Wd.tolist(), 2):
        c = (a + b) % p
        best[c] = best.get(c, 0) + 1
    if not best:
        return [(len(Wd), None)]
    mx = max(best.values())
    return [(len(Wd) - mx, c) for c, k in best.items() if k == mx]


def is_exact(p: int, m: int, E1: int, E2: int, beta: int, c: int) -> bool:
    L = 2 * m + 1
    v = (-c * pow(2, -1, p)) % p
    y = (np.arange(p, dtype=object) + m) % p
    psi = np.array([(pow(int(t), E1, p) + beta * pow(int(t), E2, p)) % p for t in y], dtype=np.int64)
    phi = ((psi + v) % p) ** 2 % p
    ZL = np.array([(-i) % p for i in range(L)])
    R = np.zeros(p, dtype=bool)
    R[phi[ZL]] = True
    zero = np.nonzero(R[phi])[0]
    return set(zero.tolist()) == set(ZL.tolist())


def search_prime(p: int, cpu_cap: float = 1800.0, exps: List[int] = EXPS) -> Dict:
    t0 = time.process_time()
    m, L = params(p)
    j = np.arange(m + 1, dtype=np.int64)
    iu = np.triu_indices(m + 1, k=1)
    betas = np.arange(1, p, dtype=np.int64)
    min_d = (m + 2, None)
    candidates = []  # (E1, E2, beta, c, d) with 8d <= 3L
    n_rows = n_repeat_rows = 0
    for E1, E2 in combinations(exps, 2):
        if time.process_time() - t0 > cpu_cap:
            return {"p": p, "m": m, "L": L, "status": "inconclusive: CPU cap",
                    "cpu_seconds": time.process_time() - t0}
        a1 = np.array([pow(int(x), E1, p) for x in j], dtype=np.int64)
        a2 = np.array([pow(int(x), E2, p) for x in j], dtype=np.int64)
        s1 = a1[iu[0]] + a1[iu[1]]
        s2 = a2[iu[0]] + a2[iu[1]]
        for k0 in range(0, len(betas), CHUNK):
            b = betas[k0:k0 + CHUNK]
            W = (a1[None, :] + b[:, None] * a2[None, :]) % p
            Ws = np.sort(W, axis=1)
            rep = (Ws[:, 1:] == Ws[:, :-1]).any(axis=1)
            S = (s1[None, :] + b[:, None] * s2[None, :]) % p
            run = _maxrun_rows(S)
            d = (m + 1) - run
            n_rows += len(b)
            for i in np.nonzero(rep)[0]:  # exact handling of rows with repeated W values
                n_repeat_rows += 1
                lst = exact_d_rowwise(W[i], p)
                d[i] = lst[0][0]
            i = int(d.argmin())
            if d[i] < min_d[0]:
                min_d = (int(d[i]), (E1, E2, int(b[i])))
            for i in np.nonzero(8 * d <= 3 * L)[0]:
                for dd, c in exact_d_rowwise(W[i], p):
                    candidates.append((E1, E2, int(b[i]), c, dd))
    exact_hits = [cand for cand in candidates
                  if cand[3] is not None and is_exact(p, m, cand[0], cand[1], cand[2], cand[3])]
    return {"p": p, "m": m, "L": L, "status": "ok", "experiment": "C2R_X1_mixing",
            "min_d_any_c": min_d[0], "argmin_E1_E2_beta": min_d[1],
            "n_candidates_8d_le_3L": len(candidates), "n_exact_hits": len(exact_hits),
            "exact_hits": exact_hits[:20], "pass": len(exact_hits) > 0,
            "n_rows": n_rows, "n_repeat_rows": n_repeat_rows,
            "cpu_seconds": time.process_time() - t0}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="experiments/raw/repair_C2R_mixing.jsonl")
    ap.add_argument("--procs", type=int, default=32)
    ap.add_argument("--cpu-cap", type=float, default=1800.0)
    ap.add_argument("--primes", type=int, nargs="*")
    a = ap.parse_args(argv)
    primes = a.primes or instance_primes()
    with ProcessPoolExecutor(max_workers=a.procs) as pool, open(a.out, "a") as fh:
        futs = [pool.submit(search_prime, p, a.cpu_cap) for p in primes]
        for f in futs:
            row = f.result()
            fh.write(json.dumps(row) + "\n")
            fh.flush()
            print(json.dumps({k: row.get(k) for k in ("p", "L", "min_d_any_c", "pass", "cpu_seconds")}),
                  flush=True)


if __name__ == "__main__":
    main()
