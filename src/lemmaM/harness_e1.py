"""HARNESS / DIAGNOSTIC for E1 (Dixon descent).  Uses the known factors only to classify
X/Y in mu_2(Z/N)/{+-1} (3-prime N), to build the nontrivial epsilon for the flip check,
and to label instances.  Operational code (descent.py) never sees the factors."""
from __future__ import annotations

import json
import random
import sys
import time
from math import gcd
from typing import Dict, List

from ..nfsfruit.gf2 import left_kernel
from ..nfsfruit.harness import phi_class
from ..nfsfruit.instances import composite
from . import descent as dsc


def crt_eps(fs: List[int], N: int) -> int:
    """epsilon = -1 mod fs[0], +1 mod the others (nontrivial square root of 1)."""
    e = 0
    for i, p in enumerate(fs):
        Mi = N // p
        ri = (p - 1) if i == 0 else 1
        e = (e + ri * Mi * pow(Mi, -1, p)) % N
    return e


def random_kernel_elements(kern: List[int], rng: random.Random, k: int) -> List[int]:
    out = []
    if not kern:
        return out
    for _ in range(k):
        S = 0
        while S == 0:
            for b in kern:
                if rng.random() < 0.5:
                    S ^= b
        out.append(S)
    return out


def run_instance(cls: str, bits: int, seed: int, cfg: Dict, mode: str) -> Dict:
    N, fs = composite(cls, bits, seed)
    p = cfg[str(bits)]
    t0 = time.time()
    rec = dict(cls=cls, bits=bits, seed=seed, mode=mode, N=str(N), omega=len(fs))
    try:
        if mode == "pure":
            sysm = dsc.DescentSystem(N, p["B"], p["B1"], p["J"])
            rng = random.Random(f"lemmaM-desc|{N}|{seed}|{mode}")
            target = sysm.ncols + cfg["excess"]
            while len(sysm.rows) < target:
                if not sysm.add_descent(rng, cfg["max_trials"]):
                    break
        elif mode == "plain":   # Dixon with Farey reconstruction, no descent (B1 = B)
            sysm = dsc.DescentSystem(N, p["B"], p["B"], p["J"])
            rng = random.Random(f"lemmaM-desc|{N}|{seed}|{mode}")
            for _ in range(cfg["plain_n"]):
                if not sysm.add_descent(rng, cfg["max_trials"]):
                    break
            rec.update(nrel=len(sysm.rows), counters=dict(sysm.counters),
                       secs=round(time.time() - t0, 2))
            return rec
        else:  # mixed: plain QS relations up to ncols - n_mixD, then n_mixD descent rels
            sysm = dsc.DescentSystem(N, p["B"], p["B1"], p["J"])
            nq = sysm.ncols - cfg["n_mixD"] + cfg["excess"]
            sysm.add_qs(p["qs_M"], max(0, nq))
            rng = random.Random(f"lemmaM-desc|{N}|{seed}|{mode}")
            for _ in range(cfg["n_mixD"]):
                if not sysm.add_descent(rng, cfg["max_trials"]):
                    break
    except dsc.FactorFound as e:
        rec.update(accidental_factor=True, g=e.g)
        return rec
    rank, kern = left_kernel(sysm.rows)
    rec.update(nrel=len(sysm.rows), ncols=sysm.ncols, nD=sysm.kind.count("D"),
               nQ=sysm.kind.count("Q"), rank=rank, dimK=len(kern), counters=dict(sysm.counters))
    rng2 = random.Random(f"lemmaM-e1-kernel|{N}|{seed}|{mode}")
    Ss = random_kernel_elements(kern, rng2, cfg["kernel_samples"])
    if mode == "mixed":   # also sample the pure-QS dependency space (rows of kind Q only)
        qidx = [i for i, k in enumerate(sysm.kind) if k == "Q"]
        _, qk = left_kernel([sysm.rows[i] for i in qidx])
        qk_full = []
        for b in qk:
            S = 0
            for t, i in enumerate(qidx):
                if (b >> t) & 1:
                    S |= 1 << i
            qk_full.append(S)
        Ss += random_kernel_elements(qk_full, rng2, cfg["kernel_samples"])
        rec["dimK_Q"] = len(qk)
    eps = crt_eps(fs, N)
    res = dict(D=[0, 0], Q=[0, 0], nonsquare=0, flip_ok=0, flip_fail=0, cls_counts={})
    Dmask = 0
    for i, k in enumerate(sysm.kind):
        if k == "D":
            Dmask |= 1 << i
    for S in Ss:
        X, Y = sysm.congruence(S)
        if (X * X - Y * Y) % N != 0:
            res["nonsquare"] += 1
            continue
        g = gcd((X - Y) % N, N)
        fruitful = 1 < g < N
        c = phi_class(X, Y, N, fs)
        res["cls_counts"][str(c)] = res["cls_counts"].get(str(c), 0) + 1
        key = "D" if S & Dmask else "Q"
        res[key][0] += 1
        res[key][1] += int(fruitful)
        if key == "D":
            # flip check: rebuild X with r_j -> eps*r_j for the lowest D-relation j in S.
            # The stored X_j = r_j W_j and the descent output depends on R = r_j^2 only.
            j = (S & Dmask & -(S & Dmask)).bit_length() - 1
            r = sysm.hidden_r[j]
            R2 = (eps * r) ** 2 % N
            again = sysm.split_and_descend(R2)
            W_j = sysm.X[j] * pow(r, -1, N) % N
            same = again is not None and again[0] == W_j and again[1] == sysm.E[j]
            X2 = X * pow(sysm.X[j], -1, N) * (eps * r % N) * again[0] % N if same else None
            c2 = phi_class(X2, Y, N, fs) if same else None
            if same and c is not None and c2 == c ^ 1:
                res["flip_ok"] += 1
            else:
                res["flip_fail"] += 1
    rec.update(res=res, secs=round(time.time() - t0, 2))
    return rec


def _job(a):
    return run_instance(*a)


def main():
    cfg = json.load(open(sys.argv[1]))
    out = sys.argv[2]
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    jobs = [(c, b, s, cfg, m) for m in cfg["modes"] for b in cfg["bits"]
            for c in cfg["classes"] for s in range(cfg["seeds"])]
    from multiprocessing import Pool
    with Pool(nproc) as pool, open(out, "w") as fh:
        for rec in pool.imap_unordered(_job, jobs):
            fh.write(json.dumps(rec) + "\n")
            fh.flush()


if __name__ == "__main__":
    main()
