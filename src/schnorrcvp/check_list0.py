"""Certification spot-check (codex schnorrcvp_2 #19): list0 (mpfr enumeration, float target) vs fplll proved CVP.
Each proved call runs in a child process with a 60 s alarm.  Usage: .venv/bin/python -m src.schnorrcvp.check_list0"""
import json
import multiprocessing as mpc
import random

from src.common.instances import balanced_semiprime

from .lattice import build_lattice, exact_cvp, list_cvp, reduce_basis

VARS = [("S21perm", "perm", None), ("S21sqrt", "sqrtperm", None), ("A0.5", "sqrtlnpperm", 0.5),
        ("A1.0", "sqrtlnpperm", 1.0), ("A1.25", "sqrtlnpperm", 1.25)]


def one(args, q):
    bits, s, n, (lab, wt, c), ls = args
    N = balanced_semiprime(bits, s).N
    B, t, d, P, meta = build_lattice(N, n, c_exp=c, weight=wt, schnorr_scale=(c is None), seed=ls)
    reduce_basis(B)
    L0 = list_cvp(B, t, 20)[0]
    pv = exact_cvp(B, t)
    d0 = sum((a - b) ** 2 for a, b in zip(L0, t)); dp = sum((a - b) ** 2 for a, b in zip(pv, t))
    q.put({"bits": bits, "seed": s, "n": n, "variant": lab, "lat_seed": ls, "dist_list0": str(d0), "dist_proved": str(dp),
           "agree": d0 == dp})


if __name__ == "__main__":
    rng = random.Random("schnorrcvp|check_list0")
    out = open("experiments/raw/schnorrcvp_E2_list0check.jsonl", "w")
    res = {"agree": 0, "disagree": 0, "timeout": 0}
    for _ in range(100):
        args = (rng.choice([20, 24, 28, 32, 36, 40]), rng.randrange(3), rng.choice([12, 20, 30]), rng.choice(VARS), rng.randrange(10))
        q = mpc.Queue(); p = mpc.Process(target=one, args=(args, q)); p.start(); p.join(60)
        if p.is_alive():
            p.terminate(); res["timeout"] += 1
            out.write(json.dumps({"args": str(args), "timeout": True}) + "\n"); continue
        r = q.get() if not q.empty() else None
        if r is None:
            res["timeout"] += 1; continue
        res["agree" if r["agree"] else "disagree"] += 1
        out.write(json.dumps(r) + "\n"); out.flush()
    print(res)
