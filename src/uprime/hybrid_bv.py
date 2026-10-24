"""uprime E1: toy check of the hybrid restricted BV (Lemma 8.5H; notes/uprime_prereg.md).  DIAGNOSTIC.

Values: y2-smooth n <= x (exact, segmented trial division by the primes <= y2 on numpy int64 arrays).
Moduli t = l * s with l prime in (Y, 2Y] (l > y2) and s in (Q/(2l), Q/l], split into families
  F1 = {s y1-smooth}   (the hybrid family: B-smooth modulus cofactors),
  Fall = {all s}        (here s <= Q/Y <= y2, so Fall is also the y2-smooth-cofactor family of Lemma 8.5' at y = y2),
  Frough = Fall \\ F1.
Per modulus: Psi_t (smooth n coprime to t), mu_t = Psi_t/phi(t), the residue histogram over reduced classes,
  maxE(t) = max_{(r,t)=1} |count(r) - mu_t|,  dev(t) = maxE/mu_t,  bad(t) = [dev > 1/2]  (paper2 def:good),
  V(t) = sum_r (count - mu)^2/(phi(t) mu)  (variance ratio against Poisson).
Hybrid normalised sum  H(F) = sum_{t in F} maxE(t) / (Psi(x, y2) #F / Q).
Integer arithmetic for all counts; floats only for the summary ratios (statistics).

Run: .venv/bin/python -m src.uprime.hybrid_bv <out.jsonl>
"""
from __future__ import annotations

import json
import math
import sys
import time
from math import gcd
from typing import Dict, List

import numpy as np

CONFIG = dict(x=2 ** 28, y2=200, y1s=[7, 13], Y=256, Qs=[2 ** 12, 2 ** 13, 2 ** 14], seg=2 ** 23)


def primes_upto(n: int) -> List[int]:
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def smooth_upto(x: int, y: int, seg: int = 2 ** 23) -> np.ndarray:
    """All y-smooth integers in [1, x], exactly (int64)."""
    ps = primes_upto(y)
    out = []
    for lo in range(1, x + 1, seg):
        hi = min(x, lo + seg - 1)
        arr = np.arange(lo, hi + 1, dtype=np.int64)
        for p in ps:
            pk = p
            while pk <= hi:
                start = (-lo) % pk
                arr[start::pk] //= p
                pk *= p
        out.append(np.arange(lo, hi + 1, dtype=np.int64)[arr == 1])
    return np.concatenate(out)


def is_smooth(n: int, y: int) -> bool:
    for p in primes_upto(y):
        while n % p == 0:
            n //= p
    return n == 1


def phi(n: int) -> int:
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def modulus_stats(vals: np.ndarray, t: int) -> Dict:
    h = np.bincount((vals % t).astype(np.int64), minlength=t)
    red = np.array([gcd(r, t) == 1 for r in range(t)])
    cnt = h[red]
    ph = int(red.sum())
    psi_t = int(cnt.sum())
    mu = psi_t / ph
    dev_abs = np.abs(cnt - mu)
    maxE = float(dev_abs.max())
    V = float(((cnt - mu) ** 2).sum() / (ph * mu))
    return {"t": t, "phi": ph, "psi_t": psi_t, "maxE": maxE, "dev": maxE / mu, "bad": bool(maxE / mu > 0.5),
            "V": V, "z": maxE / math.sqrt(2 * mu * math.log(max(ph, 3)))}


def run(out: str, cfg: Dict = CONFIG) -> None:
    t0 = time.time()
    vals = smooth_upto(cfg["x"], cfg["y2"], cfg["seg"])
    psi = int(vals.size)
    ells = [p for p in primes_upto(2 * cfg["Y"]) if p > cfg["Y"]]
    assert min(ells) > cfg["y2"]
    with open(out, "w") as fh:
        fh.write(json.dumps({"header": True, "cfg": cfg, "Psi": psi, "secs_sieve": time.time() - t0}) + "\n")
        for Q in cfg["Qs"]:
            assert Q * Q <= cfg["x"]
            for ell in ells:
                for s in range(max(1, Q // (2 * ell) + 1), Q // ell + 1):
                    t = ell * s
                    if not (Q // 2 < t <= Q):
                        continue
                    st = modulus_stats(vals, t)
                    st.update({"Q": Q, "ell": ell, "s": s,
                               "smooth": {str(y1): is_smooth(s, y1) for y1 in cfg["y1s"]}})
                    fh.write(json.dumps(st) + "\n")
            fh.flush()


def summarize(path: str) -> Dict:
    rows = [json.loads(l) for l in open(path)]
    hdr = rows[0]
    psi = hdr["Psi"]
    cfg = hdr["cfg"]
    rows = rows[1:]
    out = {"Psi": psi, "cells": []}
    for Q in cfg["Qs"]:
        RQ = [r for r in rows if r["Q"] == Q]
        allH = sum(r["maxE"] for r in RQ) / (psi * len(RQ) / Q)
        for y1 in cfg["y1s"]:
            F1 = [r for r in RQ if r["smooth"][str(y1)]]
            Fr_ = [r for r in RQ if not r["smooth"][str(y1)]]
            H1 = sum(r["maxE"] for r in F1) / (psi * len(F1) / Q)
            med = lambda xs: float(np.median(xs)) if xs else float("nan")
            out["cells"].append({
                "Q": Q, "y1": y1, "nF1": len(F1), "nFall": len(RQ), "nFrough": len(Fr_),
                "bad_F1": sum(r["bad"] for r in F1), "bad_all": sum(r["bad"] for r in RQ),
                "H_F1": H1, "H_all": allH, "H_ratio": H1 / allH,
                "medV_F1": med([r["V"] for r in F1]), "medV_all": med([r["V"] for r in RQ]),
                "medV_ratio": med([r["V"] for r in F1]) / med([r["V"] for r in RQ]),
                "medz_F1": med([r["z"] for r in F1]), "medz_all": med([r["z"] for r in RQ]),
                "maxdev_F1": max(r["dev"] for r in F1),
            })
    return out


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "summary":
        print(json.dumps(summarize(sys.argv[2]), indent=1))
    else:
        run(sys.argv[1])
