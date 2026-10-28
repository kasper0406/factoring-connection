"""HARNESS.  Pre-registered experiments (notes/cov38b_prereg.md, addendum A).
Run: PYTHONPATH=. .venv/bin/python -m src.cov38b.experiments e1|e2|e3"""
from __future__ import annotations

import json
import os
import sys
import time
from math import log, sqrt, isqrt

import numpy as np
import sympy

from src.cov38b import design as ds
from src.cov38b import fourier as fo

RAW = "experiments/raw"


def _out(name):
    os.makedirs(RAW, exist_ok=True)
    return os.path.join(RAW, f"cov38b_{name}.jsonl")


def _write(path, rec):
    with open(path, "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(json.dumps({k: v for k, v in rec.items() if k != "hist"}), flush=True)


def qs_e1():
    out = []
    for Q0 in (10 ** 5, 10 ** 6, 10 ** 7, 10 ** 8, 10 ** 9):
        for f in ((1, 1), (13, 10), (17, 10)):
            out.append((Q0, int(sympy.nextprime(Q0 * f[0] // f[1]))))
    return out


def e1():
    path = _out("e1")
    if os.path.exists(path):
        os.remove(path)
    for Q0, q in qs_e1():
        t = time.time()
        r = ds.worst_case(q)
        r.update({"Q0": Q0, "secs": round(time.time() - t, 1)})
        _write(path, r)
        if Q0 <= 10 ** 7:
            t = time.time()
            r = ds.worst_case(q, all_residues=True)
            r.update({"Q0": Q0, "secs": round(time.time() - t, 1)})
            _write(path, r)


def e2():
    path = _out("e2")
    if os.path.exists(path):
        os.remove(path)
    for Q0 in (10 ** 5, 10 ** 6, 10 ** 7):
        q = int(sympy.nextprime(Q0))
        Y = ds.y_of(q)
        L = ds.band(Y)
        J = len(L)
        K = fo.kloosterman_all(q, L)
        kappa = float(K[1:].max())
        _write(path, {"part": "kappa", "q": q, "Y": Y, "J": J, "kappa": kappa,
                      "kappa_over_sqrtJlnq": kappa / sqrt(J * log(q)),
                      "parseval_rms": float(sqrt((K[1:] ** 2).mean()))})
        for Hp in (1, 2, 4, 8, 16, 32):
            k0 = q // (8 * Hp ** 3)
            if k0 < 16:
                continue
            for shape in ("sharp", "tri"):
                st = fo.stats(fo.farey(q, Hp, k0, shape))
                _write(path, {"part": "wiener", "q": q, "Hp": Hp, "k0": k0, "shape": shape, **st})
        # certificate scan at the design
        found = {"farey": None, "arc": None}
        H = 4
        while H <= 2 * Y and (found["farey"] is None or found["arc"] is None):
            k0 = min(q * H // (4 * Y), (q - 1) // 2)
            Hp = max(1, H // 4)
            for shape in ("farey", "arc"):
                if found[shape] is not None:
                    continue
                g = fo.farey(q, Hp if shape == "farey" else 1, k0, "tri")
                st = fo.stats(g)
                c = fo.certificate(J, kappa, st, q)
                _write(path, {"part": "cert", "q": q, "Y": Y, "J": J, "H": H, "Hp": Hp, "k0": k0,
                              "shape": shape, "mass": st["mass"], "W": st["W"], "cert": c})
                if c > 0:
                    found[shape] = H
            H += 4
        _write(path, {"part": "HET", "q": q, "Y": Y, "J": J, "sqrtY": isqrt(Y),
                      "HET_farey": found["farey"], "HET_arc": found["arc"], "Hmax_scanned": H - 4})


def e3():
    from src.dchworst.experiments import sample_pairs
    from src.coverage38.explore import lstar_vec
    path = _out("e3")
    if os.path.exists(path):
        os.remove(path)
    for n in (48, 64, 80):
        P, Q = sample_pairs(n, 200000, seed0=3838 + n, procs=8)
        Pf = np.sqrt(P.astype(np.float64) * Q.astype(np.float64))
        X = Pf ** 0.4
        Bmax = int(2 * X.max() ** 0.5) + 4
        sc = lstar_vec(P, Q, Bmax)
        rec = {"n": n, "pairs": int(len(P)), "max_Lstar_over_P25": float((sc / X).max())}
        for G in (2, 4, 8, 16, 32, 64):
            rec[f"frac_{G}"] = float((sc <= X / G).mean())
        _write(path, rec)


if __name__ == "__main__":
    {"e1": e1, "e2": e2, "e3": e3}[sys.argv[1]]()
