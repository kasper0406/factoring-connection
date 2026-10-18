"""sieve2d pre-registered experiments (notes/sieve2d_prereg.md).

  PYTHONPATH=. .venv/bin/python -m src.sieve2d.experiments exh     # E1 exhaustive n=20..30
  PYTHONPATH=. .venv/bin/python -m src.sieve2d.experiments samp    # E2/E3 sampled n=48..96
  PYTHONPATH=. .venv/bin/python -m src.sieve2d.experiments summarize

HARNESS code (uses p, q).  Raw JSONL: experiments/raw/sieve2d_{exh,samp}.jsonl.
"""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from math import exp, log
from typing import Dict, List, Sequence, Tuple

import numpy as np

from src.dcltheorem.params import icbrt, pow2_frac, primes_upto, weight_set

from .coverage import (BandTables, class_table, degenerate_mask, exhaustive_integer_pairs,
                       exhaustive_prime_pairs, inverse_table, lam_nominal, model_integer,
                       model_prime)
from .sampling import sample_integer_pairs, sample_prime_pairs

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "experiments", "raw")
CFG = os.path.join(ROOT, "experiments", "configs", "sieve2d.json")

CONFIG = {
    "exh_n": [20, 22, 24, 26, 28, 30],
    "exh_int_n": [20, 22, 24, 26],
    "exh_H": list(range(2, 13)),
    "sep_H": list(range(2, 11)),
    "samp_n": [48, 64, 80, 96],
    "samp_pairs": 1000000,
    "samp_seed": 20261018,
    "samp_H": [2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 28, 32],
    "samp_blocks": {"48": [0, 1, 2, 3], "64": [0, 1, 2, 3], "80": [0], "96": [0]},
}


def thm_Y(n: int, H: int) -> int:
    """Theorem-shaped band base Y = ceil(floor(2^(n/6)) / floor(H^(1/3)))
    (src/dcltheorem/params.theorem_design convention)."""
    return max(3, -(-pow2_frac(n, 1, 6) // max(1, icbrt(H))))


def band(Y: int, k: int) -> List[int]:
    lo, hi = (1 << k) * Y, (1 << (k + 1)) * Y
    return [l for l in primes_upto(hi) if lo < l <= hi and l % 4 == 3]


def eval_cells(p: np.ndarray, q: np.ndarray, bandl: Sequence[int],
               Vs: Sequence[Sequence[Tuple[int, int]]], want_M: bool = True):
    """For each weight set V in Vs: M(p,q) array (cover multiplicity) over the band."""
    Ms = [np.zeros(p.shape, dtype=np.int64) for _ in Vs]
    for l in bandl:
        inv = inverse_table(l)
        pm = p % l
        qm = q % l
        qz = qm == 0
        r = np.where(qz, l, (pm * inv[qm]) % l)
        both = qz & (pm == 0)
        anyboth = bool(both.any())
        for i, V in enumerate(Vs):
            C = class_table(l, V)
            c = C[r]
            if anyboth:
                c = np.where(both, len(V), c)
            Ms[i] += c
    return Ms


def model_cell(bandl, V, kind: str) -> Dict:
    tabs = [BandTables(l=l, C=class_table(l, V), inv=None, nV=len(V)) for l in bandl]
    return model_prime(tabs) if kind == "prime" else model_integer(tabs)


def cell_record(n, kind, fam, k, H, Y, bandl, V, M, ndeg=None) -> Dict:
    mod = model_cell(bandl, V, kind)
    lam = lam_nominal(bandl, len(V))
    npairs = int(M.size)
    unc = int((M == 0).sum())
    Mf = M.astype(float)
    mu = Mf.mean()
    var = Mf.var()
    m4 = ((Mf - mu) ** 4).mean()
    # theorem statistic: mean of (1 - M/lam)^4  (Q of the note)
    Qstat = float((((lam - Mf) / lam) ** 4).mean()) if lam > 0 else None
    Qmodel = None
    if lam > 0:
        # E_model (lam - M)^4 / lam^4 with model mean mu_m (== lam for separated bands)
        mm = mod["mu"]
        dl = lam - mm
        # E(lam - M)^4 = dl^4 + 6 dl^2 var - 4 dl m3 + m4   (m3, m4 central, model exact)
        Qmodel = (dl ** 4 + 6 * dl * dl * mod["var"] - 4 * dl * mod["m3"] + mod["m4"]) / lam ** 4
    rec = {"n": n, "kind": kind, "family": fam, "block": k, "H": H, "Y": Y,
           "band_lo": bandl[0] if bandl else None, "band_hi": bandl[-1] if bandl else None,
           "K": len(bandl), "nV": len(V), "separated": bool(bandl) and bandl[0] > H * H,
           "lam": lam, "lam_eff": -log(mod["f_model"]) if mod["f_model"] > 0 else None,
           "f_model": mod["f_model"], "pairs": npairs, "uncovered": unc,
           "E_model": npairs * mod["f_model"],
           "mean_M": mu, "var_M": var, "m4_M": m4,
           "model_mu": mod["mu"], "model_var": mod["var"], "model_m4": mod["m4"],
           "Qstat": Qstat, "Qmodel": Qmodel}
    if ndeg is not None:
        rec["degenerate_removed"] = ndeg
    return rec


# ---------------------------------------------------------------------------
def _exh_job(args):
    n, kind = args
    t0 = time.time()
    if kind == "prime":
        p, q = exhaustive_prime_pairs(n)
        ndeg = None
    else:
        p, q = exhaustive_integer_pairs(n)
    out = []
    # family thm: theorem-shaped Y, all 4 blocks
    byY: Dict[int, List[int]] = {}
    for H in CONFIG["exh_H"]:
        byY.setdefault(thm_Y(n, H), []).append(H)
    jobs = [("thm", Y, Hs, k) for Y, Hs in byY.items() for k in range(4)]
    # family sep: Y = H^2 + 1 (smallest separated band), blocks 0, 1
    for H in CONFIG["sep_H"]:
        for k in (0, 1):
            jobs.append(("sep", H * H + 1, [H], k))
    for fam, Y, Hs, k in jobs:
        bl = band(Y, k)
        if not bl:
            continue
        Vs = [weight_set(H) for H in Hs]
        if kind == "int":
            Ms = eval_cells(p, q, bl, Vs)
            for H, V, M in zip(Hs, Vs, Ms):
                dm = degenerate_mask(p, q, V)
                out.append(cell_record(n, kind, fam, k, H, Y, bl, V, M[~dm], int(dm.sum())))
        else:
            Ms = eval_cells(p, q, bl, Vs)
            for H, V, M in zip(Hs, Vs, Ms):
                rec = cell_record(n, kind, fam, k, H, Y, bl, V, M)
                if n == 30 and fam == "thm" and k == 0:
                    rec["qdisp"] = q_dispersion(q, M)
                out.append(rec)
    for r in out:
        r["secs_job"] = round(time.time() - t0, 1)
    return out


def q_dispersion(q: np.ndarray, M: np.ndarray) -> Dict:
    """P6: dispersion of the per-q uncovered fraction (q with >= 200 partners)."""
    unc = (M == 0)
    qs, inv, cnt = np.unique(q, return_inverse=True, return_counts=True)
    u = np.bincount(inv, weights=unc.astype(float))
    sel = cnt >= 200
    if sel.sum() < 10:
        return {"nq": int(sel.sum())}
    f = unc.mean()
    if f <= 0 or f >= 1:
        return {"nq": int(sel.sum()), "f": f}
    z2 = (u[sel] - cnt[sel] * f) ** 2 / (cnt[sel] * f * (1 - f))
    fr = u[sel] / cnt[sel]
    return {"nq": int(sel.sum()), "f": f, "D": float(z2.mean()),
            "fq_p01": float(np.quantile(fr, 0.01)), "fq_p99": float(np.quantile(fr, 0.99)),
            "fq_max": float(fr.max())}


def run_exh():
    jobs = [(n, "prime") for n in CONFIG["exh_n"]] + [(n, "int") for n in CONFIG["exh_int_n"]]
    path = os.path.join(RAW, "sieve2d_exh.jsonl")
    with ProcessPoolExecutor(max_workers=10) as ex, open(path, "w") as fh:
        for res in ex.map(_exh_job, jobs):
            for r in res:
                fh.write(json.dumps(r) + "\n")
            fh.flush()
            print(res[0]["n"], res[0]["kind"], len(res), "cells", flush=True)


# ---------------------------------------------------------------------------
def _samp_job(args):
    n, kind, Y, Hs, k = args
    t0 = time.time()
    cnt = CONFIG["samp_pairs"]
    if kind == "prime":
        p, q = sample_prime_pairs(n, cnt, CONFIG["samp_seed"])
    else:
        p, q = sample_integer_pairs(n, cnt, CONFIG["samp_seed"])
    bl = band(Y, k)
    Vs = [weight_set(H) for H in Hs]
    out = []
    if kind == "int":
        Ms = eval_cells(p, q, bl, Vs)
        for H, V, M in zip(Hs, Vs, Ms):
            dm = degenerate_mask(p, q, V)
            out.append(cell_record(n, kind, "thm", k, H, Y, bl, V, M[~dm], int(dm.sum())))
    else:
        Ms = eval_cells(p, q, bl, Vs)
        for H, V, M in zip(Hs, Vs, Ms):
            out.append(cell_record(n, kind, "thm", k, H, Y, bl, V, M))
    for r in out:
        r["secs_job"] = round(time.time() - t0, 1)
    return out


def run_samp():
    jobs = []
    for n in CONFIG["samp_n"]:
        byY: Dict[int, List[int]] = {}
        for H in CONFIG["samp_H"]:
            Y = thm_Y(n, H)
            if Y > H * H:                      # theorem regime only (separated classes)
                byY.setdefault(Y, []).append(H)
        for Y, Hs in byY.items():
            for k in CONFIG["samp_blocks"][str(n)]:
                for kind in ("prime", "int"):
                    jobs.append((n, kind, Y, Hs, k))
    path = os.path.join(RAW, "sieve2d_samp.jsonl")
    with ProcessPoolExecutor(max_workers=12) as ex, open(path, "w") as fh:
        for res in ex.map(_samp_job, jobs):
            for r in res:
                fh.write(json.dumps(r) + "\n")
            fh.flush()
            r = res[0]
            print(r["n"], r["kind"], "Y", r["Y"], "k", r["block"], [x["H"] for x in res],
                  r["secs_job"], flush=True)


def main():
    os.makedirs(RAW, exist_ok=True)
    cmd = sys.argv[1] if len(sys.argv) > 1 else "config"
    if cmd == "config":
        with open(CFG, "w") as fh:
            json.dump(CONFIG, fh, indent=1)
        print(json.dumps(CONFIG))
    elif cmd == "exh":
        run_exh()
    elif cmd == "samp":
        run_samp()
    elif cmd == "summarize":
        from .analyze import summarize
        summarize()


if __name__ == "__main__":
    main()
