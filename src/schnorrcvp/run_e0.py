"""E0 (pre-registered): H_S model total fac-relation counts E_tot(bits, n) and n_min(bits).  MODEL numerics.

Usage: .venv/bin/python -m src.schnorrcvp.run_e0  -> experiments/raw/schnorrcvp_E0.jsonl
"""
from __future__ import annotations

import json
import time

import mpmath as mp

from .core import first_primes
from .psi import expected_relations

BITS = [64, 128, 200, 256, 400, 512, 800, 1024]


def etot(bits, n, K=50, lnv_max=None, h=0.1):
    lnN = bits * mp.log(2)
    if lnv_max is None:
        lnv_max = 1.5 * float(lnN) + 50
    r = expected_relations(first_primes(n), lnN, K=K, h=h, lnv_max=lnv_max)
    prof = r["profile"]
    peak = max(prof, key=lambda t: t[1])[0] if prof else None
    return r["E"], peak


def grid_n(lo=5):
    n = lo
    while True:
        yield n
        n += 5 if n < 200 else 25


def main():
    out = open("experiments/raw/schnorrcvp_E0.jsonl", "w")
    # sensitivity checks at Schnorr's points
    for bits, n in [(400, 47), (800, 95)]:
        e50, pk = etot(bits, n, K=50)
        e100, _ = etot(bits, n, K=100)
        lnN = bits * float(mp.log(2))
        e_wide, _ = etot(bits, n, K=50, lnv_max=3 * lnN)
        rec = {"kind": "sensitivity", "bits": bits, "n": n, "log10E_K50": float(mp.log10(e50)),
               "log10E_K100": float(mp.log10(e100)), "log10E_wide": float(mp.log10(e_wide)), "peak_lnv": pk}
        out.write(json.dumps(rec) + "\n"); out.flush(); print(rec, flush=True)
    for bits in BITS:
        t0 = time.time()
        last = None
        for n in grid_n():
            e, pk = etot(bits, n)
            rec = {"kind": "scan", "bits": bits, "n": n, "pn": first_primes(n)[-1],
                   "log10E": float(mp.log10(e)) if e > 0 else None, "peak_lnv": pk}
            out.write(json.dumps(rec) + "\n"); out.flush()
            if e >= n + 1:
                rec2 = {"kind": "nmin", "bits": bits, "n_min": n, "p_nmin": first_primes(n)[-1],
                        "prev_n": last, "t": round(time.time() - t0, 1)}
                out.write(json.dumps(rec2) + "\n"); out.flush(); print(rec2, flush=True)
                break
            last = n


if __name__ == "__main__":
    main()
