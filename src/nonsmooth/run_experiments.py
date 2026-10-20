"""HARNESS: experiments E1 (operational CM/anomalous factorer) and E2 (coverage census).

Pre-registered in notes/nonsmooth_prereg.md. Usage:
  .venv/bin/python -m src.nonsmooth.run_experiments e1
  .venv/bin/python -m src.nonsmooth.run_experiments e2
Raw JSONL: experiments/raw/nonsmooth_e{1,2}.jsonl. Summary: experiments/analysis/nonsmooth_summary.json.
"""
from __future__ import annotations

import json
import os
import random
import sys
import time

import gmpy2

from src.nonsmooth import census, cheng

RAW = "experiments/raw"
ANA = "experiments/analysis"


def special_prime(D, bits, rng):
    """Prime p with 4p - 1 = D s^2, s odd, p of about `bits` bits (harness)."""
    sb = max(2, (bits - D.bit_length() + 2) // 2)
    while True:
        s = rng.getrandbits(sb) | 1 | (1 << (sb - 1))
        if (D * s * s + 1) % 4:
            continue
        p = (D * s * s + 1) // 4
        if gmpy2.is_prime(p):
            return int(p)


def random_prime(rng, bits):
    while True:
        x = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if gmpy2.is_prime(x):
            return int(x)


def e1(seeds=10, sizes=(128, 256, 512), controls=60):
    rows = []
    for n in sizes:
        for D in cheng.D_LIST:
            for seed in range(seeds):
                rng = random.Random(f"nonsmooth-e1|special|{n}|{D}|{seed}")
                p = special_prime(D, n // 2, rng)
                q = random_prime(rng, n // 2)
                N = p * q
                t = time.time()
                f, DD, c = cheng.cheng_factor(N, seed=seed)       # sees N only
                dt = time.time() - t
                cls = None
                if f is not None:
                    # STRUCTURAL iff the factor found is the special prime and the hit is a Z-gcd
                    cls = "STRUCTURAL" if (f == p and DD == D and c.log and c.log[-1][0] == "Z") \
                        else "ACCIDENTAL"
                rows.append({"exp": "E1", "arm": "special", "n": n, "D": D, "seed": seed,
                             "success": f in (p, q), "class": cls, "D_found": DD,
                             "tries": c.tries, "mul": c.mul, "gcd": c.gcd, "sec": dt})
        for seed in range(controls):
            rng = random.Random(f"nonsmooth-e1|control|{n}|{seed}")
            p = random_prime(rng, n // 2)
            q = random_prime(rng, n // 2)
            N = p * q
            t = time.time()
            f, DD, c = cheng.cheng_factor(N, seed=seed)
            dt = time.time() - t
            cls = None
            if f is not None:
                sp = any((4 * r - 1) % D == 0 and gmpy2.is_square((4 * r - 1) // D)
                         for r in (p, q) for D in cheng.D_LIST)
                cls = "STRUCTURAL" if sp else "ACCIDENTAL"
            rows.append({"exp": "E1", "arm": "control", "n": n, "seed": seed,
                         "success": f in (p, q), "class": cls, "tries": c.tries,
                         "mul": c.mul, "gcd": c.gcd, "sec": dt})
    return rows


def summarize_e1(rows):
    out = {}
    for arm in ("special", "control"):
        for n in sorted({r["n"] for r in rows}):
            rs = [r for r in rows if r["arm"] == arm and r["n"] == n]
            out[f"{arm}|{n}"] = {
                "count": len(rs),
                "success": sum(r["success"] for r in rs),
                "structural": sum(r["class"] == "STRUCTURAL" for r in rs),
                "accidental": sum(r["class"] == "ACCIDENTAL" for r in rs),
                "mean_tries": sum(r["tries"] for r in rs) / len(rs),
                "max_tries": max(r["tries"] for r in rs),
                "mean_mul": sum(r["mul"] for r in rs) / len(rs),
                "mean_sec": sum(r["sec"] for r in rs) / len(rs),
            }
    return out


def e2(count=20000, sizes=(24, 32, 40, 48)):
    rows = []
    for b in sizes:
        rows += census.sample(b, count, 0)
    return rows


def _write(name, rows):
    os.makedirs(RAW, exist_ok=True)
    with open(f"{RAW}/{name}.jsonl", "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")


def _update_summary(key, val):
    os.makedirs(ANA, exist_ok=True)
    path = f"{ANA}/nonsmooth_summary.json"
    s = {}
    if os.path.exists(path):
        with open(path) as fh:
            s = json.load(fh)
    s[key] = val
    with open(path, "w") as fh:
        json.dump(s, fh, indent=1, sort_keys=True)


if __name__ == "__main__":  # pragma: no cover
    which = sys.argv[1]
    if which == "e1":
        rows = e1()
        _write("nonsmooth_e1", rows)
        _update_summary("E1", summarize_e1(rows))
        print(json.dumps(summarize_e1(rows), indent=1))
    elif which == "e2":
        rows = e2()
        _write("nonsmooth_e2", rows)
        s = census.summarize(rows)
        slopes = {str(th): census.fit_slope(s, th) for th in (0.5, 0.75, 0.9)}
        # prefactor c_theta = frac / p^{-(1-theta)/2}, using p ~ 2^(b - 0.5)
        pref = {}
        for th in (0.5, 0.75, 0.9):
            pref[str(th)] = {b: (v["frac_le_p_theta"][str(th)]["frac"]
                                 * 2 ** ((int(b) - 0.5) * (1 - th) / 2))
                             for b, v in s.items()}
        _update_summary("E2", {"by_bits": s, "slopes": slopes, "prefactor": pref})
        print(json.dumps({"by_bits": s, "slopes": slopes, "prefactor": pref}, indent=1))
