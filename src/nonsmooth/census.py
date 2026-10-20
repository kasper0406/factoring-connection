"""DIAGNOSTIC: coverage of the anomalous-CM paradigm over random primes.

For a prime p the CM/anomalous method (cheng.py; Cheng 2002 in general form) needs the
discriminant -D, D = squarefree part of 4p - 1 (the trace-1 CM order), and costs at least
deg H_D = h(-D) ring operations over Z/N. This module samples primes p (it generates them, so it
may factor 4p - 1) and records D. Exact integers only; logs are for the summary fit.
"""
from __future__ import annotations

import json
import math
import random

import gmpy2
from sympy import factorint


def squarefree_part(m):
    f = factorint(m)
    d = 1
    for r, e in f.items():
        if e % 2:
            d *= r
    return d


def random_prime(rng, bits):
    while True:
        x = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if gmpy2.is_prime(x):
            return x


def sample(bits, count, seed):
    rng = random.Random(f"nonsmooth-census|{bits}|{seed}")
    rows = []
    for _ in range(count):
        p = random_prime(rng, bits)
        D = squarefree_part(4 * p - 1)
        rows.append({"bits": bits, "p": p, "D": D, "D_mod8": D % 8})
    return rows


def summarize(rows, thetas=(0.25, 0.5, 0.75, 0.9)):
    out = {}
    by = {}
    for r in rows:
        by.setdefault(r["bits"], []).append(r)
    for b, rs in sorted(by.items()):
        ratios = sorted(math.log(r["D"]) / math.log(r["p"]) for r in rs)
        n = len(rs)
        fr = {}
        for th in thetas:
            # exact integer comparison D <= p^theta  <=>  D^den <= p^num  (theta rational)
            num, den = {0.25: (1, 4), 0.5: (1, 2), 0.75: (3, 4), 0.9: (9, 10)}[th]
            k = sum(1 for r in rs if r["D"] ** den <= r["p"] ** num)
            fr[str(th)] = {"count": k, "frac": k / n}
        out[str(b)] = {
            "n": n,
            "median_logD_over_logp": ratios[n // 2],
            "q10_logD_over_logp": ratios[n // 10],
            "all_D_3mod8": all(r["D_mod8"] == 3 for r in rs),
            "min_D": min(r["D"] for r in rs),
            "frac_le_p_theta": fr,
        }
    return out


def fit_slope(summary, theta):
    """Least-squares slope of log2(frac) against bits (float; fit only). Skips zero counts."""
    pts = []
    for b, s in summary.items():
        c = s["frac_le_p_theta"][str(theta)]
        if c["count"] > 0:
            pts.append((int(b), math.log2(c["frac"])))
    if len(pts) < 2:
        return None
    mx = sum(x for x, _ in pts) / len(pts)
    my = sum(y for _, y in pts) / len(pts)
    sxx = sum((x - mx) ** 2 for x, _ in pts)
    sxy = sum((x - mx) * (y - my) for x, y in pts)
    return sxy / sxx


if __name__ == "__main__":  # pragma: no cover
    import sys
    rows = []
    for b in (24, 32, 40, 48):
        rows += sample(b, int(sys.argv[1]) if len(sys.argv) > 1 else 200, 0)
    s = summarize(rows)
    print(json.dumps(s, indent=1))
    for th in (0.5, 0.75, 0.9):
        print(th, fit_slope(s, th))
