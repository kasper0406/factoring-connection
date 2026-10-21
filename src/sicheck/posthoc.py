"""POST-HOC analysis of E1 (disclosed; NOT pre-registered). HARNESS.

Motivation: the pre-registered P1/P4 failed because almost all smooth residues are inherited multiples
(u = g u', v = g v' with g = gcd(u, k) > 1), i.e. duplicates of a relation from k/g. The relevant H_S question is then
about *primitive* outputs (gcd(u, v) = 1), for which v is coprime to every prime of u (those primes do not divide kN).

PH1: O/E on primitive outputs with the naive null p_null(v).
PH2: O/E on primitive outputs with the coprimality-corrected null p_cop(v; u) = density of y-smooth integers among the
     integers coprime to rad(u) in the same window [v - D, v + D] (exact counts).
PH3: distinct fresh relations per interval (the rate that actually feeds linear algebra).
Usage: .venv/bin/python -m src.sicheck.posthoc
"""
from __future__ import annotations

import json
import math
import os
from collections import defaultdict
from math import gcd

import numpy as np

from src.sicheck import core
from src.sicheck.diag import smooth_numbers_upto

RAW = "experiments/raw"


class CoprimeNull:
    def __init__(self, hmax, primes):
        self.primes = primes
        self.vals = smooth_numbers_upto(int(hmax * 1.06) + 400, primes)
        self.cache = {}

    def p_cop(self, v: int, rad_primes) -> float:
        D = max(200, v // 20)
        a, b = max(1, v - D), v + D
        i0 = np.searchsorted(self.vals, a, side="left")
        i1 = np.searchsorted(self.vals, b, side="right")
        win = self.vals[i0:i1]
        R = math.prod(rad_primes)
        if win.size:
            sm = int(np.count_nonzero(np.gcd(win, R) == 1)) if R < 2 ** 62 else sum(1 for w in win.tolist() if gcd(w, R) == 1)
        else:
            sm = 0
        # integers in [a, b] coprime to R: inclusion-exclusion over rad_primes (few primes)
        tot = 0
        ps = list(rad_primes)
        for mask in range(1 << len(ps)):
            d = 1
            bits = 0
            for j, p in enumerate(ps):
                if mask >> j & 1:
                    d *= p
                    bits += 1
            tot += (-1) ** bits * (b // d - (a - 1) // d)
        return sm / tot if tot else 0.0


def main():
    rows = [json.loads(l) for l in open(os.path.join(RAW, "sicheck_E1.jsonl"))]
    out = []
    by = defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0])  # O, Enaive, Vn, Ecop, Vc, nprim, fresh, intervals
    nulls = {}
    for r in rows:
        primes = core.primes_upto(r["y"])
        key = (r["y"], r["h"])
        if key not in nulls:
            nulls[key] = CoprimeNull(r["h"], primes)
        cn = nulls[key]
        N = r["N"]
        acc = by[(r["cell"], r["constructor"])]
        acc[7] += r["K"]
        for k, v, vs, pn, g, dup in zip(r["k"], r["v"], r["vs"], r["pn"], r["g"], r["dup"]):
            if vs and not dup:
                acc[6] += 1
            if g != 1:
                continue
            u = k * N + v
            rp = [p for p in primes if u % p == 0]
            pc = cn.p_cop(v, rp)
            acc[0] += vs; acc[1] += pn; acc[2] += pn * (1 - pn); acc[3] += pc; acc[4] += pc * (1 - pc); acc[5] += 1
    allc = defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0])
    for (cell, c), a in sorted(by.items()):
        for i in range(8):
            allc[c][i] += a[i]
        print(_line(cell, c, a))
        out.append({"cell": cell, "constructor": c, "O_prim": a[0], "E_naive_prim": a[1], "E_cop_prim": a[3],
                    "n_prim": a[5], "fresh": a[6], "intervals": a[7]})
    for c, a in allc.items():
        print(_line("ALL", c, a))
    # PH3: fresh relations versus the naive H_S expectation over all outputs
    print("PH3: fresh distinct relations / naive expected smooth outputs (sum of p_null over all outputs)")
    naive = defaultdict(float)
    for r in rows:
        naive[(r["cell"], r["constructor"])] += sum(r["pn"])
    for (cell, c), a in sorted(by.items()):
        print(f"  {cell:4s} {c:8s} fresh {a[6]:5d}  naive E {naive[(cell, c)]:8.1f}  efficiency {a[6] / naive[(cell, c)]:.3f}")
        for o in out:
            if o["cell"] == cell and o["constructor"] == c:
                o["E_naive_all"] = naive[(cell, c)]
    with open(os.path.join(RAW, "sicheck_E1_posthoc.jsonl"), "w") as f:
        for o in out:
            f.write(json.dumps(o) + "\n")


def ph2b():
    """PH2b (after codex sicheck_2 #4-5): every output whose reduced pair (u/g, v/g) has not occurred before among
    ALL earlier outputs of the same constructor (smooth or not) -- an outcome-independent selection that contains
    every fresh relation. Nulls: p_cop over rad(u/g) (as PH2) and, as a sensitivity check, over rad(u/g * k/g)."""
    rows = [json.loads(l) for l in open(os.path.join(RAW, "sicheck_E1.jsonl"))]
    nulls = {}
    acc = defaultdict(lambda: [0, 0.0, 0.0, 0, 0])  # O, E_u, E_uk, n, fresh_from_nonprimitive
    for r in rows:
        primes = core.primes_upto(r["y"])
        key = (r["y"], r["h"])
        if key not in nulls:
            nulls[key] = CoprimeNull(r["h"], primes)
        cn = nulls[key]
        N = r["N"]
        seen = set()
        a = acc[(r["cell"], r["constructor"])]
        for k, v, vs, g in zip(r["k"], r["v"], r["vs"], r["g"]):
            u = k * N + v
            red = (u // g, v // g)
            if red in seen:
                continue
            seen.add(red)
            up, vp, kp = u // g, v // g, k // g
            ru = [p for p in primes if up % p == 0]
            ruk = [p for p in primes if up % p == 0 or kp % p == 0]
            a[0] += vs; a[1] += cn.p_cop(vp, ru); a[2] += cn.p_cop(vp, ruk); a[3] += 1
            a[4] += int(vs and g > 1)
    tot = defaultdict(lambda: [0, 0.0, 0.0, 0, 0])
    print("PH2b: first-seen reduced pairs (all outputs); O/E vs coprime null over rad(u/g) and over rad(u/g*k/g)")
    for (cell, c), a in sorted(acc.items()):
        for i in range(5):
            tot[c][i] += a[i]
        print(f"  {cell:4s} {c:8s} n={a[3]:6d} O={a[0]:5d} E_u={a[1]:8.1f} O/E_u={a[0] / a[1]:.3f} "
              f"E_uk={a[2]:8.1f} O/E_uk={a[0] / a[2]:.3f} fresh-from-g>1={a[4]}")
    for c, a in tot.items():
        print(f"  ALL  {c:8s} n={a[3]:6d} O={a[0]:5d} O/E_u={a[0] / a[1]:.3f} O/E_uk={a[0] / a[2]:.3f} "
              f"fresh-from-g>1={a[4]}")


def _line(cell, c, a):
    O, En, Vn, Ec, Vc, npr, fresh, K = a
    zn = (O - En) / math.sqrt(Vn) if Vn else float("nan")
    zc = (O - Ec) / math.sqrt(Vc) if Vc else float("nan")
    return (f"{cell:4s} {c:8s} primitive outputs {npr:6d}  O={O:5d}  naive E={En:7.1f} (O/E {O / En if En else 0:.3f}, "
            f"z {zn:+.1f})  coprime-null E={Ec:7.1f} (O/E {O / Ec if Ec else 0:.3f}, z {zc:+.1f})  "
            f"fresh relations {fresh} / {K} intervals")


if __name__ == "__main__":
    import sys
    if sys.argv[1:] == ["ph2b"]:
        ph2b()
    else:
        main()
