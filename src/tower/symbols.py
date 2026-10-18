"""Part 3 checks: residue symbols, Redei-symbol identities, and the quartic-reciprocity question.  DIAGNOSTIC.

S1  classical identity, p = 1 (mod 8):  (-4/p)_8 = (2/p)_4 * psi(p),  psi(p) = (-1)^((p-1)/8).
    (Frobenius of sqrt(1+i) = 2^(1/4) zeta_16 over Q(zeta_8); this is [-1,2,p] = [-1,2,p] in Stevenhagen's notation.)
S2  single-prime prototypes of the Redei computation in notes/tower_results.md section 3 (Stevenhagen Thm 4.1 +
    linearity Thm 7.9 + Def 7.1), checked against the brute-force class groups of src.tower.classgroup:
      p = 1 (mod 8):  r8 Cl(-4p) = [(-4/p)_8 = 1]           (Stevenhagen eq. (28), Barrucand-Cohn)
                      r8 Cl(-8p) = [(2/p)_4 = 1]            (R8 = [p,-8,2] = [p,-1,2] + [p,2,2] = (2/p)_4)
      p = 7 (mod 8):  r8 Cl(-8p) = [p = 15 (mod 16)]       (R8 = [-p,8,2] = [-p,2,2] = psi_2(-p))
    and, for p = 1 (mod 8), [T in Cl(-4p)^4] XOR [T' in Cl(-8p)^4] = [p = 9 (mod 16)] (the prototype of T3/T4).
S3  quartic reciprocity question: for pairs with p = q = 1 (mod 8), (p/q) = 1, p = q (mod 16), is the hidden bit
    h = [psi(p) = -1] related to (p/q)_4 (q/p)_4, (2/p)_4 (2/q)_4, (-4/p)_8 (-4/q)_8 ?  Contingency tables.

Run:  .venv/bin/python -m src.tower.symbols --bound 2000 --out experiments/raw/tower128_symbols.json
"""
from __future__ import annotations

import argparse
import json
import time
from collections import Counter, defaultdict
from typing import Dict

from src.repair.verify_r4 import legendre, odd_primes_below
from src.tower.classgroup import two_sylow, rank
from src.tower.run_tower import ambiguous_prime_form, depth_info


def power_symbol(a: int, p: int, k: int) -> int:
    """(a/p)_k = a^((p-1)/k) mod p in {+1, -1}, for k | p - 1 and a a (k/2)-th power residue mod p.  Asserts +-1."""
    assert (p - 1) % k == 0
    t = pow(a % p, (p - 1) // k, p)
    assert t in (1, p - 1), (a, p, k, t)
    return 1 if t == 1 else -1


def psi_pm(p: int) -> int:
    """psi(p) for p = +-1 (mod 8), as +-1: +1 iff p = +-1 (mod 16)."""
    return 1 if p % 16 in (1, 15) else -1


def s2_single(p: int) -> Dict[str, object]:
    out: Dict[str, object] = {"p": p}
    if p % 8 == 1:
        S4 = two_sylow(-4 * p, return_elems=True)
        S8 = two_sylow(-8 * p, return_elems=True)
        out["r8_-4p"], out["r8_-8p"] = rank(S4, 3), rank(S8, 3)
        out["ok_-4p"] = (rank(S4, 3) == 1) == (power_symbol(-4, p, 8) == 1)
        out["ok_-8p"] = (rank(S8, 3) == 1) == (power_symbol(2, p, 4) == 1)
        z4 = depth_info(-4 * p, ambiguous_prime_form(-4 * p, 2), S4)["depth"] >= 2
        z8 = depth_info(-8 * p, ambiguous_prime_form(-8 * p, 2), S8)["depth"] >= 2
        out["ok_depth_xor"] = (z4 ^ z8) == (p % 16 == 9)
    elif p % 8 == 7:
        S8 = two_sylow(-8 * p)
        out["r8_-8p"] = rank(S8, 3)
        out["ok_-8p_7mod8"] = (rank(S8, 3) == 1) == (p % 16 == 15)
    return out


def run(bound: int = 2000, single_bound: int = 20000) -> Dict[str, object]:
    t0 = time.time()
    res: Dict[str, object] = {"bound": bound, "single_bound": single_bound}
    P1 = [p for p in odd_primes_below(single_bound) if p % 8 == 1]
    res["S1_checked"] = len(P1)
    res["S1_fail"] = [p for p in P1 if power_symbol(-4, p, 8) != power_symbol(2, p, 4) * psi_pm(p)]
    singles = [s2_single(p) for p in odd_primes_below(single_bound) if p % 8 in (1, 7)]
    s2 = {}
    for key in ("ok_-4p", "ok_-8p", "ok_depth_xor", "ok_-8p_7mod8"):
        vals = [r[key] for r in singles if key in r]
        s2[key] = {"checked": len(vals), "fail": sum(not v for v in vals),
                   "first_fail": [r["p"] for r in singles if key in r and not r[key]][:5]}
    res["S2"] = s2
    # S3
    tabs: Dict[str, Counter] = defaultdict(Counter)
    P = odd_primes_below(bound)
    n = 0
    for i, p in enumerate(P):
        if p % 8 != 1:
            continue
        for q in P[i + 1:]:
            if q % 8 != 1 or p % 16 != q % 16 or legendre(p, q) != 1:
                continue
            n += 1
            h = int(p % 16 == 9)
            cell = (p * q) % 128
            tabs["(p/q)_4(q/p)_4"][(cell, power_symbol(p, q, 4) * power_symbol(q, p, 4), h)] += 1
            tabs["(2/p)_4(2/q)_4"][(cell, power_symbol(2, p, 4) * power_symbol(2, q, 4), h)] += 1
            tabs["(-4/p)_8(-4/q)_8"][(cell, power_symbol(-4, p, 8) * power_symbol(-4, q, 8), h)] += 1
            tabs["(p/q)_4"][(cell, power_symbol(p, q, 4), h)] += 1
    from math import log2

    def cmi(tab: Counter) -> float:
        tot = sum(tab.values())
        c_, cx, ch = Counter(), Counter(), Counter()
        for (c, x, h), k in tab.items():
            c_[c] += k
            cx[(c, x)] += k
            ch[(c, h)] += k
        return sum(k * log2(k * c_[c] / (cx[(c, x)] * ch[(c, h)])) for (c, x, h), k in tab.items()) / tot

    def deterministic(tab: Counter) -> bool:
        seen = defaultdict(set)
        for (c, x, h) in tab:
            seen[(c, x)].add(h)
        return all(len(v) == 1 for v in seen.values())

    import random

    def cmi_null(tab: Counter, nperm: int = 200) -> float:
        rows = [(c, x, h) for (c, x, h), k in tab.items() for _ in range(k)]
        rng = random.Random("tower128-S3")
        by = defaultdict(list)
        for c, x, h in rows:
            by[c].append((x, h))
        vals = []
        for _ in range(nperm):
            t = Counter()
            for c, lst in by.items():
                xs = [x for x, _ in lst]
                rng.shuffle(xs)
                for x, (_, h) in zip(xs, lst):
                    t[(c, x, h)] += 1
            vals.append(cmi(t))
        vals.sort()
        return {"mean": round(sum(vals) / len(vals), 5), "p95": round(vals[int(0.95 * (len(vals) - 1))], 5)}

    res["S3_pairs"] = n
    res["S3"] = {k: {"cmi_bits_given_N_mod_128": round(cmi(t), 5), "h_function_of_it": deterministic(t),
                     "null": cmi_null(t),
                     "table_(x,h)": {f"{x},{h}": sum(v for (c, xx, hh), v in t.items() if xx == x and hh == h)
                                     for x in (1, -1) for h in (0, 1)}}
                 for k, t in tabs.items()}
    res["seconds"] = round(time.time() - t0, 1)
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bound", type=int, default=2000)
    ap.add_argument("--single_bound", type=int, default=20000)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    res = run(args.bound, args.single_bound)
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(res, fh, sort_keys=True, indent=1)
            fh.write("\n")
    print(json.dumps(res, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
