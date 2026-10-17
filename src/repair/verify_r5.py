"""Checks of two computational claims in the GPT-6.1-sol report review (codex/r5_report_review.out).  DIAGNOSTIC.

Exact integer arithmetic only.  Everything here knows p and q.

Item 3 (bit tower).  For N = pq = 1 (mod 8), p = eps + 4u, q = eps + 4v, b = [2 is a non-square mod N]:
    uv mod 4 = b * (1 + 2*[N = 1 (mod 16)])  (mod 4),
hence s = p + q mod 64 is already a function of (N mod 64, eps, QR_2), with no class-group rank needed:
    s = eps(N+1) - 16*eps*(uv mod 4)   (mod 64).
Also a diagnostic for the next level, s mod 128 (needs uv mod 8): is it a function of
(N mod 128, eps, QR_2), and if not, of (N mod 128, eps, QR_2, r4 and r8 of Cl(-4N) and Cl(-8N))?

Item 10 (C-2R).  p = 257, m = 8, E = 70, v = 186: phi(x) = ((x + m)^E + v)^2 is an exact fold of Z_17 with outer
degree d = 7 < 9 = ceil(L/2), although (p - 1) does not divide k(E - 1) for any k <= 12.

Run:  .venv/bin/python -m src.repair.verify_r5 --bound 1000 --out experiments/raw/verify_r5.json
"""
from __future__ import annotations

import argparse
import json
import time
from collections import defaultdict
from math import gcd
from typing import Dict, List, Tuple

from src.repair import verify_r4 as V


# ---------------------------------------------------------------------------
# item 3
# ---------------------------------------------------------------------------
def two_power_torsion(D: int) -> Dict[str, int]:
    """r2, r4, r8 of Cl(D) by brute force: |G[2^j]| for j = 1, 2, 3."""
    one = V.identity_form(D)
    n2 = n4 = n8 = 0
    for f in V.reduced_forms(D):
        f2 = V.compose(f, f)
        f4 = V.compose(f2, f2)
        f8 = V.compose(f4, f4)
        n2 += f2 == one
        n4 += f4 == one
        n8 += f8 == one
    e2, e4, e8 = (x.bit_length() - 1 for x in (n2, n4, n8))
    assert (n2, n4, n8) == (1 << e2, 1 << e4, 1 << e8)
    return {"r2": e2, "r4": e4 - e2, "r8": e8 - e4}


def item3(bound: int = 1000, class_groups: bool = True) -> Dict[str, object]:
    P = V.odd_primes_below(bound)
    recs = []
    for i, p in enumerate(P):
        for q in P[i + 1:]:
            N = p * q
            if N % 8 != 1:
                continue
            eps = 1 if p % 4 == 1 else -1
            u, v = (p - eps) // 4, (q - eps) // 4
            b = 0 if (V.legendre(2, p) == 1 and V.legendre(2, q) == 1) else 1
            rec = {"p": p, "q": q, "N": N, "eps": eps, "b": b, "s": p + q,
                   "uv_claim": (u * v) % 4 == (b * (1 + 2 * (N % 16 == 1))) % 4,
                   "s64_formula": (p + q) % 64 == (eps * (N + 1) - 16 * eps * ((u * v) % 4)) % 64}
            if class_groups:
                g4, g8 = two_power_torsion(-4 * N), two_power_torsion(-8 * N)
                rec["cg"] = (g4["r4"], g8["r4"], g4["r8"], g8["r8"])
            recs.append(rec)

    def ambiguity(keyf, mod):
        groups = defaultdict(set)
        for r in recs:
            groups[keyf(r)].add(r["s"] % mod)
        amb = [k for k, vals in groups.items() if len(vals) > 1]
        return {"keys": len(groups), "ambiguous_keys": len(amb),
                "pairs_in_ambiguous_keys": sum(1 for r in recs if keyf(r) in set(amb))}

    out: Dict[str, object] = {
        "bound": bound,
        "pairs_N_1_mod_8": len(recs),
        "uv_claim_fail": sum(not r["uv_claim"] for r in recs),
        "s64_formula_fail": sum(not r["s64_formula"] for r in recs),
        "s64_given_N64_eps_QR2": ambiguity(lambda r: (r["N"] % 64, r["eps"], r["b"]), 64),
        "s128_given_N128_eps_QR2": ambiguity(lambda r: (r["N"] % 128, r["eps"], r["b"]), 128),
    }
    if class_groups:
        out["s128_given_N128_eps_QR2_r4r8"] = ambiguity(
            lambda r: (r["N"] % 128, r["eps"], r["b"]) + r["cg"], 128)
    return out


# ---------------------------------------------------------------------------
# item 10
# ---------------------------------------------------------------------------
def fold(p: int, m: int, E: int, v: int) -> Dict[str, object]:
    """Exact brute force over all of F_p for phi(x) = ((x+m)^E + v)^2 and the interval Z_L, L = 2m + 1."""
    L = 2 * m + 1
    phi = [pow((pow((x + m) % p, E, p) + v) % p, 2, p) for x in range(p)]
    ZL = {(-i) % p for i in range(L)}
    image = {phi[x] for x in ZL}
    zero_set = {x for x in range(p) if phi[x] in image}
    twist = [k for k in range(1, 13) if (k * (E - 1)) % (p - 1) == 0]
    return {"p": p, "m": m, "L": L, "E": E, "v": v, "d": len(image), "exact": zero_set == ZL,
            "one_layer_degree": m + 1, "three_L_over_8": 3 * L / 8, "passes_3L8": 8 * len(image) <= 3 * L,
            "gcd_E_pm1": gcd(E, p - 1), "character_exponent": (2 * E - 2) % (p - 1) == 0,
            "twist_orders_le_12": twist}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bound", type=int, default=1000)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    t0 = time.time()
    res = {"item3": item3(args.bound), "item10": fold(257, 8, 70, 186)}
    res["seconds"] = round(time.time() - t0, 1)
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(res, fh, sort_keys=True, indent=1)
            fh.write("\n")
    print(json.dumps(res, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
