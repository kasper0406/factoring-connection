"""Part 2: extended class-group key for every pair p < q < bound with N = pq = 1 (mod 8).  DIAGNOSTIC.

For each pair and each public discriminant D in DISCS (all negative, all = 0,1 mod 4, primitive forms):
  * the exact 2-Sylow type of Cl(D) (hence r2, r4, r8, r16, ...), via src.tower.classgroup.two_sylow;
  * for each public prime l | D/N (l in {2, 3, 5, 7}) whose ambiguous form (l, b, c) is primitive: its 2-depth,
    the largest j with [form] in Cl(D)^(2^j) (capped at e = v2(h)), and its order (1 or 2).
Pairs with p <= 7 are excluded so that gcd(k, N) = 1 for every multiplier k.  Discriminants -N, -2N, -5N, -6N are
not = 0,1 (mod 4) when N = 1 (mod 8) and are skipped (as are their -4kN counterparts' duplicates).

Run:  .venv/bin/python -m src.tower.run_tower --bound 2000 --workers 20 --out experiments/raw/tower128_invariants.jsonl
      gzip -9 experiments/raw/tower128_invariants.jsonl
Config: experiments/configs/tower128.json
"""
from __future__ import annotations

import argparse
import json
import time
from math import gcd
from multiprocessing import Pool
from typing import Dict, List, Optional, Tuple

from src.repair.verify_r4 import compose, identity_form, odd_primes_below, reduce_form
from src.tower.classgroup import two_sylow
from src.tower.formula import hidden_bit, omega_exp, psi_exp

DISCS: Dict[str, Tuple[int, Tuple[int, ...]]] = {
    # name: (multiplier k with D = -k N, public primes dividing k)
    "-4N": (4, (2,)), "-8N": (8, (2,)), "-12N": (12, (2, 3)), "-20N": (20, (2, 5)),
    "-24N": (24, (2, 3)), "-28N": (28, (2, 7)), "-3N": (3, (3,)), "-7N": (7, (7,)),
}


def ambiguous_prime_form(D: int, l: int) -> Optional[Tuple[int, int, int]]:
    """The form (l, b, c) of discriminant D with l | b (ramified prime l | D), if primitive; else None."""
    for b in range(0, 2 * l):
        if (b - D) % 2 == 0 and (b * b - D) % (4 * l) == 0 and b % l == 0:
            c = (b * b - D) // (4 * l)
            if gcd(gcd(l, b), c) == 1:
                return reduce_form((l, b, c))
    return None


def depth_info(D: int, x: Tuple[int, int, int], syl: Dict[str, object]) -> Dict[str, int]:
    """2-depth of a class x of order <= 2: max j <= e with x in Cl^(2^j).  Uses the 2-Sylow S (x lies in S)."""
    one = identity_form(D)
    if x == one:
        return {"order": 1, "depth": syl["e"]}
    S = syl["elems"]
    j, layer = 0, S
    while True:
        nxt = {compose(y, y) for y in layer}
        if x not in nxt:
            return {"order": 2, "depth": j}
        j += 1
        layer = nxt


def invariants(pq: Tuple[int, int]) -> Dict[str, object]:
    p, q = pq
    N = p * q
    eps = 1 if p % 4 == 1 else -1
    rec: Dict[str, object] = {
        "p": p, "q": q, "N": N, "eps": eps, "QR2": int(pow(2, (p - 1) // 2, p) == 1),  # (2/p) = (2/q) = +1
        "N_mod_128": N % 128, "h_bit": hidden_bit(p, q), "s_mod_128": (p + q) % 128,
        "psi_p": psi_exp(p), "psi_q": psi_exp(q), "omega_p": omega_exp(p), "omega_q": omega_exp(q),
    }
    cg = {}
    for name, (k, ls) in DISCS.items():
        D = -k * N
        syl = two_sylow(D, seed=0, return_elems=True)
        d = {"h": syl["h"], "ranks": syl["ranks"], "type": syl["type"]}
        for l in ls:
            x = ambiguous_prime_form(D, l)
            if x is not None:
                di = depth_info(D, x, syl)
                d[f"amb{l}"] = di
        cg[name] = d
    rec["cg"] = cg
    return rec


def pairs(bound: int, pmin: int = 11) -> List[Tuple[int, int]]:
    P = [p for p in odd_primes_below(bound) if p >= pmin]
    return [(p, q) for i, p in enumerate(P) for q in P[i + 1:] if (p * q) % 8 == 1]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bound", type=int, default=2000)
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    t0 = time.time()
    todo = pairs(args.bound)
    todo.sort(key=lambda t: -t[0] * t[1])  # big first for load balance
    with Pool(args.workers) as pool, open(args.out, "w") as fh:
        for i, rec in enumerate(pool.imap_unordered(invariants, todo, chunksize=8)):
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
            if (i + 1) % 1000 == 0:
                print(f"{i + 1}/{len(todo)}  {time.time() - t0:.0f}s", flush=True)
    print(json.dumps({"bound": args.bound, "pairs": len(todo), "seconds": round(time.time() - t0, 1)}))


if __name__ == "__main__":
    main()
