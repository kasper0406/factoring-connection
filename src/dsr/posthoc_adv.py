"""POST HOC (labelled; not pre-registered): the box obstruction (Prop BOX in notes/dsr_results.md).

(i)  For every E2 instance, D_box's worst-case coverage m* (S_can) versus Psi(N/l, B) + 1, where l is the least
     prime > B: the target a = l is missed by every multiplier w < N/l, since l*w < N is its own least residue.
(ii) Adversarial DSR targets: scan admissible 24-bit balanced semiprimes for primes q <= B whose target q^H has a
     non-smooth least residue r <= N / 2^12; for each found, compare D_box's first hit on that target with
     Psi(N/r, B) (a lower bound on the first hit index) and with the random-walk control.
  PYTHONPATH=. .venv/bin/python -m src.dsr.posthoc_adv > experiments/raw/dsr_posthoc_adv.jsonl
"""
from __future__ import annotations

import json
import math
import random

from src.common.instances import balanced_semiprime
from src.dsr import core
from src.dsr.harness import A_CAN, r_walk


def _next_prime(x: int) -> int:
    y = x + 1
    while any(y % d == 0 for d in range(2, math.isqrt(y) + 1)):
        y += 1
    return y


def part_i() -> None:
    for line in open("experiments/raw/dsr_e2.jsonl"):
        r = json.loads(line)
        if r["list"] != "D_box" or r["set"] != "S_can":
            continue
        N, B = r["N"], r["B"]
        _, P = core.factor_base(N, A_CAN)
        ell = _next_prime(B)
        psi = int(core.smooth_mask(N // ell + 1, P)[: (N - 1) // ell + 1].sum())  # smooth w with w*ell < N
        print(json.dumps({"part": "i", "bits": r["bits"], "seed": r["seed"], "ell": ell, "Psi_N_over_ell": psi,
                          "m_star": r["m_star"], "cap": r["cap"],
                          "bound_holds": (r["m_star"] is None and psi + 1 > r["cap"]) or (
                              r["m_star"] is not None and r["m_star"] >= psi + 1),
                          "gap": None if r["m_star"] is None else r["m_star"] - psi - 1}))


def part_ii(n_scan: int = 4000, thresh_bits: int = 12) -> None:
    found = 0
    for seed in range(10_000, 10_000 + n_scan):
        N = balanced_semiprime(24, seed).N
        try:
            B, P = core.factor_base(N, A_CAN)
        except ValueError:
            continue
        if min(f for f in balanced_semiprime(24, seed).factors) <= 3 * math.log(N) ** 2:
            continue  # preprocessing bound (harness check only)
        H = 8 * math.ceil(math.log2(N))
        for q in P:
            r = pow(q, H, N)
            if r > N >> thresh_bits:
                continue
            m = r
            for p in P:
                while m % p == 0:
                    m //= p
            if m == 1:
                continue  # r itself smooth: a trivial hit, not adversarial
            found += 1
            if found > 12:
                continue
            S = core.s_can(N, P)
            lb = int(core.smooth_mask((N - 1) // r + 1, P).sum())  # number of smooth w <= (N-1)//r
            cap = 200_000
            fb = core.first_hits(N, S, [r], core.d_box(P), cap)[r]
            t = core.walk_length(N, B)
            fw = core.first_hits(N, S, [r], r_walk(P, t, random.Random(f"dsr|adv|{seed}|{q}")), cap)[r]
            print(json.dumps({"part": "ii", "seed": seed, "N": N, "q": q, "B": B, "r": r, "largest_cofactor": m,
                              "lower_bound_Psi(N/r)": lb, "D_box_first_hit": fb, "R_walk_first_hit": fw,
                              "delta": round(float(S.sum()) / N, 4)}), flush=True)
    print(json.dumps({"part": "ii_summary", "scanned": n_scan, "pairs_found": found, "thresh_bits": thresh_bits}))


if __name__ == "__main__":
    part_i()
    part_ii()
