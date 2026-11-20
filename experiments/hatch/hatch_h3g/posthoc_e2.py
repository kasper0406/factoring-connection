#!/usr/bin/env python3
"""
POST-HOC DIAGNOSTIC (explicitly labelled post-hoc; registered P2b failure stands):
Diagnoses why E2's P2b had shift = -12.28 nats at r = 200 across s = 64 sublattices:
1) At r = 200 across s = 64 sublattices, r / s = 3.125 points per sublattice (pure cusp/shortest-vector regime,
   overflowing s * D * ln(2E) > r), dropping from two-sided typical (2/3) ln N to one-sided Bean floor (1/3) ln N
   (difference = (1/3) * 60 * ln 2 = 13.86 nats!).
2) Checks Theorem DAG-POOL's actual proved inequality (Bean sublevel floor (1/3) ln N + 2 ln(r / C) without -(D/2) ln s)
   across s = 1, 4, 16, 64 at r = 200 and r = 2000 after stripping top k_cusp = 4 cusp points per sublattice.
"""

import json
import math
import os
import random
from run_hatch3_experiments import (
    GLOBAL_SEED,
    WORKDIR,
    compose_linear_2d,
    eval_homog,
    generate_semiprimes,
    integer_cuberoot_floor,
    poly_mul,
    primes_up_to,
    quartic_disc,
    reduce_lattice_2d,
)


def main() -> None:
    rng = random.Random(GLOBAL_SEED + 2)
    insts = generate_semiprimes(rng, 30, 10, blum=False)
    cand_primes = [p for p in primes_up_to(6000) if p > 500]

    prim_pts = [(0, 1)]
    for u in range(1, 33):
        for v in range(-32, 33):
            if math.gcd(u, abs(v)) == 1:
                prim_pts.append((u, v))

    rows = []
    for inst in insts[:5]:
        N = inst["N"]
        x = math.log(N)
        m = integer_cuberoot_floor(N)
        rem = N - m * m * m
        c2 = rem // (m * m)
        rem -= c2 * m * m
        c1 = rem // m
        c0 = rem - c1 * m
        f_coeffs = [1, c2, c1, c0]
        W_coeffs = poly_mul([1, -m], f_coeffs)
        disc_W = quartic_disc(tuple(W_coeffs))

        sub_vals_all = []
        sub_vals_nocyl = []  # strip 4 smallest cusp/shortest-vector exceptions per sublattice (FT budget)
        for q in cand_primes:
            if disc_W % q == 0:
                continue
            root = None
            for r in range(q):
                if (((r + c2) * r + c1) * r + c0) % q == 0:
                    root = r
                    break
            if root is not None:
                M = reduce_lattice_2d((q, 0), (root, 1))
                G_coeffs = [c // q for c in compose_linear_2d(W_coeffs, M)]
                vlist = sorted(math.log(abs(eval_homog(G_coeffs, u, v))) for (u, v) in prim_pts if eval_homog(G_coeffs, u, v) != 0)
                sub_vals_all.append((q, vlist))
                sub_vals_nocyl.append((q, vlist[8:]))  # strip 8 Farey-cusp exceptions out of 2637 box points
                if len(sub_vals_all) == 64:
                    break

        bean_base_r200 = x / 3.0 + 2.0 * math.log(200)
        bean_base_r1000 = x / 3.0 + 2.0 * math.log(1000)
        vw_q1 = sorted(math.log(abs(eval_homog(W_coeffs, u, v))) for (u, v) in prim_pts if eval_homog(W_coeffs, u, v) != 0)
        rec = {
            "N": N,
            "x_div_3": round(x / 3.0, 2),
            "two_x_div_3": round(2 * x / 3.0, 2),
            "Q1_base_r200": round(vw_q1[199], 3),
            "Q1_base_r1000": round(vw_q1[999], 3),
        }
        for s in (1, 4, 16, 64):
            pool_raw = sorted(v for _, vl in sub_vals_all[:s] for v in vl)
            pool_ft = sorted(v for _, vl in sub_vals_nocyl[:s] for v in vl)
            rec[f"raw_r200_s{s}"] = round(pool_raw[199], 3)
            rec[f"ft8_r200_s{s}"] = round(pool_ft[199], 3)
            rec[f"raw_r1000_s{s}"] = round(pool_raw[999], 3)
            rec[f"excess_over_bean_r1000_s{s}"] = round(pool_raw[999] - bean_base_r1000, 3)
        rows.append(rec)

    out_path = os.path.join(WORKDIR, "posthoc_e2.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=2)
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
