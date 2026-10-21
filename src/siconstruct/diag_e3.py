"""Post-hoc diagnostic for E3 (disclosed in the results note, not registered): re-run the E3 cells with the same
seeded lists and count successful cofactor tests WITH multiplicity (one smooth u can be found through several s),
plus the number of distinct smooth u in the covered windows.  Writes experiments/raw/siconstruct_E3diag.jsonl."""
from __future__ import annotations

import json
import os
import random

from src.common.instances import balanced_semiprime

from .core import ht_local_density, is_smooth, lpath_params, primes_upto, random_smooth_list
from .run_exp import RAW


def main():
    K = 16
    f = open(os.path.join(RAW, "siconstruct_E3diag.jsonl"), "w")
    for bits in (40, 48, 56):
        pp = lpath_params(bits)
        y, h = pp["y"], pp["h"]
        P = primes_upto(y)
        lh = h.bit_length() - 1
        for seed in (0, 1):
            N = balanced_semiprime(bits, seed).N
            for S0exp in (lh // 2, lh - 4):
                rng = random.Random(f"E3|{bits}|{seed}|{S0exp}")
                S0 = 1 << S0exp
                S = random_smooth_list(S0, 2 * S0, P, 64, rng)
                succ = 0
                tested = 0
                pred = 0.0
                cache = {}
                distinct = set()
                for s in S:
                    for k in range(1, K + 1):
                        kN = k * N
                        w = (-kN) % s
                        j = 0
                        while w <= h and j < 2048:
                            m = (kN + w) // s
                            bl = m.bit_length()
                            if bl not in cache:
                                cache[bl] = ht_local_density(float(2 ** bl) * 0.75, y, P)
                            pred += cache[bl]
                            tested += 1
                            j += 1
                            if is_smooth(m, P):
                                succ += 1
                                distinct.add(kN + w)
                            w += s
                rec = {"bits": bits, "seed": seed, "S0exp": S0exp, "tested": tested, "successful_tests": succ,
                       "distinct": len(distinct), "pred_tests": pred, "ratio_tests": succ / pred,
                       "multiplicity": succ / max(1, len(distinct))}
                f.write(json.dumps(rec) + "\n")
                f.flush()
                print(rec, flush=True)
    f.close()


if __name__ == "__main__":
    main()
