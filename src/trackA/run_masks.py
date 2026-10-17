"""HARNESS: Track A (e) mask cost measurements.

For balanced semiprimes, apply prefix / interval / class / bit-pattern masks
at increasing parameter and record rank, q_mul_dimsum and wall time; check
the prefix test against the factorial-unit predicate (L3) using the factors.
Usage: .venv/bin/python -m src.trackA.run_masks --out experiments/raw/trackA_masks.jsonl
"""
from __future__ import annotations

import argparse
import time
from math import gcd

from ..common import instances as I
from ..common.metrics import Metrics, append_jsonl
from . import masks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="experiments/raw/trackA_masks.jsonl")
    a = ap.parse_args()
    for bits in (24, 32, 40, 48):
        inst = I.balanced_semiprime(bits, 0)
        N = inst.N
        p = inst.factors[0]
        for e in range(4, 21, 2):
            R = 1 << e
            for kind in ("prefix", "class", "bit"):
                m = Metrics()
                t0 = time.perf_counter()
                if kind == "prefix":
                    v = masks.prefix(N, R, m)
                    nz = any(v)
                    truth = R > p  # L3: prefix of length R is zero iff R <= p
                    extra = {"prefix_nonzero": nz, "harness_R_gt_p": truth, "agree": nz == truth}
                elif kind == "class":
                    val = masks.residue_class(N, 3 % R, R, 5, m)
                    extra = {"value_gcd": gcd(val, N)}
                else:
                    val = masks.bit_pattern(N, e - 1, 5, m)
                    extra = {"value_gcd": gcd(val, N)}
                rec = m.to_record(experiment="trackA_masks", N=N, bits=bits, kind=kind,
                                  rank=R, log2_rank=e, log2_p=p.bit_length(),
                                  wall_s=time.perf_counter() - t0, **extra)
                append_jsonl(a.out, rec)
                print(bits, kind, e, rec["counts"].get("q_mul_dimsum"), round(rec["wall_s"], 4),
                      extra, flush=True)


if __name__ == "__main__":
    main()
