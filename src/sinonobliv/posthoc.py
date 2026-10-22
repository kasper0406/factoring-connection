"""POST-HOC (disclosed; written after E1's P1 failed in 4 cells): local-residue correction for consecutive pieces.

Model: a y-smooth integer has residue distribution q_p(0) = a_p (measured: fraction of smooth t in the range
divisible by p), q_p(c) = (1 - a_p)/(p - 1) for c != 0, independently across primes. For r consecutive pieces
t, t+1, ..., t+r-1 the joint/independent ratio at p is  C_p = p^{r-1} * sum_{c mod p} prod_j q_p(c + j),
which is 1 when q_p is uniform. Prediction: R ~ prod_{p <= y} C_p. Run:
  .venv/bin/python -m src.sinonobliv.posthoc
"""
from __future__ import annotations

import json
import os

import numpy as np

from .core import primes_upto, smooth_mask_window
from .run_exp import RAW, _t_range


def local_factor(t0: int, T: int, y: int, r: int) -> float:
    m = smooth_mask_window(t0, T, y)
    ts = np.nonzero(m)[0].astype(np.int64) + t0
    n = len(ts)
    prod = 1.0
    for p in primes_upto(y):
        a = float(np.count_nonzero(ts % p == 0)) / n
        q = np.full(p, (1.0 - a) / (p - 1))
        q[0] = a
        s = 0.0
        for c in range(p):
            v = 1.0
            for j in range(r):
                v *= q[(c + j) % p]
            s += v
        prod *= p ** (r - 1) * s
    return prod


def main():
    recs = [json.loads(l) for l in open(os.path.join(RAW, "sinonobliv_E1.jsonl")) if '"family"' in l and "P2" not in l]
    out = open(os.path.join(RAW, "sinonobliv_E1_posthoc.jsonl"), "w")
    for rec in recs:
        if rec["family"] not in ("cp2", "cp3"):
            continue
        r = 2 if rec["family"] == "cp2" else 3
        C = local_factor(rec["t0"], rec["T"], rec["y"], r)
        o = {"family": rec["family"], "bits": rec["bits"], "y": rec["y"], "R": rec["R"], "local_factor": C,
             "R_corrected": rec["R"] / C}
        print(json.dumps(o))
        out.write(json.dumps(o) + "\n")
    out.close()




def e2b_seed_spread(n_seeds: int = 6):
    """POST-HOC (disclosed; after P4/P5 failed only in the C1 Johnson cell): rerun that cell with fresh seeds to
    measure the seed-to-seed spread of the heavy-tailed (d, z) counts. Run: python -m src.sinonobliv.posthoc e2b"""
    import math as _m
    from .core import (DensityCurve, divisors_from_factored, johnson_log_threshold, lcm_upto_factored,
                       value_of_factored)
    from .msample import joint_tuples_indep, predicted_smooth_first, predicted_tuples
    X, h, B, y = 1 << 36, 1 << 10, 32, 1 << 8
    fac = lcm_upto_factored(B)
    P = value_of_factored(fac)
    divs = divisors_from_factored(fac)
    sig = johnson_log_threshold(_m.log(h + 1), _m.log(P))
    curve = DensityCurve(y, 2.0, 37.0, width=1 << 18)
    pt = predicted_tuples(divs, h, 1, sig)
    ps = predicted_smooth_first(divs, h, 1, sig, X, curve)
    out = open(os.path.join(RAW, "sinonobliv_E2B_posthoc.jsonl"), "w")
    for s in range(n_seeds):
        rng = np.random.default_rng(90000 + s)
        res = joint_tuples_indep(X, h, P, B, 1, sig, y, 20000, rng)
        o = {"seed": 90000 + s, **res, "ratio_tuples": res["tuples"] / (20000 * pt),
             "ratio_z1_smooth": res["tuples_z1_smooth"] / (20000 * ps)}
        print(json.dumps(o))
        out.write(json.dumps(o) + "\n")
    out.close()


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "e2b":
        e2b_seed_spread()
    else:
        main()
