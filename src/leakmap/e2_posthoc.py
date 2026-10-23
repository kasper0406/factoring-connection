"""E2 post-hoc (NOT pre-registered): the archimedean CF hit (A1, cf64_n_sq_even).

For each instance recompute the first 64 CF steps with convergents A_i mod N; at every even i with Q_i = r^2 test
whether gcd(A_{i-1} - r, N) is proper (a congruence-of-squares split from the square form; Lehman/SQUFOF
mechanism).  Then re-test the A1 correlation of cf64_n_sq_even on the instances where NO such square splits N
("non-split" subpopulation), with the same strata.
"""
import gzip
import json
import sys
from math import gcd

import gmpy2
import numpy as np

from .cfgenus import strat_z, cluster_z, labels


def squares_split(N, T=64):
    a0 = int(gmpy2.isqrt(N))
    P, Q, a = 0, 1, a0
    Am2, Am1 = 1, a0 % N           # A_{-1}, A_0
    nsq = nsplit = 0
    for i in range(1, T + 1):
        P = a * Q - P
        Q = (N - P * P) // Q
        # here Am1 = A_{i-1}
        if i % 2 == 0 and gmpy2.is_square(Q):
            r = int(gmpy2.isqrt(Q))
            nsq += 1
            g = gcd((Am1 - r) % N, N)
            if 1 < g < N:
                nsplit += 1
        a = (a0 + P) // Q
        Am2, Am1 = Am1, (a * Am1 + Am2) % N
    return nsq, nsplit


def main(path):
    recs = [json.loads(l) for l in gzip.open(path, "rt")]
    rows = []
    for r in recs:
        nsq, nsp = squares_split(r["N"])
        assert nsq == r["cf64_n_sq_even"], (nsq, r["cf64_n_sq_even"])
        rows.append((r, nsq, nsp, labels(r["p"], r["q"])["A1"]))
    y = np.array([t[3] for t in rows], float)
    x = np.array([t[1] for t in rows], float)
    Ns = [t[0]["N"] for t in rows]
    s = np.array([(Nv % 256) * 64 + (Nv >> (Nv.bit_length() - 6)) % 64 for Nv in Ns])
    cp = np.array([t[0]["p"] for t in rows]); cq = np.array([t[0]["q"] for t in rows])
    anysq = x > 0
    anysplit = np.array([t[2] > 0 for t in rows])
    out = dict(file=path, count=len(rows), frac_with_square=float(anysq.mean()),
               frac_square_events_splitting=float(sum(t[2] for t in rows) / max(1, sum(t[1] for t in rows))),
               frac_instances_split=float(anysplit.mean()),
               P_A1_given_split=float(y[anysplit].mean()) if anysplit.any() else None,
               P_A1_given_square_nosplit=float(y[anysq & ~anysplit].mean()) if (anysq & ~anysplit).any() else None,
               P_A1_given_nosquare=float(y[~anysq].mean()),
               z_all=strat_z(x, y, s), zcr_all=cluster_z(x, y, s, cp, cq))
    m = ~anysplit
    out["z_nonsplit_subpop"] = strat_z(x[m], y[m], s[m])
    out["zcr_nonsplit_subpop"] = cluster_z(x[m], y[m], s[m], cp[m], cq[m])
    # q/p of the splitting instances
    ratios = [t[0]["q"] / t[0]["p"] for t in rows if t[2] > 0]
    out["median_q_over_p_splitting"] = float(np.median(ratios)) if ratios else None
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    res = [main(p) for p in sys.argv[1:]]
    json.dump(res, open("experiments/analysis/leakmap_e2_posthoc_A1.json", "w"), indent=1)
