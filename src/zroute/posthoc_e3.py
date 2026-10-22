"""POST HOC (labelled; not pre-registered): why E3's orbit hit counts are over-dispersed.
Hypothesis: consecutive targets lambda, 2*lambda share split structure (a split u/v of lambda gives the representation
2u/v of 2 lambda, often still the Gauss-reduced one), so hits come in runs.  Measures P(hit_{j+1} | hit_j) versus
the unconditional hit rate, and the fraction of hit pairs whose split is (2u, v) or (u, v/2) of the previous one."""
import json, sys
from src.common.instances import balanced_semiprime
from src.detGRH.core import FactorBase, fb_bound
from src.zroute import core

def main(count=50, J=2000, bits=40):
    tot = hits = hh = derived = 0
    for seed in range(count):
        N = balanced_semiprime(bits, seed).N
        fb = FactorBase.make(N, fb_bound(N))
        lam = pow(2, N.bit_length() + 1, N)
        prev = None
        for j in range(J):
            (u, v), _ = core.gauss_reduce(N, lam)
            h = v != 0 and fb.vector(u) is not None and fb.vector(v) is not None
            tot += 1
            if h:
                hits += 1
            if prev is not None and prev[0]:
                if h:
                    hh += 1
                    pu, pv = prev[1]
                    if {(abs(u), abs(v))} & {(abs(2 * pu), abs(pv)), (abs(pu), abs(pv) // 2)}:
                        derived += 1
            prev = (h, (u, v))
            lam = lam * 2 % N
    print(json.dumps({"rate": hits / tot, "P_hit_given_prev_hit": hh / max(1, hits), "consecutive_hit_pairs": hh,
                      "derived_from_prev_split": derived}))

if __name__ == "__main__":
    main()
