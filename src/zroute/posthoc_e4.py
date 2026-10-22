"""POST HOC (labelled; not pre-registered): why E4's tries per prime are heavy-tailed (P4d failed).
Hypothesis: for a target a with Gauss-reduced (u0, v0), the targets a*w (small w) have representations derived from
(u0, v0) -- (w u0, v0), or (u0, v0/w) when w | v0 -- so the split outcomes of a*w are inherited from a's lattice and
are strongly correlated across w.  For one 40-bit instance: tries per prime, and for primes needing > 50 tries, the
fraction of failed w whose reduced vector is (w u0, v0)/g or (u0, v0/w) up to sign."""
import json
from math import gcd
from src.common.instances import balanced_semiprime
from src.detGRH.core import FactorBase, fb_bound
from src.zroute import core

def main(seed=3, bits=40, budget=20000):
    N = balanced_semiprime(bits, seed).N
    fb = FactorBase.make(N, fb_bound(N))
    ws = [w for w in range(1, budget + 1) if fb.vector(w) is not None]
    K = max(sum(fb.vector(w)[1:]) for w in ws)
    H = K + (2 * N).bit_length() + 2
    tries, derived, checked = {}, 0, 0
    for q in fb.primes:
        a = pow(q, H, N)
        (u0, v0), _ = core.gauss_reduce(N, a)
        for t, w in enumerate(ws):
            (u, v), _ = core.gauss_reduce(N, a * w % N)
            ok = v != 0 and fb.vector(u) is not None and fb.vector(v) is not None
            if ok:
                tries[q] = t + 1
                break
            if t < 200:
                checked += 1
                g = gcd(w * u0, v0)
                cands = {(abs(w * u0) // g, abs(v0) // g)}
                if v0 % w == 0:
                    cands.add((abs(u0), abs(v0) // w))
                if (abs(u), abs(v)) in cands:
                    derived += 1
    tv = sorted(tries.values())
    print(json.dumps({"N": N, "primes": len(fb.primes), "complete": len(tv), "median": tv[len(tv)//2], "max": tv[-1],
                      "n_over_50": sum(1 for x in tv if x > 50), "mean": sum(tv)/len(tv),
                      "failed_w_checked": checked, "failed_w_derived_from_a": derived}))

if __name__ == "__main__":
    main()
