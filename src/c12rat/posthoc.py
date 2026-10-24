"""POST HOC (not pre-registered): E2's sanity ratio with the correct expectation Psi_D(I)/phi(D) (smooth numbers in
I coprime to D, divided by phi(D)); the registered P2b used (Z/phi(D))*density and omitted the coprimality factor.
Usage: python -m src.c12rat.posthoc raw_e2.jsonl out.json"""
import json
import sys

import numpy as np

from src.pc3 import toy as pc3toy

rows = [json.loads(l) for l in open(sys.argv[1])]
X, y = 50_000_000, 50
smooth, _ = pc3toy.build_tables(X, y, y + 1)
rat = []
for r in rows:
    D, X0, Z = r["D"], r["X0"], r["Z"]
    n = np.arange(X0 + 1, X0 + Z + 1)
    sm = smooth[X0 + 1:X0 + Z + 1] & (np.gcd(n, D) == 1)
    phiD = D * np.prod([1 - 1 / p for p in pc3toy.factor_small(D)])
    e = sm.sum() / phiD
    if e > 0:
        rat.append(r["lhs"] / e)
rat = np.array(rat)
out = dict(n=len(rat), median=float(np.median(rat)), q10=float(np.quantile(rat, .1)), q90=float(np.quantile(rat, .9)))
json.dump(out, open(sys.argv[2], "w"), indent=1)
print(out)
