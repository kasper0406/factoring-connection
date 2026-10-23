"""POST HOC analysis of E3 (not pre-registered): joint vs product goodness within height classes of the fiber family.

The fiber fields F + G*L0*H have heights growing with max|h_i|; height heterogeneity induces positive association
between all targets (Jensen), which the pre-registered null did not remove.  Here J and Pi are computed within each
class max|h_i| = const (the H list is recovered from the deterministic enumeration order used for d = 3) and summed.
Run: .venv/bin/python -m src.multidescent.posthoc <e3d3.jsonl>
"""
from __future__ import annotations

import json
import sys

from src.multidescent.experiments import _all_H, CONFIGS, fiber_from_h


def classwise(cols, classes):
    J = Pi = 0.0
    for c in sorted(set(classes)):
        idx = [i for i, k in enumerate(classes) if k == c]
        if not idx:
            continue
        good = [[1 if col[i] > 0 else 0 for i in idx] for col in cols]
        n = len(idx)
        j = sum(all(g[t] for g in good) for t in range(n)) / n
        p = 1.0
        for g in good:
            p *= sum(g) / n
        J += n * j
        Pi += n * p
    return J, Pi


def main():
    rows = [json.loads(l) for l in open(sys.argv[1])]
    cfg = CONFIGS["E3d3"]
    Hs = _all_H(3, cfg["Hh"])
    out = {}
    for key in ("m<=2", "m>=3"):
        out[key] = {"alg": [0.0, 0.0], "null": [0.0, 0.0], "parents": 0}
    for r in rows:
        cs = list(r["counts"].values())
        ns = list(r["null_counts"].values())
        if not cs or not cs[0] or len(ns) != len(cs) or any(len(c) != len(cs[0]) for c in ns):
            continue
        # reconstruct which H were kept (fields that passed the checks, in enumeration order) is not stored;
        # use the stored length: when all 81 passed, classes follow the enumeration order
        if len(cs[0]) != len(Hs):
            continue
        classes = [max(abs(x) for x in h) for h in Hs]
        k = "m<=2" if len(cs) <= 2 else "m>=3"
        a = classwise(cs, classes)
        b = classwise(ns, classes)
        out[k]["alg"][0] += a[0]; out[k]["alg"][1] += a[1]
        out[k]["null"][0] += b[0]; out[k]["null"][1] += b[1]
        out[k]["parents"] += 1
    for k, v in out.items():
        ra = v["alg"][0] / v["alg"][1] if v["alg"][1] else None
        rn = v["null"][0] / v["null"][1] if v["null"][1] else None
        print(k, "parents", v["parents"], "alg J/Pi within height classes", ra and round(ra, 3),
              "null", rn and round(rn, 3), "alg/null", ra and rn and round(ra / rn, 3))


if __name__ == "__main__":
    main()
