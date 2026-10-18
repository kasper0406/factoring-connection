"""Analysis for Addendum A (E1b, E2, E3).
Usage: .venv/bin/python -m src.nfsfruit.analyze_add experiments/raw/nfsfruit_add.jsonl"""
import json
import sys
from collections import defaultdict
from math import sqrt

from .analyze_e1 import p_zero

SEMI = ("blum", "mixed", "both1")


def cluster_z(groups, p):
    """groups: list of per-instance lists of 0/1; cluster-robust z for mean = p."""
    allv = [v for g in groups for v in g]
    n = len(allv)
    if n == 0 or p in (0.0, 1.0):
        return n, sum(allv), float("nan")
    obs = sum(allv)
    # variance of the sum using cluster totals
    var = sum((sum(g) - len(g) * p) ** 2 for g in groups)
    if var == 0:
        var = n * p * (1 - p)
    return n, obs, (obs - n * p) / sqrt(var)


def main():
    recs = [json.loads(l) for l in open(sys.argv[1])]
    st = defaultdict(int)
    for r in recs:
        st[(r.get("exp"), r.get("status"))] += 1
    print("## status", dict(st))
    errs = [r for r in recs if r.get("status") == "error"]
    if errs:
        print("errors:", errs[:5])
    # ---------------- E1b
    print("\n## E1b prefix nullity: deficient-rate at k vs H_rand (cluster-robust z)")
    e1b = [r for r in recs if r.get("exp") == "e1b" and r.get("status") == "ok"]
    fams = sorted({r["family"] for r in e1b})
    for f in fams:
        for grp, classes, w1 in (("semi", SEMI, 1), ("tri", ("tri",), 2)):
            line = f"  {f:10s} {grp:4s}"
            for k in range(1, 7):
                groups = []
                for r in e1b:
                    if r["family"] == f and r["cls"] in classes:
                        g = [o[k - 1] for o in r["prefix"] if o[k - 1] is not None]
                        if g:
                            groups.append(g)
                n, obs, zz = cluster_z(groups, p_zero(k, w1))
                line += f" | k={k}: {obs}/{n} (exp {p_zero(k, w1):.3f}) z={zz:+.1f}"
            print(line)
    print("  per-class (semiprimes, k=1..3 pooled, standardized):")
    for f in fams:
        line = f"  {f:10s}"
        for c in SEMI:
            E = V = O = 0.0
            for r in e1b:
                if r["family"] == f and r["cls"] == c:
                    for o in r["prefix"]:
                        for k in (1, 2, 3):
                            if o[k - 1] is not None:
                                p = p_zero(k, 1)
                                E += p
                                V += p * (1 - p)
                                O += o[k - 1]
            line += f" | {c}: {int(O)} vs {E:.1f} z={(O - E) / sqrt(V) if V else float('nan'):+.1f}"
        print(line)
    # ---------------- E2
    print("\n## E2 Jacobi-row identity c(S) = Phi(S) (Blum)")
    e2 = [r for r in recs if r.get("exp") == "e2" and r.get("status") == "ok"]
    for f in sorted({r["family"] for r in e2}):
        rs = [r for r in e2 if r["family"] == f]
        ag = sum(r["agree"] for r in rs)
        di = sum(r["disagree"] for r in rs)
        c0 = [r["sub_c0"] for r in rs]
        c0_dims = sum(d for d, _ in c0)
        c0_bad = sum(1 for d, rk in c0 if rk != 0)
        c1 = [r["sub_c1"] for r in rs]
        print(f"  {f:9s} instances={len(rs)} basis checks agree={ag} disagree={di} | "
              f"c=0 subset: total dim {c0_dims}, instances with rank Phi>0: {c0_bad} | "
              f"c=1 subset: total dim {sum(d for d, _ in c1)}, rank>0 in {sum(1 for d, rk in c1 if rk)}")
    # ---------------- E3
    print("\n## E3 square-multiplier linear sieve")
    e3 = [r for r in recs if r.get("exp") == "e3" and r.get("status") == "ok"]
    print(f"  instances ok={len(e3)}, hidden_ok={sum(r['hidden_ok'] for r in e3)}, hidden_bad={sum(r['hidden_bad'] for r in e3)}")
    for grp, classes, p in (("semi", SEMI, 0.5), ("tri", ("tri",), 0.75)):
        rs = [r for r in e3 if r["cls"] in classes]
        has_odd = sum(1 for r in rs if r.get("odd_basis", 0) > 0)
        for lab in ("odd", "even"):
            v = [r[lab + "_fruit"] for r in rs if r.get(lab + "_fruit") is not None]
            n = len(v)
            o = sum(v)
            zz = (o - n * p) / sqrt(n * p * (1 - p)) if n else float("nan")
            print(f"  {grp:4s} {lab:4s}: fruitful {o}/{n} = {o / n if n else float('nan'):.3f} (exp {p}) z={zz:+.1f}")
        print(f"  {grp:4s} instances with an odd-weight dependency: {has_odd}/{len(rs)}")
    for c in SEMI + ("tri",):
        rs = [r for r in e3 if r["cls"] == c]
        v = [r["odd_fruit"] for r in rs if r.get("odd_fruit") is not None]
        print(f"    {c:6s} odd fruitful {sum(v)}/{len(v)}")


if __name__ == "__main__":
    main()
