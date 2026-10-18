"""Analysis for E1 (pre-registered in notes/nfsfruit_prereg.md).
Usage: .venv/bin/python -m src.nfsfruit.analyze_e1 experiments/raw/nfsfruit_e1.jsonl"""
import json
import sys
from collections import defaultdict
from math import sqrt

SEMI = ("blum", "mixed", "both1")


def p_zero(k, w1):
    """P(random linear map F2^k -> F2^{w1} has rank < w1)."""
    if w1 == 1:
        return 2.0 ** (-k)
    full = 1.0
    for i in range(w1):
        full *= 1 - 2.0 ** (i - k)
    return 1 - full


def z(obs, n, p):
    if n == 0 or p in (0, 1):
        return float("nan")
    return (obs - n * p) / sqrt(n * p * (1 - p))


def main():
    recs = [json.loads(l) for l in open(sys.argv[1])]
    out = []
    pr = out.append
    st = defaultdict(int)
    for r in recs:
        st[(r["family"], r.get("status"))] += 1
    pr("## P0 status counts")
    for k in sorted(st):
        pr(f"  {k[0]:10s} {k[1]:24s} {st[k]}")
    ok = [r for r in recs if r.get("status") == "ok"]
    pr(f"  lin_bad total = {sum(r.get('lin_bad', 0) for r in ok)}, lin_ok total = {sum(r.get('lin_ok', 0) for r in ok)}")
    pr(f"  nonsquare_basis total = {sum(r.get('nonsquare_basis', 0) for r in recs if 'nonsquare_basis' in r)}")
    pr(f"  errors: {[ (r['family'], r['cls'], r['bits'], r['seed'], r.get('error')) for r in recs if r.get('status') == 'error']}")
    pr("\n## P1/P2 rank Phi = w-1 on instances with dim K0 >= 10")
    fams = sorted({r["family"] for r in ok})
    for f in fams:
        rs = [r for r in ok if r["family"] == f]
        elig = [r for r in rs if r["dim_K0"] >= 10]
        bad = [r for r in elig if r["phi_rank"] < (2 if r["cls"] == "tri" else 1)]
        fr = sum(r["fruitful_basis"] for r in rs if r["cls"] in SEMI)
        tot = sum(r["dim_K0"] for r in rs if r["cls"] in SEMI)
        dims = sorted(r["dim_K0"] for r in rs)
        pr(f"  {f:10s} ok={len(rs):3d} eligible={len(elig):3d} rank-deficient={len(bad):3d} "
           f"dimK0 median={dims[len(dims)//2] if dims else None} min={dims[0] if dims else None} "
           f"basis fruitful (semiprimes) {fr}/{tot} = {fr / tot if tot else float('nan'):.3f}")
        for r in bad:
            pr(f"     DEFICIENT: {r['cls']} {r['bits']} seed={r['seed']} dimK0={r['dim_K0']} rank={r['phi_rank']}")
    pr("\n## P3 subset statistics: P(rank Phi|sub < w-1) by k = dim K0(sub)")
    for f in fams:
        for grp, classes, w1 in (("semi", SEMI, 1), ("tri", ("tri",), 2)):
            line = f"  {f:10s} {grp:4s}"
            for k in (1, 2, 3, 4):
                n = obs = 0
                for r in ok:
                    if r["family"] == f and r["cls"] in classes:
                        for d, rk in r.get("subsets", []):
                            if d == k:
                                n += 1
                                obs += int(rk < w1)
                p = p_zero(k, w1)
                line += f" | k={k}: {obs}/{n} exp {p:.3f} z={z(obs, n, p):+.1f}"
            pr(line)
    pr("\n## P4 per class (semiprimes), pooled k=1..4 standardized: sum(obs-exp)/sqrt(sum var)")
    for f in fams:
        line = f"  {f:10s}"
        for c in SEMI:
            E = V = O = 0.0
            for r in ok:
                if r["family"] == f and r["cls"] == c:
                    for d, rk in r.get("subsets", []):
                        if 1 <= d <= 4:
                            p = p_zero(d, 1)
                            E += p
                            V += p * (1 - p)
                            O += int(rk < 1)
            line += f" | {c}: obs {int(O)} exp {E:.1f} z={(O - E) / sqrt(V) if V else float('nan'):+.1f}"
        pr(line)
    pr("\n## P5 NFS: dim ker(parity only) - dim K0 (binding character constraints)")
    for f in fams:
        xs = [r["dim_ker_parity_only"] - r["dim_K0"] for r in ok if r["family"] == f and "dim_ker_parity_only" in r]
        if xs:
            xs.sort()
            pr(f"  {f:10s} n={len(xs)} min={xs[0]} median={xs[len(xs)//2]} max={xs[-1]}")
    print("\n".join(out))


if __name__ == "__main__":
    main()
