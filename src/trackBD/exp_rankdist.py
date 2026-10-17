"""Experiment B3b: EXACT rank distributions of multiplication-by-sketch on
F_l[X]/(X^r - c) over all admissible weights c, for both hidden primes l in {p, q}.

HARNESS script (uses p, q).  For uniform admissible c mod N (CRT: independent
uniform residues mod p and mod q), the rank-gap probability is exactly
    P_gap = 1 - sum_j pi_{p,j} pi_{q,j},
pi_{l,j} = fraction of admissible c mod l with rank_l = j.
Admissible: Res(X^r - c, W) != 0 mod l, W = X (X+1) (X^2+X+1)  (public filter: the
rejection test gcd(Res(X^r-c, W), N) is computable from N alone).

Criterion decomposition (Theorem B1 / Cor. B2): rank_l(c) = r - deg gcd(X^r - c, F_m) over
F_l, m = N/l.  Each deficient c is attributed to
   'Fl'  : c = x^r for some x in F_l with F_m(x) = 0, x non-universal (x not in {0,-1}, x^2+x+1 != 0)
   'ext' : only roots in proper extensions of F_l.
Output one JSONL row per (N, r).

Usage: .venv/bin/python -m src.trackBD.exp_rankdist --out experiments/raw/trackBD_rankdist.jsonl
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from multiprocessing import Pool

import flint

from src.common import instances as I

W_COEFFS = [0, 1, 2, 2, 1]  # X (X+1) (X^2+X+1) = X^4 + 2X^3 + 2X^2 + X


def side(l: int, m: int, r: int):
    """Rank distribution over admissible c in F_l for f = X^r - c."""
    one = flint.nmod_poly([1], l)
    X = flint.nmod_poly([0, 1], l)
    Wp = flint.nmod_poly(W_COEFFS, l)
    # F_l-roots of F_m, non-universal
    roots = [x for x in range(l) if (pow(1 + x, m, l) - 1 - pow(x, m, l)) % l == 0
             and x not in (0, l - 1) and (x * x + x + 1) % l]
    root_c = {pow(x, r, l) for x in roots}
    dist = Counter()
    attrib = Counter()
    adm = 0
    for c in range(l):
        f = flint.nmod_poly([(-c) % l] + [0] * (r - 1) + [1], l)
        if f.resultant(Wp) == 0:
            continue
        adm += 1
        F = ((X + one).pow_mod(m, f) - one - X.pow_mod(m, f)) % f
        dg = f.degree() if F.is_zero() else f.gcd(F).degree()
        rank = r - dg
        dist[rank] += 1
        if rank < r:
            attrib["Fl" if c in root_c else "ext"] += 1
    return {"adm": adm, "dist": {str(k): v for k, v in sorted(dist.items())},
            "attrib": dict(attrib), "n_Fl_roots": len(roots)}


def run(args):
    p, q, r, tag = args
    N = p * q
    sp = side(p, q, r)
    sq = side(q, p, r)
    gap = 1.0
    num = 0
    for j in range(r + 1):
        num += sp["dist"].get(str(j), 0) * sq["dist"].get(str(j), 0)
    P_gap_num = sp["adm"] * sq["adm"] - num
    return {"N": N, "p": p, "q": q, "r": r, "tag": tag, "side_p": sp, "side_q": sq,
            "P_gap_num": P_gap_num, "P_gap_den": sp["adm"] * sq["adm"],
            "P_gap": P_gap_num / (sp["adm"] * sq["adm"])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rs", type=int, nargs="+", default=[2, 3, 4, 5, 6, 7, 8, 12, 16])
    ap.add_argument("--inst", type=int, default=40)
    ap.add_argument("--bits", type=int, nargs="+", default=[16, 20, 24])
    ap.add_argument("--procs", type=int, default=12)
    ap.add_argument("--out", default="experiments/raw/trackBD_rankdist.jsonl")
    a = ap.parse_args()
    pairs = [(101, 113, "codex"), (107, 193, "codex"), (1009, 1481, "codex"), (199, 241, "codex47959")]
    seen = set()
    for b in a.bits:
        for s in range(a.inst):
            inst = I.balanced_semiprime(b, s)
            p, q = inst.labels["factors"]
            if (p, q) not in seen:
                seen.add((p, q))
                pairs.append((p, q, f"balanced{b}"))
    jobs = [(p, q, r, t) for p, q, t in pairs for r in a.rs]
    with Pool(a.procs) as pool, open(a.out, "w") as fh:
        for i, rec in enumerate(pool.imap_unordered(run, jobs)):
            fh.write(json.dumps(rec) + "\n")
            if i % 100 == 0:
                print(f"{i}/{len(jobs)}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
