"""Experiment B2b: degree-r sketches on seeded balanced semiprimes, all readouts,
classified, with per-instance predictions from Theorem B1.

HARNESS script.  Operational per instance (N only): for r in R_SET draw K public
random c in [2, N-2] (rng "trackBD_bal|N|seed|r"), S = F_N mod (X^r - c); plus
Phi_m (m = 3..60) and 20 fixed random f.  Readouts: coefficient gcds, content
gcd, Euclid/resultant (sketch.all_readouts).

Diagnostic: classify each hit (diagnose.component_vanishing) and predict the
Euclid/rank-split rate for X^r - c with uniform c from the F_l-roots of F_m:
   pred_r = |{x^r : x in A_p}|/p + |{y^r : y in A_q}|/q,
   A_l = {x in F_l : F_m(x) = 0, x not in {0,-1}, x^2+x+1 != 0}
(extension-field roots ignored: they contribute O(r/l^2)).
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from collections import Counter
from math import gcd
from multiprocessing import Pool

import flint

from src.common import instances as I
from src.trackBD import diagnose as D
from src.trackBD.exp_tiny import random_family
from src.trackBD.sketch import all_readouts, binomial_modulus, sketch

R_SET = (1, 2, 3, 4, 6, 8, 12, 16, 32, 64, 128, 200)
ROOTS_MAX = 1 << 20


def nonuniv_roots(l: int, m: int):
    return [x for x in range(2, l - 1) if (pow(1 + x, m, l) - 1 - pow(x, m, l)) % l == 0
            and (x * x + x + 1) % l]


_VT = {}


def support_zero(N, p, r, j):
    """Coefficient j of F_N mod (p, X^r - c) is zero because no Lucas-support exponent of
    F_N mod p = F_m(X^p) folds to j (Lucas sparsity; Lehman-type q/p structure)."""
    from src.trackBD.exp_tiny import lucas_valid_t
    key = (N, p)
    if key not in _VT:
        _VT.clear()
        _VT[key] = lucas_valid_t(p, N // p) * p
    import numpy as np
    return not np.any((_VT[key] % r) == j)


def classify(N, factors, f, S, rd):
    """Return (coef_label, content_label, euclid_label) for one sketch's readouts."""
    info = D.component_vanishing(N, factors, f, S)
    labs = {}
    r = len(f) - 1
    binom = all(x == 0 for x in f[1:-1])

    def plab(p):
        inf = info[p]
        if inf["td"]:
            return "TD:" + "+".join(sorted(inf["td"]))
        if inf["full"]:
            ls = {c["label"] for c in inf["comps"] if c["root_v"]} - {"TRIVIAL"}
            return "FULL_" + ("+".join(sorted(ls)) if ls else "TRIVIAL")
        return None

    ps = sorted(info)
    if rd["coef"]:
        best = None
        for j, d in rd["coef"]:
            for p in ps:
                if d % p == 0:
                    l = plab(p) or ("SUPPORT" if binom and support_zero(N, p, r, j) else "ACC")
                    best = l if best is None or (best == "ACC" and l != "ACC") else best
        labs["coef"] = best
    ct = rd["content"]
    if 1 < ct < N:
        labs["content"] = "/".join(sorted({plab(p) or "?" for p in ps if ct % p == 0}))
    eu = rd["euclid"]
    if eu.get("factor"):
        degs = {p: info[p]["deg_gcd"] for p in ps}
        tds = sorted({t for p in ps for t in info[p]["td"]})
        if tds:
            labs["euclid"] = "TD:" + "+".join(tds)
        elif len(set(degs.values())) > 1:
            lo = min(degs.values())
            ex = set()
            for p in ps:
                if degs[p] > lo:
                    ex |= {c["label"] for c in info[p]["comps"] if c["root_v"]}
            ex.discard("TRIVIAL")
            labs["euclid"] = "SPLIT_" + ("+".join(sorted(ex)) if ex else "TRIVIAL")
        else:
            sparse = binom and any(support_zero(N, p, r, j) for p in ps for j in range(r))
            labs["euclid"] = "ACC_sparse" if sparse else "ACC"
    return labs


def run(args):
    bits, seed, K = args
    inst = I.balanced_semiprime(bits, seed)
    N = inst.N
    p, q = inst.labels["factors"]
    fam = []
    for r in R_SET:
        rng = random.Random(f"trackBD_bal|{N}|{seed}|{r}")
        for _ in range(K):
            c = rng.randint(2, N - 2)
            fam.append((f"binom{r}", binomial_modulus(r, c)))
    for m in range(3, 61):
        fam.append(("cyclo", [int(x) for x in flint.fmpz_poly.cyclotomic(m).coeffs()]))
    for i, f in enumerate(random_family()):
        fam.append(("rand", f))
    cnt = Counter()
    for name, f in fam:
        S = sketch(N, f)
        rd = all_readouts(N, f, S)
        cnt[f"{name}|sketches"] += 1
        if rd["coef"] or (1 < rd["content"] < N) or rd["euclid"].get("factor"):
            for k, v in classify(N, inst.factors, f, S, rd).items():
                cnt[f"{name}|{k}|{v}"] += 1
    rec = {"bits": bits, "seed": seed, "N": N, "p": p, "q": q, "K": K, "counts": dict(cnt),
           "null_per_coef": 1 - ((p - 1) * (q - 1) + 1) / N}
    if q <= ROOTS_MAX:
        Ap, Aq = nonuniv_roots(p, q), nonuniv_roots(q, p)
        rec["A_p"], rec["A_q"] = len(Ap), len(Aq)
        rec["pred"] = {str(r): len({pow(x, r, p) for x in Ap}) / p + len({pow(y, r, q) for y in Aq}) / q
                       for r in R_SET}
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bits", type=int, nargs="+", default=[16, 20, 24, 28, 32])
    ap.add_argument("--inst", type=int, default=100)
    ap.add_argument("--K", type=int, default=100)
    ap.add_argument("--procs", type=int, default=10)
    ap.add_argument("--out", default="experiments/raw/trackBD_sketch_bal.jsonl")
    a = ap.parse_args()
    jobs = [(b, s, a.K) for b in a.bits for s in range(a.inst)]
    with Pool(a.procs) as pool, open(a.out, "a") as fh:
        for i, rec in enumerate(pool.imap_unordered(run, jobs)):
            fh.write(json.dumps(rec) + "\n")
            if i % 20 == 0:
                print(f"{i}/{len(jobs)}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
