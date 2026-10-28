"""HARNESS for E2 (cubic NFS base + h rows descended into the same field) and E3 (x-injectivity:
two-term dependencies among NFS relations).  Factors of N are used only for labels (omega).

Usage: python -m src.fruit.harness_e2 <config.json> <out.jsonl>
"""
from __future__ import annotations

import json
import random
import sys
import time
from math import gcd, isqrt

from ..nfsfruit import families as fam
from ..nfsfruit.instances import composite
from .linalg import hybrid_rank, kernel_basis, sample_forced, sample_nonzero
from .nfs_hybrid import build_base, FactorFound
from .qs_hybrid import fruitful


def is_square(n: int) -> bool:
    return n >= 0 and isqrt(n) ** 2 == n


def two_term(ab, m, rows):
    """E3 counts: duplicate full parity rows; pairs with (a-bm)(a'-b'm) a perfect square."""
    seen = {}
    dup_pairs = 0
    in_dup = set()
    for i, r in enumerate(rows):
        if r in seen:
            for j in seen[r]:
                dup_pairs += 1
                in_dup.update((i, j))
            seen[r].append(i)
        else:
            seen[r] = [i]
    # rational squares: group by squarefree kernel of z (exact, z is B-smooth)
    return dup_pairs, len(in_dup)


def rational_square_pairs(rationals, rel_r):
    ker = {}
    pairs = 0
    members = set()
    for i, (z, fr) in enumerate(zip(rationals, rel_r)):
        k = (z < 0, tuple(sorted(p for p, e in fr.items() if e & 1)))
        lst = ker.setdefault(k, [])
        for j in lst:
            pairs += 1
            members.update((i, j))
        lst.append(i)
    return pairs, len(members)


def sqrt_branch_check(s, base_idx, desc_idx):
    """Prop SR check: square roots via every inert prime agree up to sign."""
    O = s.O
    elems = [O.linear(*s.ab[i]) for i in base_idx]
    for j in desc_idx:
        for (_side, _l, e, (a, b, _c, _f, _ar)) in s.desc[j]["kids"]:
            elems += [O.linear(a, b)] * e
    dp = O.deriv_at_alpha()
    g = O.mul(O.mul(dp, dp), O.prod(elems))
    roots = []
    for l in s.inert:
        try:
            v = O.sqrt(g, l)
        except ValueError:
            continue
        if v is not None:
            roots.append(tuple(v))
    if not roots:
        return None
    v0 = roots[0]
    neg = tuple(-x for x in v0)
    return all(v == v0 or v == neg for v in roots), len(roots)


def run_instance(cls, bits, seed, p, cfg):
    N, fs = composite(cls, bits, seed)
    t0 = time.time()
    base = dict(cls=cls, bits=bits, seed=seed, N=str(N), omega=len(fs))
    try:
        s = build_base(N, seed, p["B"], p["A"], p["Bb"], p["nchar"], p["excess"], p["B1"], cfg["J"], cfg["kmax"])
    except FactorFound as e:
        return [dict(base, accidental_factor=int(e.g))]
    recs = []
    rng = random.Random(f"fruit-e2|{N}|{seed}|rounds")
    M = len(s.ab)
    for rd in range(cfg["rounds"]):
        s.desc = []
        try:
            for _ in range(cfg["hmax"]):
                if not s.add_descent(rng, cfg["max_trials"]):
                    raise RuntimeError("descent budget exhausted")
        except FactorFound as e:
            recs.append(dict(base, round=rd, accidental_factor=int(e.g)))
            continue
        s.finalize_hybrid()
        if rd == 0:
            dp, dm = two_term(s.ab, s.m, s.rows)
            rp, rm = rational_square_pairs(s.rationals, s.rel_r)
            recs.append(dict(base, kind="E3", M=M, ncols=s.ncols, full_dup_pairs=dp, full_dup_members=dm,
                             rat_sq_pairs=rp, rat_sq_members=rm))
            if "rand" in p:
                q = p["rand"]
                try:
                    sr = fam.nfs3(N, seed, q["B"], q["A"], q["Bb"], q["Hmul"], q["nchar"], "nfs3_rand",
                                  q["excess"], q["Hshift"])
                    dp2, dm2 = two_term(sr.ab, sr.m, sr.rows)
                    rp2, rm2 = rational_square_pairs(sr.rationals, sr.rel_r)
                    recs.append(dict(base, kind="E3rand", M=len(sr.ab), ncols=sr.ncols, full_dup_pairs=dp2,
                                     full_dup_members=dm2, rat_sq_pairs=rp2, rat_sq_members=rm2))
                except fam.FactorFound as e:
                    recs.append(dict(base, kind="E3rand", accidental_factor=int(e.g)))
        kern0 = kernel_basis(s.rows)
        S0 = sample_nonzero(kern0, rng)
        b0 = None
        if S0 is not None:
            out = s.hybrid_congruence([i for i in range(M) if S0 >> i & 1], [])
            if out is not None:
                X, Y = out
                assert (X * X - Y * Y) % N == 0
                b0 = fruitful(X, Y, N)
        recs.append(dict(base, kind="E2", round=rd, h=0, M=M, dimK=len(kern0), base_uniform=b0))
        for h in cfg["hs"]:
            drows = s.drows[:h]
            hr = hybrid_rank(s.rows, drows)
            kern = kernel_basis(s.rows + drows)
            rec = dict(base, kind="E2", round=rd, M=M, ncols=s.ncols, dimK=len(kern), **hr)
            for mode in ("forced", "uniform"):
                S = sample_forced(kern, M, rng) if mode == "forced" else sample_nonzero(kern, rng)
                if S is None:
                    rec[mode] = None
                    continue
                bi = [i for i in range(M) if S >> i & 1]
                di = [j for j in range(h) if S >> (M + j) & 1]
                out = s.hybrid_congruence(bi, di)
                if out is None:
                    rec[mode] = "nonsquare"
                    continue
                X, Y = out
                assert (X * X - Y * Y) % N == 0
                rec[mode] = fruitful(X, Y, N)
                if mode == "uniform":
                    rec["uniform_uses_desc"] = bool(di)
                if mode == "forced" and h == cfg["hs"][-1]:
                    rec["sqrt_branch"] = sqrt_branch_check(s, bi, di)
            recs.append(rec)
    for r in recs:
        r["secs_instance"] = round(time.time() - t0, 2)
        r["counters"] = dict(s.counters)
    return recs


def main():
    cfg = json.load(open(sys.argv[1]))
    with open(sys.argv[2], "w") as fo:
        for bits in cfg["bits"]:
            p = cfg[str(bits)]
            for cls in cfg["classes"]:
                for seed in range(cfg["seeds"]):
                    for r in run_instance(cls, bits, seed, p, cfg):
                        fo.write(json.dumps(r) + "\n")
                    fo.flush()


if __name__ == "__main__":
    main()
