"""E1 of `prim` (pre-registered in notes/prim_prereg.md): instrumented PRIM programs.

Part A (full-r, per prime):  uniform-r zero sets Z_v(p) of every affine node; coset structure; fresh-random model;
                             KS check; Theorem U check (partner counts vs |Z_dev|).                  -> P3, P4, P5, P7
Part B (pairs):              direct evaluation mod N, p, q of every node; Jacobi symbols; model probabilities at
                             both primes; batteries, PT-X, non-split gains.                         -> P1, P2, P6, P8

Usage: python -m src.prim.e1 --b 14 [--pairs 200000]
Outputs: experiments/raw/prim_e1A_b{b}.jsonl, experiments/analysis/prim_e1_b{b}.json
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import random
import resource
import time
from collections import Counter, defaultdict
from math import gcd

import gmpy2
import numpy as np

from .prog import (LucasNode, c_p_model, control_nodes, jacobi_values, nested_nodes, summed_nodes)
from .zerosets import (L_cube, coset_union, frozen_set, kelley_bound, ks_bound, ks_threshold, max_vanishing_coset,
                       order_mod, primes_b, zero_set)

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20  # macOS: bytes


# ----------------------------------------------------------------------------- node tables
def node_tables():
    nest = nested_nodes()
    summ = summed_nodes()
    ctrl = control_nodes()
    affine = {}
    for k, (nd, _) in nest.items():
        if nd.affine():
            affine[k] = nd
    for k, (nd, _) in summ.items():
        affine[k] = nd
    for k, nd in ctrl.items():
        affine[k] = nd
    return nest, summ, ctrl, affine


# ----------------------------------------------------------------------------- Part A
def part_a(b: int, out_raw: str):
    nest, summ, ctrl, affine = node_tables()
    primes = primes_b(b)
    parr = np.array(primes, dtype=np.int64)
    half = 1 << (b // 2)
    agg = defaultdict(lambda: defaultdict(float))
    p3_viol = 0
    p4_viol = 0
    p4_checked = 0
    with open(out_raw, "w") as fh:
        for p in primes:
            n = p - 1
            r = np.arange(n, dtype=np.int64)
            L = L_cube(n)
            partners = parr[parr != p]
            rq = partners % n
            o2, o3 = order_mod(2, p), order_mod(3, p)
            rec = {"p": p, "ord2": o2, "ord3": o3, "nodes": {}}
            for name, nd in affine.items():
                Z = nd.vec(p, r) == 0
                hits = Z[rq]
                maj_zero = hits.sum() * 2 > hits.size
                Zdev = ~Z if maj_zero else Z
                dev_cnt = int((~hits).sum() if maj_zero else hits.sum())
                if dev_cnt > 2 * int(Zdev.sum()):
                    p3_viol += 1
                VL = coset_union(Z, L)
                d = {"Z": int(Z.sum()), "dev": dev_cnt, "Zdev": int(Zdev.sum()), "VL": int(VL.sum()),
                     "maj_zero": bool(maj_zero), "Cmax": max_vanishing_coset(Z)}
                if name in nest:
                    bases = nd.bases()
                    flag = any(order_mod(c, p) <= half for c in bases if c not in (0, 1, -1))
                    model = nest[name][1]
                    mv = model.prob_vec(p, r)
                    d["model"] = float(mv.sum())
                    d["flag"] = bool(flag)
                    # POST HOC (not registered): FROZEN class and the residual outside COSET u FROZEN
                    Fr = frozen_set(nd, p)
                    expl = VL | Fr
                    d["Fr_Z"] = int((Z & Fr).sum())
                    d["other_Z"] = int((Z & ~expl).sum())
                    d["other_model"] = float(mv[~expl].sum())
                if name == "C_P":
                    d["model"] = float(c_p_model().prob_vec(p, r).sum())
                if name in summ or name == "C_F":
                    t = summ[name][1] if name in summ else 2
                    Vks = coset_union(Z, ks_threshold(n))
                    res = int((Z & ~Vks).sum())
                    bnd = ks_bound(n, t)
                    d["ks_res"] = res
                    d["ks_bound"] = bnd
                    d["kelley"] = kelley_bound(n, t, d["Cmax"]) if Z.any() else 0.0
                    p4_checked += 1
                    if res > bnd or d["Z"] > d["kelley"] + 1e-9:
                        p4_viol += 1
                    d["res_L"] = int((Z & ~VL).sum())
                rec["nodes"][name] = d
            fh.write(json.dumps(rec) + "\n")
    return {"n_primes": len(primes), "P3_violations": p3_viol, "P4_violations": p4_viol,
            "P4_checked": p4_checked}


def summarize_a(b: int, raw: str):
    recs = [json.loads(l) for l in open(raw)]
    nest, summ, ctrl, affine = node_tables()
    out = {}
    for name in affine:
        Zs = np.array([r["nodes"][name]["Z"] for r in recs])
        VL = np.array([r["nodes"][name]["VL"] for r in recs])
        d = {"mean_Z": float(Zs.mean()), "max_Z": int(Zs.max()), "primes_with_VL": int((VL > 0).sum())}
        if name in nest or name == "C_P":
            mod = np.array([r["nodes"][name]["model"] for r in recs])
            if name in nest:
                fl = np.array([r["nodes"][name]["flag"] for r in recs])
            else:
                fl = np.zeros(len(recs), dtype=bool)
            d["n_flagged"] = int(fl.sum())
            d["ratio_unflagged"] = float(Zs[~fl].sum() / mod[~fl].sum())
            d["ratio_flagged"] = float(Zs[fl].sum() / mod[fl].sum()) if fl.any() and mod[fl].sum() > 0 else None
            d["unflagged_with_VL"] = int(((VL > 0) & ~fl).sum())
            d["P5_pass"] = bool(0.75 <= d["ratio_unflagged"] <= 1.33 and d["unflagged_with_VL"] == 0)
            d["max_excess_unflagged"] = float((Zs[~fl] - mod[~fl]).max())
            if name in nest:   # POST HOC
                oz = np.array([r["nodes"][name]["other_Z"] for r in recs])
                om = np.array([r["nodes"][name]["other_model"] for r in recs])
                fz = np.array([r["nodes"][name]["Fr_Z"] for r in recs])
                d["posthoc_frozen_Z_share"] = float(fz.sum() / max(Zs.sum(), 1))
                d["posthoc_other_ratio_all"] = float(oz.sum() / om.sum())
                d["posthoc_other_ratio_unflagged"] = float(oz[~fl].sum() / om[~fl].sum())
                d["posthoc_other_max_excess"] = float((oz - om).max())
                d["posthoc_other_total"] = int(oz.sum())
                d["posthoc_other_model_total"] = float(om.sum())
        if name in summ:
            res = np.array([r["nodes"][name]["res_L"] for r in recs])
            d["resL_mean"] = float(res.mean())
            d["resL_max"] = int(res.max())
            d["ks_res_max"] = int(max(r["nodes"][name]["ks_res"] for r in recs))
            d["ks_bound_min"] = float(min(r["nodes"][name]["ks_bound"] for r in recs))
            if name != "S2":
                d["P7_pass"] = bool(0.3 <= d["resL_mean"] <= 3 and d["resL_max"] <= 20)
        out[name] = d
    return out


# ----------------------------------------------------------------------------- Part B
def dyadic_z(w: np.ndarray, P: np.ndarray, Q: np.ndarray):
    """z of sum(w) with the dyadic-robust variance (observations dependent iff they share a prime), residuals
    centred at their mean (as in leakmap v2).  Returns (z_dyadic, z_iid); z_dyadic is nan if the estimate is <= 0."""
    c = w - w.mean()
    S = defaultdict(float)
    for wi, a, b_ in zip(c, P, Q):
        S[a] += wi
        S[b_] += wi
    V = sum(v * v for v in S.values()) - float((c * c).sum())
    Viid = float((c * c).sum())
    zd = float(w.sum() / np.sqrt(V)) if V > 0 else float("nan")
    return zd, float(w.sum() / np.sqrt(Viid)) if Viid > 0 else float("nan")


def crt_value(primes, seed=0):
    """c = QR mod every eps=+1 prime, QR mod every other eps=-1 prime (sorted), non-residue on the rest."""
    mods, res = [], []
    k = 0
    for x in primes:
        if x % 4 == 1:
            want = 1
        else:
            want = 1 if k % 2 == 0 else -1
            k += 1
        if want == 1:
            a = 1
        else:
            a = 2
            while gmpy2.jacobi(a, x) != -1:
                a += 1
        mods.append(x)
        res.append(a)
    from sympy.ntheory.modular import crt
    return int(crt(mods, res)[0])


def part_b(b: int, npairs: int | None):
    nest, summ, ctrl, affine = node_tables()
    primes = primes_b(b)
    all_pairs = [(p, q) for i, p in enumerate(primes) for q in primes[i + 1:]]
    if npairs is not None and npairs < len(all_pairs):
        rng = random.Random(f"primE1|{b}")
        pairs = rng.sample(all_pairs, npairs)
    else:
        pairs = all_pairs
    names = list(affine) + [k for k, (nd, _) in nest.items() if not nd.affine()] + ["L1"]
    exprs = {**affine, **{k: nd for k, (nd, _) in nest.items() if not nd.affine()}}
    jv = jacobi_values(b)
    jnames = list(jv)
    crt_c = crt_value(primes) if b <= 14 else None
    if crt_c is not None:
        jnames.append("J_CRT")
    X = 1 << (b - 1)  # as in leakmap E1: target q - p < X/8 = 2^(b-4)
    nP = len(pairs)
    ap = np.zeros((nP, len(names)), dtype=np.int8)
    aq = np.zeros((nP, len(names)), dtype=np.int8)
    mismatch = Counter()
    accidental = 0
    J = np.zeros((nP, len(jnames)), dtype=np.int8)
    nest_names = [k for k in nest]
    mp = np.zeros((nP, len(nest_names)))
    mq = np.zeros((nP, len(nest_names)))
    for i, (p, q) in enumerate(pairs):
        N = p * q
        for j, nm in enumerate(names):
            if nm == "L1":
                vN = LucasNode.ev(N, N)
                vp = LucasNode.ev(N, p)
                vq = LucasNode.ev(N, q)
                zp = vp == (0, 0)
                zq = vq == (0, 0)
                g = gcd(gcd(vN[0], vN[1]), N)
                zN = vN == (0, 0)
            else:
                ex = exprs[nm]
                vN = ex.ev(N, N)
                zp = ex.ev(N, p) == 0
                zq = ex.ev(N, q) == 0
                g = gcd(vN, N)
                zN = vN == 0
            ap[i, j] = zp
            aq[i, j] = zq
            split = 1 < g < N
            # diagonal lemma: split iff zp != zq; else zN == zp
            if split != (zp != zq) or (not split and zN != zp):
                mismatch[nm] += 1
            if split and zp == zq:
                accidental += 1
        for j, nm in enumerate(nest_names):
            model = nest[nm][1]
            mp[i, j] = model.prob_scalar(N, p)
            mq[i, j] = model.prob_scalar(N, q)
        if N % 4 == 1:
            for j, nm in enumerate(jnames):
                c = (crt_c % N) if nm == "J_CRT" else jv[nm](N)
                J[i, j] = int(gmpy2.jacobi(c, N))
    P = np.array([p for p, _ in pairs])
    Q = np.array([q for _, q in pairs])
    return dict(pairs=pairs, P=P, Q=Q, names=names, ap=ap, aq=aq, mismatch=mismatch, accidental=accidental,
                J=J, jnames=jnames, nest_names=nest_names, mp=mp, mq=mq, X=X)


def batteries(names):
    nest = [k for k in names if k.startswith("D")]
    summ = [k for k in names if k.startswith("S")]
    ctrl = [k for k in names if k.startswith("C_")]
    return {"NEST": nest, "SUM": summ, "CTRL": ctrl, "ALL": [k for k in names if k != "L1"], "LUC": ["L1"]}


def analyse_b(B):
    names = B["names"]
    P, Q = B["P"], B["Q"]
    ap, aq = B["ap"], B["aq"]
    nP = len(P)
    # per-node majority answer per prime over its partners (both roles)
    maj = {}
    for j, nm in enumerate(names):
        cnt = defaultdict(lambda: [0, 0])
        for x, z in zip(P, ap[:, j]):
            cnt[x][0] += int(z); cnt[x][1] += 1
        for x, z in zip(Q, aq[:, j]):
            cnt[x][0] += int(z); cnt[x][1] += 1
        maj[nm] = {x: (c[0] * 2 > c[1]) for x, c in cnt.items()}
    devp = np.zeros_like(ap, dtype=bool)
    devq = np.zeros_like(aq, dtype=bool)
    for j, nm in enumerate(names):
        mj = maj[nm]
        devp[:, j] = ap[:, j] != np.array([mj[x] for x in P])
        devq[:, j] = aq[:, j] != np.array([mj[x] for x in Q])
    eps = np.array([1 if p % 4 == 1 else 0 for p in P])
    r3 = np.array([1 if p % 3 == 1 else 0 for p in P])
    X = B["X"]
    arch = ((Q - P) < X // 8).astype(int)
    targets = {"eps": eps, "r3": r3, "A": arch}
    sub = {"eps": (P * Q) % 4 == 1, "r3": (P * Q) % 3 == 1, "A": np.ones(nP, dtype=bool)}
    cells = {}
    p2_fail = 0
    for bname, nodes in batteries(names).items():
        idx = [names.index(k) for k in nodes]
        tp = [tuple(row) for row in ap[:, idx]]
        tq = [tuple(row) for row in aq[:, idx]]
        split = np.array([a != c for a, c in zip(tp, tq)])
        D = devp[:, idx].any(axis=1) | devq[:, idx].any(axis=1)
        for tname, F in targets.items():
            for popname, mask in (("all", np.ones(nP, dtype=bool)), ("sub", sub[tname])):
                if popname == "sub" and tname == "A":
                    continue
                m = mask
                n_ = int(m.sum())
                s = float(split[m].mean())
                Bv = max(F[m].mean(), 1 - F[m].mean())
                groups = defaultdict(lambda: [0, 0])
                for t, f, sp in zip(itertools.compress(tp, m), F[m], split[m]):
                    if not sp:
                        groups[t][int(f)] += 1
                correct_ns = sum(max(g) for g in groups.values()) / n_
                pdns = float((D[m] & ~split[m]).mean())
                lhs = correct_ns + s
                if popname == "all":
                    rhs = Bv + 2 * s + 2 * pdns
                    ok = lhs <= rhs + 1e-12
                    if not ok:
                        p2_fail += 1
                else:
                    rhs = 1.0  # fiber form: the target is a function of the fiber point, so the bound is 1
                    ok = True
                cells[f"{bname}|{tname}|{popname}"] = {
                    "n": n_, "s": s, "B": Bv, "correct_or_split": lhs, "Pr_D_nosplit": pdns, "rhs": rhs,
                    "holds": bool(ok), "nonsplit_gain": correct_ns - (1 - s) * Bv}
    return cells, p2_fail, devp, devq


def p6(B):
    names = B["names"]
    out = {}
    tot_two = 0
    tot_two_model = 0.0
    for j, nm in enumerate(B["nest_names"]):
        jj = names.index(nm)
        zp = B["ap"][:, jj].astype(bool)
        zq = B["aq"][:, jj].astype(bool)
        mp, mq = B["mp"][:, j], B["mq"][:, j]
        one = int(zp.sum() + zq.sum())
        model = float(mp.sum() + mq.sum())
        two = int((zp & zq).sum())
        two_m = float((mp * mq).sum())
        affine = not nm.startswith("D2f")
        out[nm] = {"one_sided": one, "model": model, "ratio": one / model if model > 0 else None,
                   "two_sided": two, "two_sided_model": two_m, "affine": affine,
                   "P6_ratio_pass": (0.6 <= one / model <= 1.6) if (affine and model > 0) else None}
        if affine:
            tot_two += two
            tot_two_model += two_m
    out["_two_sided_total"] = {"observed": tot_two, "model": tot_two_model,
                               "pass": bool(tot_two <= 3 * tot_two_model + 5)}
    return out


def p8(B):
    out = {}
    P, Q = B["P"], B["Q"]
    m = (P * Q) % 4 == 1
    eps = np.where(P[m] % 4 == 1, 1.0, -1.0)
    for j, nm in enumerate(B["jnames"]):
        Jv = B["J"][m, j].astype(float)
        w = (Jv - Jv.mean()) * (eps - eps.mean())
        z, ziid = dyadic_z(w, P[m], Q[m])
        out[nm] = {"n": int(m.sum()), "split": int((Jv == 0).sum()), "mean_J_eps+": float(Jv[eps > 0].mean()),
                   "mean_J_eps-": float(Jv[eps < 0].mean()), "z": z, "z_iid": ziid,
                   "P8_pass": (abs(z) < 4) if nm != "J_CRT" else (z > 10)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--b", type=int, required=True)
    ap.add_argument("--pairs", type=int, default=None)
    ap.add_argument("--skip-a", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    raw = os.path.join(ROOT, f"experiments/raw/prim_e1A_b{a.b}.jsonl")
    res = {"b": a.b}
    if not a.skip_a:
        res["partA"] = part_a(a.b, raw)
        print("part A done", time.time() - t0, "s, rss", rss_mb(), flush=True)
    res["partA_summary"] = summarize_a(a.b, raw)
    B = part_b(a.b, a.pairs)
    print("part B eval done", time.time() - t0, "s, rss", rss_mb(), flush=True)
    res["n_pairs"] = len(B["pairs"])
    res["P1_mismatches"] = dict(B["mismatch"])
    res["P1_total_mismatches"] = int(sum(B["mismatch"].values()))
    res["accidental_splits"] = B["accidental"]
    cells, p2_fail, devp, devq = analyse_b(B)
    res["P2_cells"] = cells
    res["P2_failures"] = p2_fail
    res["P6"] = p6(B)
    res["P8"] = p8(B)
    res["peak_rss_mb"] = rss_mb()
    res["seconds"] = time.time() - t0
    out = os.path.join(ROOT, f"experiments/analysis/prim_e1_b{a.b}.json")
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1, default=lambda o: o if not isinstance(o, np.generic) else o.item())
    print(json.dumps({k: res[k] for k in ("b", "n_pairs", "P1_total_mismatches", "accidental_splits",
                                           "P2_failures", "seconds", "peak_rss_mb")}))


if __name__ == "__main__":
    main()
