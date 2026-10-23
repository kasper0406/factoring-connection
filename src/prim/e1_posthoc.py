"""POST HOC analysis of E1 (not pre-registered; see notes/prim_results.md §7 "deviations").  DIAGNOSTIC.

Registered failures of P5/P6/P7 were all traced to two things the registered fresh-random model ignores:
  * FROZEN zeros: r such that the gate outputs at r coincide with those at a small integer t (|t| <= 64), where the
    integer tower T(t) vanishes or has p as a factor (e.g. 2 + 3 = 5 for S1 at t = 1);  F1 is the weaker "layer-1
    frozen" class (only gates with a pure constant base are frozen at t), on which the node is a ONE-layer node with
    small-integer bases;
  * the shift exclusion: u = (gate) + c never equals c, and for E = N the gate is a permutation at p
    (gcd(q, p-1) = 1), so A*u^e = B may have its only solution at u = c.
This script reclassifies every zero as COSET (registered), FROZEN, F1, PIN (registered), OTHER and compares OTHER
with the shift-corrected model, in the full-r world (all r mod p-1) and in the q-world over ALL pairs (Theorem U:
the answer at p for partner q is Z_p[q mod (p-1)]).

Usage: python -m src.prim.e1_posthoc --b 14   ->  experiments/analysis/prim_e1_posthoc_b{b}.json
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from .prog import C, G, PPow, nested_nodes, summed_nodes, vpow
from .zerosets import L_cube, coset_union, frozen_set, gate_subexprs, order_mod, primes_b

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

SHIFT = {"D2a_1": 1, "D2a_2": 1, "D2a_3": 1, "D2b": 1, "D2c": 2, "D2d": 3, "D2e": None, "D3": 1}


def layer1_frozen(ex, p: int, T: int = 64) -> np.ndarray:
    n = p - 1
    r = np.arange(n, dtype=np.int64)
    gates = [g for g in gate_subexprs(ex) if isinstance(g, G) and isinstance(g.x, C)]
    F = np.zeros(n, dtype=bool)
    if not gates:
        return F
    outs = [g.vec(p, r) for g in gates]
    for t in range(-T, T + 1):
        m = np.ones(n, dtype=bool)
        for o in outs:
            m &= o == o[t % n]
        F |= m
    return F


def node_constants(ex) -> set:
    from .prog import Add, Mul, Sub
    if isinstance(ex, C):
        return {ex.c}
    if isinstance(ex, (Add, Sub, Mul)):
        return node_constants(ex.x) | node_constants(ex.y)
    if hasattr(ex, "x"):
        return node_constants(ex.x)
    return set()


def closure_set(ex, p: int) -> np.ndarray:
    """POST HOC (second round) class CLOSURE: r such that some non-constant gate input (an intermediate value)
    equals s * prod y^{k_y} with s = +-1, y over the node's integer constants (|y| >= 2) and sum |k_y| <= 2.
    On such r the outer gate's output is a monomial in lower-depth gate outputs (the tower closes), e.g.
    a^r = -(a+1) gives (a^r + 1)^r = (-a)^r = a + 1 for odd r."""
    n = p - 1
    r = np.arange(n, dtype=np.int64)
    consts = sorted({abs(c) for c in node_constants(ex) if abs(c) >= 2})
    S = {1 % p, (p - 1) % p}
    mons = [1]
    for y in consts:
        for k in (1, 2):
            for e in (k, -k):
                mons.append(pow(y, e, p))
    prods = set(mons)
    for i, a in enumerate(consts):
        for b_ in consts[i + 1:]:
            for ea in (1, -1):
                for eb in (1, -1):
                    prods.add(pow(a, ea, p) * pow(b_, eb, p) % p)
    for m in prods:
        S.add(m % p)
        S.add((-m) % p)
    Sarr = np.array(sorted(S), dtype=np.int64)
    out = np.zeros(n, dtype=bool)
    for g in gate_subexprs(ex):
        if isinstance(g.x, C):
            continue
        out |= np.isin(g.x.vec(p, r), Sarr)
    return out


def shift_model(model, p: int, r: np.ndarray, c) -> np.ndarray:
    """Pr_{u uniform on F_p^* minus {c}}[A u^e = B]  (c None: plain model)."""
    if c is None:
        return model.prob_vec(p, r)
    n = p - 1
    A = model.A.vec(p, r)
    B = model.B.vec(p, r)
    e = model.exponent_mod(None, p, r)
    d = np.gcd(e, n)
    Ainv = vpow(A, np.full(r.shape, p - 2, dtype=np.int64), p)
    t = (B * Ainv) % p
    ok = (vpow(t, n // d, p) == 1) & (t != 0)
    nsol = np.where(ok, d, 0).astype(float)
    csol = (vpow(np.full(r.shape, c % p), e, p) == t) & ok
    nsol = nsol - csol
    out = nsol / (n - 1)
    return np.where(A == 0, (B == 0).astype(float), out)


def run(b: int):
    nest = {k: v for k, v in nested_nodes().items() if v[0].affine()}
    summ = summed_nodes()
    primes = primes_b(b)
    parr = np.array(primes, dtype=np.int64)
    idx = {x: i for i, x in enumerate(primes)}
    half = 1 << (b // 2)
    res = {}
    names = list(nest) + list(summ)
    # per node: hit matrix H[i, j] = node zero at primes[i] for partner primes[j]; class matrix
    for name in names:
        nd = nest[name][0] if name in nest else summ[name][0]
        model = nest[name][1] if name in nest else None
        K = len(primes)
        H = np.zeros((K, K), dtype=bool)
        CL = np.zeros((K, K), dtype=np.int8)     # 0 none, 1 COSET, 2 FROZEN, 3 F1, 4 PIN, 5 OTHER
        Mq = np.zeros((K, K))
        full = {"Z": 0, "COSET": 0, "FROZEN": 0, "F1": 0, "OTHER": 0, "model_other": 0.0, "model_all": 0.0,
                "flag_primes": 0}
        for i, p in enumerate(primes):
            n = p - 1
            r = np.arange(n, dtype=np.int64)
            Z = nd.vec(p, r) == 0
            VL = coset_union(Z, L_cube(n))
            Fr = frozen_set(nd, p)
            F1 = layer1_frozen(nd, p)
            CLO = closure_set(nd, p)
            other = Z & ~VL & ~Fr & ~F1 & ~CLO
            full["Z"] += int(Z.sum())
            full["COSET"] += int((Z & VL).sum())
            full["FROZEN"] += int((Z & Fr & ~VL).sum())
            full["F1"] += int((Z & F1 & ~Fr & ~VL).sum())
            full["CLOSURE"] = full.get("CLOSURE", 0) + int((Z & CLO & ~F1 & ~Fr & ~VL).sum())
            full["OTHER"] += int(other.sum())
            if model is not None:
                mv = shift_model(model, p, r, SHIFT[name])
                expl = VL | Fr | F1 | CLO
                full["model_other"] += float(mv[~expl].sum())
                full["model_all"] += float(mv.sum())
            rq = parr % n
            H[i] = Z[rq]
            H[i, i] = False
            cl = np.where(VL, 1, np.where(Fr, 2, np.where(F1, 3, np.where(CLO, 6, 5))))[rq]
            pin = (np.abs(parr - p) <= half) | (rq <= 64)
            cl = np.where((cl == 5) & pin, 4, cl)
            CL[i] = np.where(H[i], cl, 0)
            if model is not None:
                mv_q = mv[rq]
                mv_q = np.where(VL | Fr | F1 | CLO, 0.0, mv)[rq]   # model mass on OTHER residues only
                mv_q = np.where(pin, 0.0, mv_q)                     # v2 (codex prim_2 #20): mask PIN as in numerator
                Mq[i] = mv_q
                Mq[i, i] = 0.0
        two = H & H.T
        np.fill_diagonal(two, False)
        cls_names = {1: "COSET", 2: "FROZEN", 3: "F1", 4: "PIN", 5: "OTHER", 6: "CLOSURE"}
        qw = {cls_names[k]: int((CL == k).sum()) for k in cls_names}
        # two-sided (non-split zero) events by the pair of classes, counted once per unordered pair
        iu = np.triu_indices(len(primes), 1)
        tw = two[iu]
        c1, c2 = CL[iu][tw], CL.T[iu][tw]
        two_cls = {}
        for a_, b_ in zip(c1, c2):
            key = "/".join(sorted((cls_names[int(a_)], cls_names[int(b_)])))
            two_cls[key] = two_cls.get(key, 0) + 1
        entry = {"full_r": full, "q_world_one_sided": qw, "q_world_two_sided": int(tw.sum()),
                 "two_sided_by_class": two_cls}
        if model is not None:
            entry["full_r_other_ratio"] = full["OTHER"] / full["model_other"] if full["model_other"] else None
            entry["q_world_other_model"] = float(Mq.sum())
            entry["q_world_other_ratio"] = qw["OTHER"] / Mq.sum() if Mq.sum() else None
            entry["q_world_two_sided_other_model"] = float((Mq * Mq.T)[iu].sum())
        res[name] = entry
        print(name, json.dumps(entry), flush=True)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--b", type=int, required=True)
    a = ap.parse_args()
    t0 = time.time()
    res = run(a.b)
    res["_seconds"] = time.time() - t0
    with open(os.path.join(ROOT, f"experiments/analysis/prim_e1_posthoc_b{a.b}.json"), "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
