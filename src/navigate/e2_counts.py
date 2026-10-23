"""E2 (pre-registered, notes/navigate_prereg.md): exact completion counts of signed smooth exponent vectors per class.
HARNESS / DIAGNOSTIC: uses p, q to build discrete-log coordinates for (Z/N)^x (Fourier transform, subgroup H, squares).

C_beta(t) = #{e in Z^FB : sum |e_i| ln p_i <= beta ln N, prod p_i^e_i = t (mod N)}  (= coprime smooth pairs u/v = t).

    .venv/bin/python -m src.navigate.e2_counts
"""
from __future__ import annotations

import json
import math
import os
import random
import time
from typing import Dict, List, Tuple

import numpy as np
import sympy

from src.common.instances import balanced_semiprime
from src.navigate.core import primes_upto

SEEDS = range(4)
BITS = 20
CONFIGS = {13: [1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0], 19: [1.0, 1.25, 1.5, 1.75, 2.0]}
STEP = 0.05          # fine budget grid (units of ln N) for the descriptive oracle-walk profile
LO = 0.5
OUT = "experiments/raw/navigate_E2.jsonl"


def dlog_table(p: int) -> Tuple[int, np.ndarray]:
    g = sympy.primitive_root(p)
    tab = np.zeros(p, dtype=np.int64)
    x = 1
    for k in range(p - 1):
        tab[x] = k
        x = x * g % p
    return g, tab


def signed_head(lns: List[float], ci: List[int], cj: List[int], P1: int, Q1: int, B: float):
    """All signed exponent vectors over the head primes with weight <= B: arrays (w, i, j, vec) sorted by w."""
    W = np.zeros(1)
    I = np.zeros(1, dtype=np.int64)
    J = np.zeros(1, dtype=np.int64)
    for k in range(len(lns)):
        parts_w, parts_i, parts_j = [W], [I], [J]
        e = 1
        while e * lns[k] <= B + 1e-9:
            for s in (1, -1):
                m = W + e * lns[k] <= B + 1e-9
                parts_w.append(W[m] + e * lns[k])
                parts_i.append((I[m] + s * e * ci[k]) % P1)
                parts_j.append((J[m] + s * e * cj[k]) % Q1)
            e += 1
        W, I, J = np.concatenate(parts_w), np.concatenate(parts_i), np.concatenate(parts_j)
    o = np.argsort(W, kind="stable")
    return W[o], I[o], J[o]


def tail_vectors(lns, ci, cj, P1, Q1, B):
    out = [(0.0, 0, 0)]
    for k in range(len(lns)):
        new = []
        for (w, i, j) in out:
            new.append((w, i, j))
            e = 1
            while w + e * lns[k] <= B + 1e-9:
                for s in (1, -1):
                    new.append((w + e * lns[k], (i + s * e * ci[k]) % P1, (j + s * e * cj[k]) % Q1))
                e += 1
        out = new
    return out


def tail_vectors_sparse(lns, ci, cj, P1, Q1, B):
    """All signed vectors over the given primes with weight <= B (depth-first, primes in increasing order)."""
    out = []
    n = len(lns)

    def rec(k, w, i, j):
        out.append((w, i, j))
        for kk in range(k, n):
            l = lns[kk]
            if w + l > B + 1e-9:
                break  # lns increasing
            e = 1
            while w + e * l <= B + 1e-9:
                for s in (1, -1):
                    rec(kk + 1, w + e * l, (i + s * e * ci[kk]) % P1, (j + s * e * cj[kk]) % Q1)
                e += 1

    rec(0, 0.0, 0, 0)
    return out


def subgroup_mask(P1: int, Q1: int, gens: List[Tuple[int, int]]) -> np.ndarray:
    m = np.zeros((P1, Q1), dtype=bool)
    m[0, 0] = True
    while True:
        old = m.sum()
        for (a, b) in gens:
            m |= np.roll(np.roll(m, a, axis=0), b, axis=1)
        if m.sum() == old:
            return m


def run_instance(seed: int, y: int, betas: List[float], exclude_nonunits: bool = False) -> Dict[str, object]:
    inst = balanced_semiprime(BITS, seed)
    N, p, q = inst.N, inst.labels["p"], inst.labels["q"]
    P1, Q1 = p - 1, q - 1
    lnN = math.log(N)
    fb = primes_upto(y)
    if exclude_nonunits:  # POST HOC repair (codex navigate_2 #12): drop factor-base primes dividing N
        fb = [l for l in fb if N % l]
    assert all(N % l for l in fb), "factor base must consist of units mod N"
    _, tp = dlog_table(p)
    _, tq = dlog_table(q)
    ci = [int(tp[l % p]) for l in fb]
    cj = [int(tq[l % q]) for l in fb]
    lns = [math.log(l) for l in fb]
    bmax = max(betas)
    B = bmax * lnN
    nlev = int(round((bmax - LO) / STEP)) + 1
    hist = np.zeros((nlev, P1 * Q1), dtype=np.int32)  # hist[k] = vectors with weight in level k (see below)
    # head = vectorised over the small primes (dense exponents); tails = Python enumeration over the large, sparse
    # primes (memory: v1 vectorised all but two primes, which fragmented to 2.4-4.9 GB RSS at K ~ 100-300).
    nh = max(1, min(len(fb) - 2, sum(1 for l in fb if l <= 47)))
    W, I, J = signed_head(lns[:nh], ci[:nh], cj[:nh], P1, Q1, B)
    tails = tail_vectors_sparse(lns[nh:], ci[nh:], cj[nh:], P1, Q1, B)
    total = 0
    for (tw, ti, tj) in tails:
        n = int(np.searchsorted(W, B - tw + 1e-9, side="right"))
        if n == 0:
            continue
        for c0 in range(0, n, 1 << 21):  # chunked to bound memory
            c1 = min(n, c0 + (1 << 21))
            w = W[c0:c1] + tw
            cls = ((I[c0:c1] + ti) % P1) * Q1 + (J[c0:c1] + tj) % Q1
            lev = np.ceil((w / lnN - LO) / STEP - 1e-9).astype(np.int64)
            lev = np.clip(lev, 0, nlev - 1)  # weights <= LO*lnN land in level 0 (counted at every budget >= LO)
            np.add.at(hist, (lev, cls), 1)
        total += n
    cum = np.cumsum(hist, axis=0, dtype=np.int32)  # cum[k] = C at budget LO + k*STEP
    Hm = subgroup_mask(P1, Q1, list(zip(ci, cj))).reshape(-1)
    ii, jj = np.divmod(np.arange(P1 * Q1), Q1)
    sq = ((ii % 2) == 0) & ((jj % 2) == 0)
    # characters trivial on H: (a*ci*Q1 + b*cj*P1) = 0 mod P1*Q1 for all generators
    A, Bc = np.divmod(np.arange(P1 * Q1, dtype=np.int64), Q1)
    triv = np.ones(P1 * Q1, dtype=bool)
    for (c1, c2) in zip(ci, cj):
        triv &= ((A * c1 * Q1 + Bc * c2 * P1) % (P1 * Q1)) == 0
    res = {"seed": seed, "y": y, "K": len(fb), "excluded_nonunits": exclude_nonunits, "N": N, "phi": P1 * Q1, "H": int(Hm.sum()), "H_has_nonsquares": bool((Hm & ~sq).any()),
           "vectors": total, "cells": []}
    for beta in betas:
        k = int(round((beta - LO) / STEP))
        C = cum[k]
        CH = C[Hm].astype(np.float64)
        mu = CH.mean()
        h = float((CH >= 1).mean())
        D = float(CH.var() / mu) if mu > 0 else None
        F = np.fft.fft2(C.reshape(P1, Q1).astype(np.float64)).reshape(-1)
        mass = np.abs(F[~triv]) ** 2
        mass = np.sort(mass)[::-1]
        cs = np.cumsum(mass)
        k75 = int(np.searchsorted(cs, 0.75 * cs[-1]) + 1) if cs[-1] > 0 else None
        cell = {"beta": beta, "mu": float(mu), "h": h, "poisson_h": float(1 - math.exp(-mu)), "D": D,
                "k75": k75, "nchars": int(mass.size), "f75": (k75 / mass.size) if k75 else None,
                "max": int(CH.max())}
        if res["H_has_nonsquares"]:
            s1, s0 = Hm & sq, Hm & ~sq
            cell.update({"h_sq": float((C[s1] >= 1).mean()), "h_nonsq": float((C[s0] >= 1).mean()),
                         "mu_sq": float(C[s1].mean()), "mu_nonsq": float(C[s0].mean())})
        res["cells"].append(cell)
    # descriptive oracle-walk profile: draw a uniform vector of height <= B and keep only its class (a size-biased
    # target); then run an oracle-guided walk that, at each step, picks uniformly among the (prime, sign) letters
    # over ALL factor-base primes whose remaining instance has a nonzero count, recording the fraction sigma of such
    # letters (budgets rounded down to the grid: conservative). No -1 letter (the E2 factor base has no sign).
    rng = random.Random(f"navigate-E2-walk|{seed}|{y}")
    prof: Dict[str, List[float]] = {}
    gi = dict(zip(fb, ci)), dict(zip(fb, cj))
    for _ in range(200):
        # sample a vector with weight <= B uniformly: rejection from head x tail
        while True:
            tw, ti, tj = tails[rng.randrange(len(tails))]
            n = int(np.searchsorted(W, B - tw + 1e-9, side="right"))
            if n and rng.random() < n / len(W):
                break
        # recover an explicit vector: resample head index and rebuild exponents is costly; instead walk on classes:
        # we only need the class trajectory and budgets, which we get by peeling letters from an explicit vector.
        hidx = rng.randrange(n)
        target_i, target_j = (int(I[hidx]) + ti) % P1, (int(J[hidx]) + tj) % Q1
        wtot = float(W[hidx] + tw)
        # explicit exponents: decode by brute force over the head primes is unnecessary for the profile; we peel
        # letters greedily from the class side by choosing, at each step, a letter with nonzero completion count
        # that keeps the walk on a solution (oracle step), recording sigma.
        ci_, cj_ = target_i, target_j
        rem = bmax
        steps = 0
        while rem > LO + 1e-9 and steps < 200 and (ci_, cj_) != (0, 0):
            k = int(math.floor((rem - LO) / STEP + 1e-9))
            cnt_here = cum[k][ci_ * Q1 + cj_]
            if cnt_here == 0:
                break
            good = []
            tot = 0
            for l in fb:
                dl = math.log(l) / lnN
                if rem - dl < LO - 1e-9:
                    continue
                k2 = int(math.floor((rem - dl - LO) / STEP + 1e-9))
                for s in (1, -1):
                    tot += 1
                    a2, b2 = (ci_ - s * gi[0][l]) % P1, (cj_ - s * gi[1][l]) % Q1
                    if cum[k2][a2 * Q1 + b2] > 0:
                        good.append((l, s, a2, b2, dl))
            if tot == 0 or not good:
                break
            prof.setdefault(f"{round(rem, 2):.2f}", []).append(len(good) / tot)
            l, s, a2, b2, dl = good[rng.randrange(len(good))]
            ci_, cj_, rem = a2, b2, rem - dl
            steps += 1
    res["walk_profile"] = {k: [float(np.mean(v)), len(v)] for k, v in prof.items()}
    return res


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        for y, betas in CONFIGS.items():
            for seed in SEEDS:
                t0 = time.time()
                r = run_instance(seed, y, betas)
                fh.write(json.dumps(r) + "\n")
                fh.flush()
                print(y, seed, r["vectors"], r["H"], r["phi"], round(time.time() - t0, 1), flush=True)


if __name__ == "__main__":
    main()
