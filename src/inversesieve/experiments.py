"""Pre-registered experiments E1-E4 (notes/inversesieve_prereg.md).  HARNESS code: E3/E4 use factors.

Run:  .venv/bin/python -m src.inversesieve.experiments e1|e2|e3|e4|all [--workers 8]
Outputs: experiments/raw/inversesieve_e{1,2,3,4}.jsonl
"""
from __future__ import annotations

import json
import os
import random
import sys
import time
from math import gcd, isqrt, log
from multiprocessing import Pool
from typing import Dict, List, Tuple

import gmpy2
import numpy as np

from .design import band, weights, forbidden_set, farey_set, candidate_interval
from .usets import U_of_N, uncovered_partners, is_pair_uncovered
from . import stats as st

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "experiments", "raw")

NS = (24, 28, 32, 36, 40)
K_TARGET = 300          # H = least height with E_model|U| <= K_TARGET (signed weights incl. (1,1))
N_PER_CONFIG = 100
NULL_DRAWS = 20
AP_DRAWS = 5
AP_SUB = 200
D_PAIR = 32


def shapes(n: int) -> Dict[str, int]:
    """Band scale log2 Y for the two frozen shapes."""
    return {"M": (n + 3) // 4 + 1, "S": n // 2 - 4}


def choose_H(n: int, Y: int, primes: List[int]) -> int:
    P = 2.0 ** (n / 2)
    I = P * (1 - 2 ** -0.5)
    for H in range(2, 200):
        if 2 * H * H >= Y:
            return H - 1
        V = weights(H, "signed")
        d = 1.0
        for l in primes:
            d *= 1 - len(V) / l
        if I * d <= K_TARGET:
            return H
    raise RuntimeError


def config(n: int, shape: str) -> dict:
    logY = shapes(n)[shape]
    Y = 1 << logY
    pr = band(Y)
    H = choose_H(n, Y, pr)
    return {"n": n, "shape": shape, "logY": logY, "Y": Y, "J": len(pr), "H": H,
            "nV": len(weights(H, "signed")), "primes": pr}


def configs():
    return [config(n, s) for n in NS for s in ("M", "S")]


# ---------------------------------------------------------------------------
# instance sampling (seeded; approx. uniform over pairs of S_n)
# ---------------------------------------------------------------------------
def q_range(p: int, n: int) -> Tuple[int, int]:
    lo = max(p + 1, -(-(1 << (n - 1)) // p))
    hi = min(2 * p, ((1 << n) - 1) // p)
    return lo, hi


def p_range(n: int) -> Tuple[int, int]:
    return isqrt(1 << (n - 2)) + 1, isqrt((1 << n) - 1)


def sample_pair(rng: random.Random, n: int) -> Tuple[int, int]:
    plo, phi = p_range(n)
    Lmax = max(q_range(p, n)[1] - q_range(p, n)[0] + 1 for p in (plo, (plo + phi) // 2, int(phi / 2 ** 0.5)))
    while True:
        p = rng.randint(plo, phi)
        if not gmpy2.is_prime(p):
            continue
        lo, hi = q_range(p, n)
        if hi < lo or rng.random() * Lmax > hi - lo + 1:
            continue
        while True:
            q = rng.randint(lo, hi)
            if gmpy2.is_prime(q):
                return p, q


def random_prime(rng: random.Random, lo: int, hi: int) -> int:
    while True:
        p = rng.randint(lo, hi)
        if gmpy2.is_prime(p):
            return p


# ---------------------------------------------------------------------------
# E1: structure of U(N)
# ---------------------------------------------------------------------------
def filter_with(x, primes, Rs):
    for l in primes:
        if x.size == 0:
            break
        m = np.zeros(l, dtype=bool); m[Rs[l]] = True
        x = x[~m[x % l]]
    return x


def survivor_trajectory(N, primes, V, a0, b0, Rs):
    x = np.arange(a0, b0 + 1, dtype=np.int64)
    traj, model, d = [], [], 1.0
    J = len(primes)
    marks = {max(1, round(J * f)) for f in (0.1, 0.25, 0.5, 0.75, 1.0)}
    for i, l in enumerate(primes, 1):
        R = Rs[l]
        m = np.zeros(l, dtype=bool); m[R] = True
        x = x[~m[x % l]]
        d *= 1 - len(R) / l
        if i in marks:
            traj.append([i, int(x.size)]); model.append((b0 - a0 + 1) * d)
    return x, traj, model


def admissible_sample(rng_np, a0, b0, Rs, primes_S, k):
    """Uniform k-subset of the integers of [a0,b0] admissible at the primes of S (null for CRT stats)."""
    x = filter_with(np.arange(a0, b0 + 1, dtype=np.int64), primes_S, Rs)
    if k > x.size:
        return x
    return np.sort(rng_np.choice(x, size=k, replace=False))


def crt_stats(U, primes_S):
    M = 1
    for l in primes_S:
        M *= l
    r = U % M
    _, c = np.unique(r, return_counts=True)
    return {"M": M, "distinct": int(len(c)), "collisions": int((c * (c - 1) // 2).sum()), "maxocc": int(c.max()) if len(c) else 0}


def e1_one(args):
    cfg, idx = args
    n, pr = cfg["n"], cfg["primes"]
    V = weights(cfg["H"], "signed")
    rng = random.Random(f"inversesieve|e1|n={n}|shape={cfg['shape']}|{idx}")
    rng_np = np.random.default_rng(n * 1000003 + idx * 7919 + (0 if cfg["shape"] == "M" else 1))
    p, q = sample_pair(rng, n)
    N = p * q                                   # only N is used below (U is operational)
    a0, b0 = candidate_interval(N)
    t0 = time.time()
    Rs = {l: forbidden_set(N, l, V) for l in pr}
    U, traj, traj_model = survivor_trajectory(N, pr, V, a0, b0, Rs)
    tU = time.time() - t0
    I = b0 - a0 + 1
    k = int(U.size)
    rec = {"n": n, "shape": cfg["shape"], "Y": cfg["Y"], "J": cfg["J"], "H": cfg["H"], "nV": len(V), "idx": idx,
           "N": N, "a0": a0, "b0": b0, "U_size": k, "E_model": traj_model[-1], "traj": traj, "traj_model": traj_model,
           "t_U": tU}
    sqrtP = isqrt(isqrt(N))                       # floor(N^{1/4}) = P^{1/2}
    nulls = [st.uniform_null(rng_np, a0, b0, k) for _ in range(NULL_DRAWS)] if k >= 2 else []
    # interval covering numbers
    cov = {}
    for tag, L in (("q", max(1, sqrtP // 4)), ("1", sqrtP), ("4", 4 * sqrtP)):
        cov[tag] = {"L": L, "obs": st.interval_cover_number(U, L),
                    "null": [st.interval_cover_number(z, L) for z in nulls]}
    rec["cover"] = cov
    # AP triples (Coppersmith residue containers p = r mod M, M >= P^{1/2})
    if k >= 3:
        sub = U if k <= AP_SUB else np.sort(rng_np.choice(U, size=AP_SUB, replace=False))
        ks = len(sub)
        rec["ap"] = {"Lmin": sqrtP, "k": ks, "obs": st.ap_triples(sub, sqrtP),
                     "null": [st.ap_triples(st.uniform_null(rng_np, a0, b0, ks), sqrtP) for _ in range(AP_DRAWS)]}
    # small-prime class occupancy (low-degree image detector)
    occ = {}
    for r in (3, 5, 7, 11, 13, 17, 19, 23):
        if k < 4 * r * log(r):
            continue
        o, c2 = st.occupancy(U, r)
        nul = [st.occupancy(z, r) for z in nulls]
        occ[str(r)] = {"obs_missing": r - o, "obs_chi2": c2,
                       "null_missing": [r - a for a, _ in nul], "null_chi2": [b for _, b in nul]}
    rec["occ"] = occ
    # gaps, pair correlation, energy
    rec["gaps"] = np.diff(U).tolist()
    rec["null_gaps"] = np.diff(nulls[0]).tolist() if nulls else []
    pc = st.pair_counts(U, D_PAIR)
    g = st.local_pair_density(N, pr, V, D_PAIR, Rs)
    rec["pair_obs"] = pc[1:].tolist()
    rec["pair_model"] = [(I - d) * g[d] for d in range(1, D_PAIR + 1)]
    rec["pair_uniform"] = [k * (k - 1) * (I - d) / (I * (I - 1)) for d in range(1, D_PAIR + 1)]
    if k >= 2:
        rec["energy"] = {"obs": st.energy_int(U), "null": [st.energy_int(z) for z in nulls[:10]]}
    # CRT concentration at 1, 2, 3 band primes vs admissible-conditioned null
    crt = {}
    for s in (1, 2, 3):
        S = pr[:s]
        nul = [crt_stats(admissible_sample(rng_np, a0, b0, Rs, S, k), S) for _ in range(10)] if k >= 2 else []
        crt[str(s)] = {"obs": crt_stats(U, S), "null": nul}
    rec["crt"] = crt
    # held-out halves (independence between disjoint halves of the band)
    h1, h2 = pr[0::2], pr[1::2]
    U1 = filter_with(np.arange(a0, b0 + 1, dtype=np.int64), h1, Rs)
    U2 = filter_with(np.arange(a0, b0 + 1, dtype=np.int64), h2, Rs)
    d1 = float(np.prod([1 - len(Rs[l]) / l for l in h1]))
    d2 = float(np.prod([1 - len(Rs[l]) / l for l in h2]))
    rec["halves"] = {"U1": int(U1.size), "U2": int(U2.size), "U": k, "I": I, "E1": I * d1, "E2": I * d2}
    # Lehman containers (N-only)
    h0 = cfg["h0"]
    cen = st.lehman_centres(N, h0)
    half = sqrtP // 2
    rec["lehman"] = {"h0": h0, "ncent": int(cen.size), "measure": st.container_measure(cen, half, a0, b0),
                     "in": int(st.in_containers(U, cen, half).sum()) if k else 0}
    return rec


def lehman_h0(n: int) -> int:
    """Least h0 whose Lehman containers (half-width P^{1/2}/2) cover >= 10% of the interval for the
    reference input N = 3 * 2^{n-2} (design-only)."""
    N = 3 << (n - 2)
    a0, b0 = candidate_interval(N)
    half = isqrt(isqrt(N)) // 2
    for h0 in range(1, 10 ** 4):
        if st.container_measure(st.lehman_centres(N, h0), half, a0, b0) >= 0.1 * (b0 - a0 + 1):
            return h0
    raise RuntimeError


def _with_h0(cfg):
    cfg = dict(cfg)
    cfg["h0"] = lehman_h0(cfg["n"])
    return cfg


def run(name, fn, jobs, workers):
    path = os.path.join(RAW, f"inversesieve_{name}.jsonl")
    if os.path.exists(path):
        os.remove(path)
    t0 = time.time()
    with Pool(workers) as pool, open(path, "w") as f:
        for i, rec in enumerate(pool.imap(fn, jobs, chunksize=1)):
            f.write(json.dumps(rec) + "\n"); f.flush()
            if (i + 1) % 50 == 0:
                print(f"{name}: {i + 1}/{len(jobs)} ({time.time() - t0:.0f}s)", flush=True)
    print(f"{name}: done {len(jobs)} records in {time.time() - t0:.0f}s -> {path}", flush=True)


def e1(workers):
    cf = [_with_h0(c) for c in configs()]
    jobs = [(c, i) for c in cf for i in range(N_PER_CONFIG)]
    run("e1", e1_one, jobs, workers)


# ---------------------------------------------------------------------------
# E2: additive / multiplicative energy of F_V(l) and R_l(N)
# ---------------------------------------------------------------------------
def e2_one(args):
    cfg, idx = args
    n, pr = cfg["n"], cfg["primes"]
    V = weights(cfg["H"], "signed")
    rng = random.Random(f"inversesieve|e1|n={n}|shape={cfg['shape']}|{idx}")   # same N as E1 record idx
    p, q = sample_pair(rng, n)
    N = p * q
    rng_np = np.random.default_rng(n * 7777 + idx)
    J = len(pr)
    out = []
    for j in sorted({0, J // 4, J // 2, (3 * J) // 4, J - 1}):
        l = pr[j]
        F = farey_set(l, V)
        R = forbidden_set(N, l, V)
        units = np.arange(1, l, dtype=np.int64)
        halfu = np.arange(1, (l - 1) // 2 + 1, dtype=np.int64)
        rnd = lambda k: rng_np.choice(units, size=k, replace=False)

        def rsym(k):            # random symmetric set S u (-S) of size k (k even): null for sets with A = -A
            S = rng_np.choice(halfu, size=k // 2, replace=False)
            return np.concatenate([S, l - S])
        out.append({"l": l, "nF": int(F.size), "nR": int(R.size),
                    "EF": st.energy_mod(F, l), "ER": st.energy_mod(R, l), "MF": st.mult_energy_mod(F, l),
                    "MR": st.mult_energy_mod(R, l),
                    "EF_null": [st.energy_mod(rnd(F.size), l) for _ in range(10)],
                    "ER_null": [st.energy_mod(rnd(R.size), l) for _ in range(10)],
                    "MF_null": [st.mult_energy_mod(rnd(F.size), l) for _ in range(10)],
                    "EF_symnull": [st.energy_mod(rsym(F.size), l) for _ in range(10)],
                    "ER_symnull": [st.energy_mod(rsym(R.size), l) for _ in range(10)],
                    "F_symmetric": bool(np.array_equal(np.sort((l - F) % l), F)),
                    "R_symmetric": bool(np.array_equal(np.sort((l - R) % l), R))})
    return {"n": n, "shape": cfg["shape"], "Y": cfg["Y"], "H": cfg["H"], "idx": idx, "N": N, "primes": out}


def e2(workers):
    jobs = [(c, i) for c in configs() for i in range(5)]
    run("e2", e2_one, jobs, workers)


# ---------------------------------------------------------------------------
# E3: true uncovered pairs vs control pairs
# ---------------------------------------------------------------------------
def dirichlet_quality(p, q, h):
    """min over 1 <= B <= h of |Bp - Aq| (A = nearest) times h / P (float report)."""
    P = (p * q) ** 0.5
    best = min(abs(B * p - round(B * p / q) * q) for B in range(1, h + 1))
    return best * h / P


def cf_jump(p, q, Bmax):
    """Convergents A/B of p/q: return (B, B') for the last convergent with B <= Bmax, and B'/B."""
    a, b = p, q
    hprev, h = 0, 1      # denominators
    while b:
        t = a // b
        a, b = b, a - t * b
        hprev, h = h, t * h + hprev
        if h > Bmax:
            return hprev, h
    return h, h


def pair_features(p, q, n, cfg, h0):
    N = p * q
    a0, b0 = candidate_interval(N)
    sqrtP = isqrt(isqrt(N))
    cen = st.lehman_centres(N, h0)
    P = 2 ** (n / 2)
    Xs = P ** (12 / 31)
    Ls = st.lstar_score(p, q, max(2, int(round(P ** 0.2))) * 2)
    B, B2 = cf_jump(p, q, max(1, int(Xs ** 0.5)))
    return {"lehman_in": bool(st.in_containers(np.array([p]), cen, sqrtP // 2)[0]),
            "Lstar": Ls, "Lstar_le_X": bool(Ls <= Xs), "dq": dirichlet_quality(p, q, 4 * cfg["H"]),
            "cfB": B, "cfB2": B2}


def e3_one(args):
    cfg, idx = args
    n, pr = cfg["n"], cfg["primes"]
    V = weights(cfg["H"], "signed")
    rng = random.Random(f"inversesieve|e3|n={n}|shape={cfg['shape']}|{idx}")
    plo, phi = p_range(n)
    recs = []
    if idx < 40:
        # uncovered pairs: fix p, sieve partners on half 1 then half 2
        p = random_prime(rng, plo, phi)
        lo, hi = q_range(p, n)
        h1, h2 = pr[0::2], pr[1::2]
        Q1 = uncovered_partners(p, h1, V, lo, hi)          # survivors of half 1
        Q = uncovered_partners(p, pr, V, lo, hi)           # survivors of the full band (subset of Q1)
        d2 = float(np.prod([1 - len(V) / l for l in h2]))
        qs = [int(x) for x in Q if gmpy2.is_prime(int(x))]
        held = {"Q1": int(Q1.size), "Q": int(Q.size), "pred_ratio": d2}
        rng.shuffle(qs)
        for q in qs[:5]:
            assert is_pair_uncovered(p, q, pr, V)
            N = p * q
            a0, b0 = candidate_interval(N)
            U = U_of_N(N, pr, V, (a0, b0))
            assert p in set(U.tolist())
            others = U[U != p]
            nn = np.diff(U)
            nnd = np.minimum(np.r_[np.inf, nn], np.r_[nn, np.inf])        # NN distance of every point
            i = int(np.searchsorted(U, p))
            pct = float((np.count_nonzero(nnd < nnd[i]) + 0.5 * (np.count_nonzero(nnd == nnd[i]) - 1)) / max(1, U.size - 1))
            sqrtP = isqrt(isqrt(N))
            near = int(np.count_nonzero(np.abs(others - p) <= sqrtP))
            near_exp = others.size * min(b0 - a0 + 1, 2 * sqrtP + 1) / (b0 - a0 + 1)
            r = {"kind": "uncovered", "n": n, "shape": cfg["shape"], "idx": idx, "p": p, "q": q,
                 "U_size": int(U.size), "nn_pct": pct, "near": near, "near_exp": near_exp, "held": held}
            r.update(pair_features(p, q, n, cfg, cfg["h0"]))
            recs.append(r)
        if not qs:
            recs.append({"kind": "uncovered_none", "n": n, "shape": cfg["shape"], "idx": idx, "p": p, "held": held})
    else:
        for _ in range(5):
            p, q = sample_pair(rng, n)
            r = {"kind": "control", "n": n, "shape": cfg["shape"], "idx": idx, "p": p, "q": q,
                 "uncovered": is_pair_uncovered(p, q, pr, V)}
            r.update(pair_features(p, q, n, cfg, cfg["h0"]))
            recs.append(r)
    return recs


def e3(workers):
    cf = [_with_h0(c) for c in configs()]
    jobs = [(c, i) for c in cf for i in range(80)]
    path = os.path.join(RAW, "inversesieve_e3.jsonl")
    if os.path.exists(path):
        os.remove(path)
    t0 = time.time()
    with Pool(workers) as pool, open(path, "w") as f:
        for i, recs in enumerate(pool.imap(e3_one, jobs, chunksize=1)):
            for rec in recs:
                f.write(json.dumps(rec) + "\n")
            f.flush()
            if (i + 1) % 100 == 0:
                print(f"e3: {i + 1}/{len(jobs)} ({time.time() - t0:.0f}s)", flush=True)
    print(f"e3 done in {time.time() - t0:.0f}s", flush=True)


# ---------------------------------------------------------------------------
# E4: dchworst e1 worst pairs
# ---------------------------------------------------------------------------
def e4_one(rec):
    n, Y, signed, J, hstar = rec["n"], rec["Y"], rec["signed"], rec["J"], rec["hstar"]
    p, q = rec["worst_pair"]
    kind = "signed_dw" if signed else "positive"
    H = hstar - 1
    pr = band(Y)[:J]
    V = weights(H, kind)
    unc = is_pair_uncovered(p, q, pr, V)
    N = p * q
    a0, b0 = candidate_interval(N)
    U = U_of_N(N, pr, V, (a0, b0))
    inU = bool(p in set(U.tolist()))
    E_model = (b0 - a0 + 1) * float(np.prod([1 - len(forbidden_set(N, l, V)) / l for l in pr]))
    nn = np.diff(U)
    nnd = np.minimum(np.r_[np.inf, nn], np.r_[nn, np.inf]) if U.size > 1 else np.array([np.inf])
    i = int(np.searchsorted(U, p))
    pct = float((np.count_nonzero(nnd < nnd[i]) + 0.5 * (np.count_nonzero(nnd == nnd[i]) - 1)) / max(1, U.size - 1)) if inU and U.size > 1 else None
    P = 2 ** (n / 2)
    out = {"n": n, "Y": Y, "signed": signed, "J": J, "hstar": hstar, "H": H, "p": p, "q": q,
           "uncovered": unc, "inU": inU, "U_size": int(U.size), "E_model": E_model, "nn_pct": pct,
           "p_eq_q_mod_l1": bool((p - q) % pr[0] == 0),
           "Lstar": st.lstar_score(p, q, max(2, int(round(P ** 0.2))) * 2),
           "dq": dirichlet_quality(p, q, max(4, 4 * H)), "lehman_in": None}
    cen = st.lehman_centres(N, lehman_h0(n))
    out["lehman_in"] = bool(st.in_containers(np.array([p]), cen, isqrt(isqrt(N)) // 2)[0])
    return out


def e4_controls(n, count=2000):
    rng = random.Random(f"inversesieve|e4ctl|n={n}")
    h0 = lehman_h0(n)
    P = 2 ** (n / 2)
    out = []
    for _ in range(count):
        p, q = sample_pair(rng, n)
        N = p * q
        cen = st.lehman_centres(N, h0)
        out.append({"kind": "e4control", "n": n, "p": p, "q": q,
                    "Lstar": st.lstar_score(p, q, max(2, int(round(P ** 0.2))) * 2),
                    "lehman_in": bool(st.in_containers(np.array([p]), cen, isqrt(isqrt(N)) // 2)[0])})
    return out


def e4(workers):
    src = [json.loads(l) for l in open(os.path.join(RAW, "dchworst_e1.jsonl"))]
    path = os.path.join(RAW, "inversesieve_e4.jsonl")
    if os.path.exists(path):
        os.remove(path)
    with Pool(workers) as pool, open(path, "w") as f:
        for rec in pool.imap(e4_one, src, chunksize=1):
            f.write(json.dumps(rec) + "\n")
        for ctl in pool.map(e4_controls, (32, 36)):
            for r in ctl:
                f.write(json.dumps(r) + "\n")
    print("e4 done", flush=True)


if __name__ == "__main__":
    what = sys.argv[1]
    workers = int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 8
    os.makedirs(RAW, exist_ok=True)
    for name, fn in (("e1", e1), ("e2", e2), ("e3", e3), ("e4", e4)):
        if what in (name, "all"):
            fn(workers)
