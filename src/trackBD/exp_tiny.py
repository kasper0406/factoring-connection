"""Experiment B1: exhaustive sketch inspection on all composites 4 <= N < LIMIT.

HARNESS script (knows factors via trial division; operational readouts are
computed by src.trackBD.sketch from N and f only).

For every composite N and every modulus f in the public family
  binom : X^r - c,  r = 1..R, c in C
  cyclo : Phi_m,    m = 1..M
  rand  : 20 fixed random monic squarefree f of degree 2..8 (seeded)
compute S = F_N mod (f, N) and the readouts (coefficient gcds, content gcd,
Euclid/resultant with nonunit leading coefficients), then classify every
prime discovered (classes in src/trackBD/diagnose.py).  Regime tag:
'cmp' if deg f < p (compressed), 'unc' otherwise (p = least prime of N).

Usage:
  .venv/bin/python -m src.trackBD.exp_tiny --limit 5000 --R 200 --M 100 --procs 24 \
      --out experiments/raw/trackBD_tiny.jsonl
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from collections import Counter
from multiprocessing import Pool

import flint
import numpy as np

from src.trackBD.sketch import binomial_modulus, coef_gcds, content_gcd, euclid_readout, sketch
from src.trackBD import diagnose as D

C_SET = (1, -1, 2, 3, 5)


def factor_small(N):
    fs, M, d = [], N, 2
    while d * d <= M:
        while M % d == 0:
            fs.append(d)
            M //= d
        d += 1
    if M > 1:
        fs.append(M)
    return fs


def random_family(count=20, seed=12345):
    rng = random.Random(f"trackBD_randf|{seed}")
    out = []
    while len(out) < count:
        d = rng.randint(2, 8)
        f = [rng.randint(-5, 5) for _ in range(d)] + [1]
        P = flint.fmpz_poly(f)
        if P.discriminant() == 0 or f[0] == 0:
            continue
        out.append(f)
    return out


def families(R, M):
    fam = []
    for r in range(1, R + 1):
        for c in C_SET:
            fam.append((("binom", r, c), binomial_modulus(r, c)))
    for m in range(1, M + 1):
        cs = [int(x) for x in flint.fmpz_poly.cyclotomic(m).coeffs()]
        fam.append((("cyclo", m), cs))
    for i, f in enumerate(random_family()):
        fam.append((("rand", i), f))
    return fam


_FAM = None
_PUB = None
_COMP_CACHE = {}


def _init(R, M):
    global _FAM, _PUB
    _FAM = families(R, M)
    _PUB = [D.public_integers(f) for _, f in _FAM]


def _comps(p, idx, f):
    key = (p, idx)
    v = _COMP_CACHE.get(key)
    if v is None:
        v = D.components(p, f)
        _COMP_CACHE[key] = v
    return v


def lucas_valid_t(p, m):
    """t in [1, m-1] with binom(m, t) != 0 mod p (Lucas)."""
    digs = []
    M = m
    while M:
        digs.append(M % p)
        M //= p
    ts = np.arange(1, m, dtype=np.int64)
    ok = np.ones(len(ts), dtype=bool)
    T = ts.copy()
    for dg in digs:
        ok &= (T % p) <= dg
        T //= p
    return ts[ok]


PREC = ["PRIMEPOW", "TD", "FULL", "SPLIT", "SUPPORT", "ACC"]


def _rank(lab):
    for i, k in enumerate(PREC):
        if lab.startswith(k):
            return i
    return len(PREC)


def _labs(comps):
    s = {c["label"] for c in comps if c["root_v"] and c["label"] not in (None, "TRIVIAL")}
    return "+".join(sorted(s)) if s else "TRIVIAL"


def analyse_N(N):
    t0 = time.perf_counter()
    factors = factor_small(N)
    pd = D.prime_data(N, factors)
    pmin = pd[0][0]
    primepow = len(pd) == 1
    agg = Counter()
    events = []
    valid_t = {p: (lucas_valid_t(p, m) * (p ** k)) for p, k, m in pd}
    for idx, (kind, f) in enumerate(_FAM):
        deg = len(f) - 1
        regime = "cmp" if deg < pmin else "unc"
        fam = kind[0]
        pre = f"{fam}|{regime}"
        S = sketch(N, f)
        info = {}
        for p, k, m in pd:
            Sp = flint.nmod_poly([x % p for x in S], p)
            comps = []
            for g, e in _comps(p, idx, f):
                root_v = (Sp % g).is_zero()
                # Theorem B1 check: g | S mod p  <=>  (1+a)^m == 1 + a^m in F_p[X]/(g)
                if root_v != D.H_vanishes(p, g, m):
                    agg[f"{pre}|THEOREM_VIOLATION"] += 1
                lab = None
                if root_v:
                    if m == 1:
                        lab = "PRIMEPOW"
                    elif D.is_universal_root_factor(p, g):
                        lab = "TRIVIAL"
                    elif e > 1:
                        lab = "RAMIFIED"
                    else:
                        lab = "FROB" if D.frob_check(p, g, m)["frob"] else "SPORADIC"
                    agg[f"{pre}|vanish_root|{lab}|d{min(g.degree(), 3)}"] += 1
                comps.append({"deg": g.degree(), "e": e, "root_v": root_v, "label": lab})
            info[p] = {"k": k, "m": m, "full": Sp.is_zero(),
                       "deg_gcd": sum(c["deg"] for c in comps if c["root_v"]),
                       "comps": comps, "td": D.td_divides(p, _PUB[idx])}
        any_td = any(info[p]["td"] for p in info)
        agg[f"{pre}|sketches"] += 1
        agg[f"{pre}|coeffs"] += deg

        def prime_full_label(p):
            inf = info[p]
            if primepow:
                return "PRIMEPOW"
            if inf["td"]:
                return "TD:" + "+".join(sorted(inf["td"]))
            return "FULL_" + _labs(inf["comps"])

        # readout 1: coefficient gcds ------------------------------------
        hits = coef_gcds(S, N)
        best = None
        for j, d in hits:
            for p, k, m in pd:
                if d % p:
                    continue
                inf = info[p]
                if primepow or inf["td"] or inf["full"]:
                    lab = prime_full_label(p)
                elif kind[0] == "binom" and not np.any((valid_t[p] % kind[1]) == j):
                    lab = "SUPPORT"
                else:
                    lab = "ACC"
                agg[f"{pre}|coefhit|{lab}"] += 1
                if best is None or _rank(lab) < _rank(best):
                    best = lab
        if best is not None:
            agg[f"{pre}|coef|{best}"] += 1
            if _rank(best) >= 2 and best != "ACC":
                events.append({"N": N, "f": list(kind), "readout": "coef", "class": best, "regime": regime})
        # readout 2: content gcd -----------------------------------------
        ct = content_gcd(S, N)
        if 1 < ct < N:
            labs = [prime_full_label(p) for p, k, m in pd if ct % p == 0]
            lab = min(labs, key=_rank)
            agg[f"{pre}|content|{lab}"] += 1
            if _rank(lab) >= 2:
                events.append({"N": N, "f": list(kind), "readout": "content", "class": lab,
                               "d": ct, "regime": regime})
        # readout 3: Euclid / resultant ----------------------------------
        eu = euclid_readout(f, S, N)
        if eu.get("factor"):
            degs = {p: info[p]["deg_gcd"] for p in info}
            if primepow:
                lab = "PRIMEPOW"
            elif any_td:
                lab = "TD:" + "+".join(sorted({t for p in info for t in info[p]["td"]}))
            elif len(set(degs.values())) > 1:
                lo = min(degs.values())
                extra = set()
                for p in info:
                    if degs[p] > lo:
                        extra |= {c["label"] for c in info[p]["comps"] if c["root_v"]}
                extra.discard("TRIVIAL")
                lab = "SPLIT_" + ("+".join(sorted(extra)) if extra else "TRIVIAL")
            else:
                lab = "ACC"
            agg[f"{pre}|euclid|{lab}"] += 1
            if lab.startswith("SPLIT"):
                events.append({"N": N, "f": list(kind), "readout": "euclid", "class": lab,
                               "d": eu["factor"], "regime": regime,
                               "deg_gcd": {str(p): v for p, v in degs.items()}})
    return {"N": N, "factors": factors, "agg": dict(agg), "events": events,
            "wall_s": round(time.perf_counter() - t0, 3)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=5000)
    ap.add_argument("--start", type=int, default=4)
    ap.add_argument("--R", type=int, default=200)
    ap.add_argument("--M", type=int, default=100)
    ap.add_argument("--procs", type=int, default=24)
    ap.add_argument("--out", default="experiments/raw/trackBD_tiny.jsonl")
    a = ap.parse_args()
    Ns = [N for N in range(a.start, a.limit) if len(factor_small(N)) > 1]
    t0 = time.time()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with Pool(a.procs, initializer=_init, initargs=(a.R, a.M)) as pool, open(a.out, "w") as fh:
        for i, rec in enumerate(pool.imap_unordered(analyse_N, Ns, chunksize=4)):
            fh.write(json.dumps(rec) + "\n")
            if i % 500 == 0:
                print(f"{i}/{len(Ns)} {time.time() - t0:.0f}s", file=sys.stderr, flush=True)
    print(f"done {len(Ns)} composites in {time.time() - t0:.0f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
