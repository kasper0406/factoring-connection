"""Pre-registered experiments of agent below16 (notes/below16_prereg.md).

  PYTHONPATH=. .venv/bin/python -m src.below16.experiments e1     # T1/T2 bound checks (exhaustive S_n)
  PYTHONPATH=. .venv/bin/python -m src.below16.experiments e2     # multi-base reduction checks
  PYTHONPATH=. .venv/bin/python -m src.below16.experiments e4     # deterministic QS statistics
  PYTHONPATH=. .venv/bin/python -m src.below16.experiments summarize
"""
from __future__ import annotations

import json
import math
import random
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import gmpy2
import numpy as np

from . import avgbound as AB
from . import multibase as MB
from . import qs_pilot as QS

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "experiments" / "raw"
CONFIG = {
    "e1": {"ns": [20, 22, 24, 26, 28, 30], "seed": 20261018, "n_random": 150, "n_oracle": 40,
           "t2_ns": [20, 22, 24, 26], "n_t2": 8, "dcl_H": 3},
    "e2": {"n": 22, "seed": 20261019, "trials": 60,
           "basesets": [[2, 3, 6], [2, 4, 8], [2, 3, 12, 18], [6, 10, 15], [2, 9, 12]]},
    "e4": {"ns": [40, 48, 56, 64, 72, 80], "seed": 20261020, "blum": 400, "general": 200,
           "blum_80": 200, "general_80": 100},
}


# --------------------------------------------------------------------------
# E1
# --------------------------------------------------------------------------

def _rand_window(rng, center):
    t = rng.choice(["interval", "ap", "gap", "rand"])
    if t == "interval":
        W = rng.randint(1, 64)
        U = center - rng.randint(0, W)
        return t, frozenset(range(U, U + W))
    if t == "ap":
        d, s = rng.randint(2, 9), rng.randint(1, 32)
        U = center - rng.randint(0, d * s)
        return t, frozenset(U + d * i for i in range(s))
    if t == "gap":
        d1, d2 = rng.randint(1, 5), rng.randint(6, 40)
        s1, s2 = rng.randint(1, 6), rng.randint(1, 6)
        U = center - rng.randint(0, d1 * s1 + d2 * s2)
        return t, frozenset(U + d1 * i + d2 * j for i in range(s1) for j in range(s2))
    span, s = rng.randint(8, 200), rng.randint(1, 40)
    U = center - rng.randint(0, span)
    return t, frozenset(rng.sample(range(U, U + span), min(s, span)))


def _small_primes(lo, hi):
    return [int(l) for l in AB.primes_upto(hi) if l >= lo]


def random_t1(rng, n, p, q):
    xmin, _, _ = AB.class_var_range(n, "p")
    Mmax = max(3, min(40, math.isqrt(xmin - 1)))
    M = rng.randint(1, Mmax)
    divs = [t for t in range(1, M + 1) if M % t == 0]
    K = M * rng.choice(divs)
    if K >= xmin:
        K = M
    D = rng.randint(1, max(1, M * M))
    A = 0 if rng.random() < 0.15 else rng.choice([a for a in range(-6, 7) if a != 0])
    B = rng.choice([b for b in range(-6, 7) if b != 0])
    side = rng.choice(["p", "q"])
    x, y = (p, q) if side == "p" else (q, p)
    vals = (B * x + A * y) // D
    R, S, kinds = {}, {}, Counter()
    for kap in range(K):
        if math.gcd(kap, K) != 1 or rng.random() < 0.5:
            continue
        R[kap] = frozenset(rng.sample(range(M), min(M, rng.randint(1, 3))))
        center = int(vals[rng.randrange(len(vals))])
        kind, S[kap] = _rand_window(rng, center)
        kinds[kind] += 1
    g = AB.GenT1(A, B, M, K, D, side, R, S)
    return g, dict(kinds)


def oracle_t1(rng, n, p, q, s):
    """Per key: the class and the s values maximising coverage (stress test)."""
    xmin, _, _ = AB.class_var_range(n, "p")
    pr = _small_primes(3, max(5, min(60, math.isqrt(xmin - 1))))
    M = rng.choice(pr)
    K = M * (M if rng.random() < 0.5 and M * M < xmin else 1)
    D = rng.choice([M, M * M, max(1, M * M // 2)])
    A = rng.choice([1, 2, 3, -1, -2])
    B = rng.choice([1, 2, 3, 5])
    side = rng.choice(["p", "q"])
    x, y = (p, q) if side == "p" else (q, p)
    kap = ((p * q) % K).tolist()
    cls = (x % M).tolist()
    val = ((B * x + A * y) // D).tolist()
    groups = defaultdict(Counter)
    for k, c, v in zip(kap, cls, val):
        groups[(k, c)][v] += 1
    best = {}
    for (k, c), cnt in groups.items():
        tot = sum(m for _, m in cnt.most_common(s))
        if k not in best or tot > best[k][0]:
            best[k] = (tot, c, frozenset(v for v, _ in cnt.most_common(s)))
    R = {k: frozenset({c}) for k, (_, c, _) in best.items()}
    S = {k: ss for k, (_, _, ss) in best.items()}
    return AB.GenT1(A, B, M, K, D, side, R, S)


def oracle_t2(rng, n, p, q, W):
    pr = _small_primes(3, 30)
    M = rng.choice(pr)
    K = M
    D = rng.choice([M, M * M])
    A, B = rng.choice([1, 2, 3]), rng.choice([1, 2, 3, 4])
    Delta = AB.t2_delta(n, K)
    g = AB.GenT2(A, B, M, K, D, Delta)
    kap, b = AB.keys_T2(g, n, p, q)
    val = (B * p + A * q) // D
    cls = p % M
    groups = defaultdict(list)
    for k, bb, c, v in zip(kap.tolist(), b.tolist(), cls.tolist(), val.tolist()):
        groups[(k, bb, c)].append(v)
    best = {}
    for (k, bb, c), vs in groups.items():
        vs.sort()
        j, top, arg = 0, 0, vs[0]
        for i in range(len(vs)):
            while vs[i] - vs[j] >= W:
                j += 1
            if i - j + 1 > top:
                top, arg = i - j + 1, vs[j]
        key = (k, bb)
        if key not in best or top > best[key][0]:
            best[key] = (top, c, arg)
    g.R = {key: frozenset({c}) for key, (_, c, _) in best.items()}
    g.UW = {key: (u, W) for key, (_, _, u) in best.items()}
    return g


def run_e1():
    cfg = CONFIG["e1"]
    rng = random.Random(cfg["seed"])
    out = open(RAW / "below16_e1.jsonl", "w")
    for n in cfg["ns"]:
        t0 = time.time()
        p, q = AB.population(n)
        Sn = len(p)
        P = 2.0 ** (n / 2)
        recs = []
        # (a) DCL
        H = cfg["dcl_H"]
        Y = math.ceil(2 ** (n / 6) / H ** (1 / 3))
        for ksq in (False, True):
            for g in AB.dcl_generators(n, H, Y, K_sq=ksq):
                if g.K >= AB.class_var_range(n, "p")[0]:
                    continue
                recs.append(("dcl_K" + ("l2" if ksq else "l"), g, {}))
        # (b) random
        for _ in range(cfg["n_random"]):
            g, kinds = random_t1(rng, n, p, q)
            recs.append(("random", g, kinds))
        # (c) oracle-local
        for i in range(cfg["n_oracle"]):
            s = [1, 4, 16, 64][i % 4]
            recs.append((f"oracle_s{s}", oracle_t1(rng, n, p, q, s), {}))
        for fam, g, kinds in recs:
            c = AB.count_T1(g, p, q)
            core = AB.bound_T1_core(g, n)
            cf = AB.bound_T1_closed(g, n)
            cond = AB.t1_conditions(g, n)
            out.write(json.dumps({"exp": "e1_t1", "n": n, "family": fam, "A": g.A, "B": g.B,
                                  "M": g.M, "K": g.K, "D": g.D, "side": g.side, "kinds": kinds,
                                  "count": c, "core": core, "closed": cf, "cond": cond,
                                  "frac": c / Sn, "Sn": Sn}) + "\n")
        # (d) T2
        if n in cfg["t2_ns"]:
            for i in range(cfg["n_t2"]):
                W = [1, 4, 16, 64][i % 4]
                g = oracle_t2(rng, n, p, q, W)
                c = AB.count_T2(g, n, p, q)
                core = AB.bound_T2_core(g, n)
                cf = AB.bound_T2_closed(g, n)
                out.write(json.dumps({"exp": "e1_t2", "n": n, "A": g.A, "B": g.B, "M": g.M,
                                      "K": g.K, "D": g.D, "W": W, "Delta": g.Delta,
                                      "count": c, "core": core, "closed": cf,
                                      "frac": c / Sn, "Sn": Sn}) + "\n")
        out.flush()
        print(n, "done", round(time.time() - t0, 1), "s", flush=True)
    out.close()


# --------------------------------------------------------------------------
# E2
# --------------------------------------------------------------------------

def run_e2():
    cfg = CONFIG["e2"]
    rng = random.Random(cfg["seed"])
    n = cfg["n"]
    p, q = AB.population(n)
    out = open(RAW / "below16_e2.jsonl", "w")
    for bases in cfg["basesets"]:
        primes, E = MB.prime_exponents(bases)
        k = len(bases)
        for t in range(cfg["trials"]):
            forms = [(rng.randint(-3, 3), rng.randint(-3, 3)) for _ in range(k)]
            W = [rng.randint(1, 6) for _ in range(k)]
            i0 = rng.randrange(len(p))
            c = [b * int(p[i0]) + a * int(q[i0]) - rng.randint(0, w - 1)
                 for (b, a), w in zip(forms, W)]
            cov = MB.covered_multibase(p, q, forms, c, W, E)
            rec = {"exp": "e2", "bases": bases, "forms": forms, "W": W, "n_cov": int(cov.sum()),
                   "box": int(np.prod(W)), "img": len(MB.image_of_box(E, W)), "coords": []}
            for j in range(len(primes)):
                cj, sj = MB.covered_coordinate(p, q, forms, c, W, E, j)
                nonconst = any(int(E[i, j]) * (abs(forms[i][0]) + abs(forms[i][1])) for i in range(k))
                rec["coords"].append({"prime": primes[j], "n_cov": int(cj.sum()), "Sj": sj,
                                      "superset": bool(np.all(cj | ~cov)), "nonconst": bool(nonconst)})
            out.write(json.dumps(rec) + "\n")
    out.close()


# --------------------------------------------------------------------------
# E4
# --------------------------------------------------------------------------

def sample_pair(rng, n, blum):
    lo, hi = 1 << (n - 1), 1 << n
    while True:
        p = rng.randrange(math.isqrt(lo // 2), math.isqrt(hi))
        q = rng.randrange(p + 1, 2 * p + 1)
        if not (lo <= p * q < hi):
            continue
        if blum and (p % 4 != 3 or q % 4 != 3):
            continue
        if gmpy2.is_prime(p) and gmpy2.is_prime(q):
            return p, q


def run_e4():
    cfg = CONFIG["e4"]
    out = open(RAW / "below16_e4.jsonl", "w")
    for n in cfg["ns"]:
        rng = random.Random(f"{cfg['seed']}|{n}")
        nb = cfg["blum_80"] if n == 80 else cfg["blum"]
        ng = cfg["general_80"] if n == 80 else cfg["general"]
        t0 = time.time()
        for blum, cnt in ((True, nb), (False, ng)):
            for i in range(cnt):
                p, q = sample_pair(rng, n, blum)
                N = p * q
                r = QS.qs_factor(N, QS.default_config(n))       # sees only N
                f = r["factor"]
                ok = f in (p, q)
                deps = r.get("deps", [])
                rule = None
                if blum and deps:
                    rule = sum((d["nontriv"] == (d["jacobi"] == -1)) for d in deps)
                prefix = []
                for j in (1, 2, 3, 4, 5, 6, 8, 10):
                    if len(deps) >= j:
                        prefix.append([j, all(not d["nontriv"] for d in deps[:j])])
                out.write(json.dumps({"exp": "e4", "n": n, "blum": blum, "i": i,
                                      "ok": ok, "how": r["how"], "fb": r.get("fb_size"),
                                      "rels": r.get("n_rels"), "kdim": r.get("kdim"),
                                      "X0": r.get("X0_final"),
                                      "n_nontriv": sum(d["nontriv"] for d in deps),
                                      "n_deps": len(deps), "jacobi_rule_agree": rule,
                                      "prefix_all_trivial": prefix,
                                      "class": ("STRUCTURAL" if r["how"] == "dependency" else
                                                "ACCIDENTAL" if r["how"] == "trial" else str(r["how"]))})
                          + "\n")
            out.flush()
        print(n, "done", round(time.time() - t0, 1), "s", flush=True)
    out.close()


# --------------------------------------------------------------------------

def summarize():
    res = {}
    e1 = [json.loads(l) for l in open(RAW / "below16_e1.jsonl")]
    viol_core = [r for r in e1 if r["count"] > r["core"]]
    viol_cf = [r for r in e1 if r.get("cond", True) and r["core"] > r["closed"] + 1e-6]
    viol_cf2 = [r for r in e1 if r["count"] > r["closed"] + 1e-6 and r.get("cond", True)]
    fam = defaultdict(list)
    for r in e1:
        fam[(r["exp"], r.get("family", "t2"), r["n"])].append(r)
    res["e1"] = {"n_records": len(e1), "count>core": len(viol_core), "core>closed(cond)": len(viol_cf),
                 "count>closed(cond)": len(viol_cf2),
                 "n_cond_false": sum(1 for r in e1 if r.get("cond") is False),
                 "families": {f"{k[0]}|{k[1]}|{k[2]}": {
                     "k": len(v), "max_count/core": max(r["count"] / max(1, r["core"]) for r in v),
                     "max_count/closed": max(r["count"] / max(1e-9, r["closed"]) for r in v),
                     "median_count/closed": float(np.median([r["count"] / max(1e-9, r["closed"]) for r in v]))}
                     for k, v in sorted(fam.items())}}
    try:
        e2 = [json.loads(l) for l in open(RAW / "below16_e2.jsonl")]
        bad = sum(1 for r in e2 for c in r["coords"] if not c["superset"])
        res["e2"] = {"targets": len(e2), "superset_failures": bad,
                     "img<=box": all(r["img"] <= r["box"] for r in e2),
                     "Sj<=img": all(c["Sj"] <= r["img"] for r in e2 for c in r["coords"])}
    except FileNotFoundError:
        pass
    try:
        e4 = [json.loads(l) for l in open(RAW / "below16_e4.jsonl")]
        by = defaultdict(list)
        for r in e4:
            by[(r["n"], r["blum"])].append(r)
        e4s = {}
        for (n, bl), v in sorted(by.items()):
            deps = sum(r["n_deps"] for r in v)
            nt = sum(r["n_nontriv"] for r in v)
            pref = {}
            for j in (1, 2, 3, 4, 5, 6, 8, 10):
                xs = [x[1] for r in v for x in r["prefix_all_trivial"] if x[0] == j]
                if xs:
                    pref[j] = {"m": len(xs), "all_trivial": sum(xs), "expected": len(xs) * 2.0 ** (-j)}
            sm = [r["rels"] for r in v if r["rels"] is not None]
            e4s[f"{n}|{'blum' if bl else 'general'}"] = {
                "inst": len(v), "success": sum(r["ok"] for r in v) / len(v),
                "classes": dict(Counter(r["class"] for r in v)),
                "deps": deps, "nontriv_rate": nt / max(1, deps),
                "jacobi_rule_agree": (sum(r["jacobi_rule_agree"] or 0 for r in v), deps) if bl else None,
                "kdim_min": min((r["kdim"] or 0) for r in v),
                "prefix": pref,
                "rels_mean": float(np.mean(sm)) if sm else None}
        res["e4"] = e4s
    except FileNotFoundError:
        pass
    p = ROOT / "experiments" / "analysis" / "below16_summary.json"
    p.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1)[:6000])


if __name__ == "__main__":
    cmd = sys.argv[1]
    {"e1": run_e1, "e2": run_e2, "e4": run_e4, "summarize": summarize}[cmd]()
