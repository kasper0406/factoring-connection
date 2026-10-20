"""Experiments E1-E3 of notes/floorhard_prereg.md. HARNESS (brute-force oracles; costs are model calls).

Usage: PYTHONPATH=. .venv/bin/python -m src.floorhard.run_experiments {e1,e2,e3,summary}
Raw JSONL: experiments/raw/floorhard_e{1,2,3}.jsonl; summary: experiments/analysis/floorhard_summary.json
"""
from __future__ import annotations

import json
import os
import random
import sys
import time
from math import gcd

import gmpy2

from src.common.instances import random_prime
from src.nonsmooth import forms
from src.floorhard import core

RAW = "experiments/raw"
ANA = "experiments/analysis"


def _dump(name, rows):
    os.makedirs(RAW, exist_ok=True)
    with open(f"{RAW}/floorhard_{name}.jsonl", "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")


def _semiprime(rng, bits, cls):
    """cls 'm3': pq = 3 mod 4; 'pq33': p = q = 3 mod 4. Balanced, p != q, both > 3."""
    while True:
        p = random_prime(rng, bits // 2)
        q = random_prime(rng, bits - bits // 2)
        if p == q or min(p, q) <= 3 or (p * q) % 3 == 0:
            continue
        m = p * q
        if cls == "m3" and m % 4 == 3:
            return p, q
        if cls == "pq33" and p % 4 == 3 and q % 4 == 3:
            return p, q


def e1():
    rows = []
    # E1a
    fails = checked = 0
    for m in range(3, 3000, 2):
        if not core.squarefree_odd(m):
            continue
        for a in (1, 2, 3, 5, 7, 11):
            if gcd(a, m) != 1:
                continue
            checked += 1
            if core.full_period_brute(m, a) != core.full_period_rhs(m, a):
                fails += 1
                rows.append(dict(part="E1a_fail", m=m, a=a))
    rows.append(dict(part="E1a", checked=checked, fails=fails))
    print("E1a", checked, fails, flush=True)
    # E1b
    rng = random.Random(20261020)
    for i in range(30):
        cls = "m3" if i < 15 else "pq33"
        bits = 16 + (i % 15) * 6 // 14        # 16..22
        p, q = _semiprime(rng, bits, cls)
        m = p * q
        t0 = time.time()
        M, hM, calls = core.decode_class_number_semiprime(m, rng, amax=1000)  # as run (prereg); theorem samples mod M
        true_h = forms.class_number(-M)
        g = core.factor_from_class_number(m, M, hM, rng) if hM is not None else None
        r = dict(part="E1b", cls=cls, bits=m.bit_length(), m=m, M=M, h_decoded=hM, h_true=true_h,
                 h_ok=(hM == true_h), calls=calls, factor=g,
                 factored=(g is not None and m % g == 0 and 1 < g < m), secs=round(time.time() - t0, 2))
        rows.append(r)
        print(r, flush=True)
    _dump("e1", rows)


def e2():
    rows = []
    rng = random.Random(777)
    for i in range(40):
        w = [rng.randint(1, 6) for _ in range(3)]
        T = sum(x for x in w if rng.random() < 0.5) if i % 2 == 0 else rng.randint(0, sum(w))
        inst = core.subset_sum_to_qc(w, T)
        truth = core.subset_sum_count(w, T)
        fs = core.qc_count_floorsums(inst, oracle=core.floor_sum_np)
        rows.append(dict(part="E2a", w=w, T=T, b=inst["b"], truth=truth, floorsum_count=fs, ok=(fs == truth)))
    print("E2a", sum(r["ok"] for r in rows), "/", len(rows), flush=True)
    for i in range(200):
        k = 4 + i % 7
        w = [rng.randint(1, 2 ** k) for _ in range(k)]
        T = sum(x for x in w if rng.random() < 0.5) if i % 2 == 0 else rng.randint(0, sum(w))
        inst = core.subset_sum_to_qc(w, T)
        truth = core.subset_sum_count(w, T)
        c = core.qc_count_crt(inst)
        rows.append(dict(part="E2b", k=k, w=w, T=T, bbits=(inst["b"].bit_length() if inst else 0), truth=truth, crt_count=c,
                         ok=(c == truth)))
    print("E2b", sum(r["ok"] for r in rows if r["part"] == "E2b"), "/ 200", flush=True)
    _dump("e2", rows)


def e3():
    rows = []
    rng = random.Random(4242)
    for i in range(60):
        j = 2 + i % 3
        b = rng.randint(2, 200)
        A = rng.randint(1, 5)
        cp = rng.randint(1, 3 * b)
        emb = core.TaylorEmbedding(A, b, cp, j)
        orc = core.PrefixOracle(cp)
        ok = emb.H >= cp
        for r in range(b):
            want = sum(1 for h in range(cp) if (A * h * h - r) % b == 0)
            if emb.window_count(orc, r, 0, cp) != want:
                ok = False
                break
        rows.append(dict(part="E3a", j=j, A=A, b=b, cp=cp, t=emb.t, Nbits=emb.N.bit_length(), H=emb.H, ok=ok))
    print("E3a", sum(r["ok"] for r in rows), "/ 60", flush=True)
    for i in range(20):
        b = rng.randint(3, 120)
        A = rng.randint(1, 5)
        cp = rng.randint(1, 2 * b)
        emb = core.find_prime_embedding(A, b, cp, 2, rng)
        ok = emb is not None and gmpy2.is_prime(emb.N) and emb.H >= cp
        if ok:
            orc = core.PrefixOracle(cp)
            for r in range(b):
                want = sum(1 for h in range(cp) if (A * h * h - r) % b == 0)
                if emb.window_count(orc, r, 0, cp) != want:
                    ok = False
                    break
        rows.append(dict(part="E3b", A=A, b=b, cp=cp, e=(emb.e if emb else None),
                         Nbits=(emb.N.bit_length() if emb else None), ok=ok))
    print("E3b", sum(r["ok"] for r in rows if r["part"] == "E3b"), "/ 20", flush=True)
    done = 0
    while done < 40:
        m = rng.randrange(15, 20000) | 1
        f = [p for p in range(3, 200, 2) if m % p == 0]
        # need >= 2 distinct primes; reject prime powers / primes (checked with sympy-free test)
        if gmpy2.is_prime(m):
            continue
        if any(gmpy2.iroot(m, k)[1] for k in range(2, 15)):
            # perfect powers: bisection precondition may fail; skip (prereg: >= 2 distinct primes)
            # (p^k has only 2 roots of 1)
            continue
        # m with a single prime factor and exponent: excluded by perfect-power test above
        g, calls, N = core.factor_via_taylor_class(m, 2)
        ok = g is not None and 1 < g < m and m % g == 0
        bound = 2 * (m - 1).bit_length() + 4
        rows.append(dict(part="E3c", m=m, factor=g, calls=calls, bound=bound, Nbits=N.bit_length(),
                         ok=ok and calls <= bound))
        done += 1
    print("E3c", sum(r["ok"] for r in rows if r["part"] == "E3c"), "/ 40", flush=True)
    _dump("e3", rows)


def e2c():
    """POST-HOC (after referee round): deterministic prime-power variant and out-of-range targets."""
    rows = []
    rng = random.Random(9090)
    for i in range(300):
        k = 1 + i % 10
        w = [rng.randint(1, 2 ** min(k, 8)) for _ in range(k)]
        T = sum(x for x in w if rng.random() < 0.5) if i % 3 else rng.randint(-5, sum(w) + 5)
        inst = core.subset_sum_to_qc(w, T, deterministic=True)
        truth = core.subset_sum_count(w, T) if T >= 0 else 0
        c = core.qc_count_crt(inst)
        rows.append(dict(part="E2c", k=k, w=w, T=T, truth=truth, crt_count=c, ok=(c == truth)))
    print("E2c", sum(r["ok"] for r in rows), "/", len(rows), flush=True)
    _dump("e2c", rows)


def summary():
    out = {}
    for e in ("e1", "e2", "e2c", "e3"):
        path = f"{RAW}/floorhard_{e}.jsonl"
        if not os.path.exists(path):
            continue
        rows = [json.loads(l) for l in open(path)]
        by = {}
        for r in rows:
            by.setdefault(r["part"], []).append(r)
        for part, rs in by.items():
            if part == "E1a":
                out[part] = rs[0]
            elif part == "E1b":
                calls = sorted(r["calls"] for r in rs)
                out[part] = dict(n=len(rs), h_ok=sum(r["h_ok"] for r in rs),
                                 factored=sum(r["factored"] for r in rs),
                                 median_calls=calls[len(calls) // 2], max_calls=calls[-1])
            elif part == "E1a_fail":
                out[part] = len(rs)
            else:
                d = dict(n=len(rs), ok=sum(bool(r["ok"]) for r in rs))
                if part == "E3c":
                    d["max_calls"] = max(r["calls"] for r in rs)
                    d["max_calls_over_log2m"] = max(r["calls"] / (r["m"].bit_length()) for r in rs)
                if part == "E2a":
                    d["with_solutions"] = sum(r["truth"] > 0 for r in rs)
                if part == "E2b":
                    d["with_solutions"] = sum(r["truth"] > 0 for r in rs)
                    d["max_bbits"] = max(r["bbits"] for r in rs)
                out[part] = d
    os.makedirs(ANA, exist_ok=True)
    json.dump(out, open(f"{ANA}/floorhard_summary.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    {"e1": e1, "e2": e2, "e2c": e2c, "e3": e3, "summary": summary}[sys.argv[1]]()
