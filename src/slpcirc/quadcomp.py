"""E4: fibres of compositions of x^2 + c_i (c_i in {-2..2}) mod p, by backward preimage propagation.

Usage: python -m src.slpcirc.quadcomp OUT.jsonl
"""
from __future__ import annotations
import json, math, random, statistics, sys, time
from sympy import primerange

CS = (-2, -1, 0, 1, 2)


def sqrt_table(p):
    t = {}
    for x in range(p):
        t.setdefault(x * x % p, []).append(x)
    return t


def preimage(S, c, p, sq):
    out = set()
    for s in S:
        out.update(sq.get((s - c) % p, ()))
    return out


def fibre(seq, p, sq):
    """seq = (c_1..c_k) with F = f_k o ... o f_1; returns |F^{-1}(0)|."""
    T = {0}
    for c in reversed(seq):
        T = preimage(T, c, p, sq)
        if not T:
            return 0
    return len(T)


def all_fibres_depth(k, p, sq):
    """dict seq -> fibre size for all 5^k sequences (DFS over suffixes, shared)."""
    res = {}
    def rec(suffix, T):
        if len(suffix) == k:
            res[tuple(reversed(suffix))] = len(T)
            return
        for c in CS:
            T2 = preimage(T, c, p, sq) if T else set()
            rec(suffix + [c], T2)
    rec([], {0})
    return res


def main(out):
    t0 = time.time()
    rng = random.Random(17)
    primes = sorted(rng.sample(list(primerange(2 ** 12 + 1, 2 ** 13)), 60))
    sqs = {p: sqrt_table(p) for p in primes}
    # (a) exhaustive k = 8
    table = {}
    for p in primes:
        for seq, n in all_fibres_depth(8, p, sqs[p]).items():
            table.setdefault(seq, []).append(n)
    med = {s: statistics.median(v) for s, v in table.items()}
    thr = {p: p / math.log2(p) ** 2 for p in primes}  # prereg threshold p/(log2 p)^2
    def heavy(vals):
        return sum(1 for p, n in zip(primes, vals) if n >= thr[p])
    top_med = sorted(med.items(), key=lambda kv: -kv[1])[:15]
    max_heavy_a = max(heavy(v) for v in table.values())
    # (b) random k = 39
    rng2 = random.Random(39)
    tot = 0; cnt = 0; max_heavy_b = 0; maxfib = 0; best_b = None
    for _ in range(10000):
        seq = tuple(rng2.choice(CS) for _ in range(39))
        vals = [fibre(seq, p, sqs[p]) for p in primes]
        tot += sum(vals); cnt += len(vals)
        h = heavy(vals)
        if h > max_heavy_b:
            max_heavy_b = h
        if max(vals) > maxfib:
            maxfib = max(vals); best_b = (seq, vals)
    res = {"exp": "E4", "primes": primes,
           "a_max_median": top_med[0][1], "a_top_median": [[list(s), m, table[s][:10]] for s, m in top_med],
           "a_max_heavy_count": max_heavy_a,
           "b_mean_fibre": tot / cnt, "b_max_heavy_count": max_heavy_b, "b_max_fibre": maxfib,
           "b_max_fibre_seq": [list(best_b[0]), best_b[1]] if best_b else None,
           "P4a_pass": top_med[0][1] <= 16, "P4b_pass": 0.5 <= tot / cnt <= 2,
           "P4c_pass": max(max_heavy_a, max_heavy_b) < 30, "seconds": round(time.time() - t0, 1)}
    with open(out, "a") as fh:
        fh.write(json.dumps(res) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k not in ("primes",)})[:3000])


if __name__ == "__main__":
    main(sys.argv[1])
