"""HARNESS / DIAGNOSTIC experiments for track C mechanisms (owner factC).

This module knows the factors (labels) to classify outcomes.  It is never
imported by operational code.

Usage: .venv/bin/python -m src.trackC.harness_mechanisms --exp E1 [--out ...]
Experiments
  E1  random-shift transfer (M1): hit rate of random probes conditioned on the
      target interval's predicate value, vs accidental-gcd null L/p + L/q.
  E2  power sums (M2): EXACT count over all shifts A mod p of zeros of
      S_e(A, A+L) mod p, split by whether [A, A+L) contains a multiple of p.
  E3  truncated F_N (M3): cost of Z(N,K) vs BGS U(K), and agreement.
  E4  exhaustive check of the C-finite (LRS) obstruction for small p, R.
  E5  exhaustive check of the one-layer power-map obstruction.
  E6  check of the 'least element outside a proper subgroup < sqrt(p)+1' lemma.
  E7  end-to-end least-factor search via U(K) (BGS) on all families.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import random
import time
from math import gcd, isqrt

import gmpy2

from src.common import instances as I
from src.common.metrics import Metrics, append_jsonl
from src.trackC.baselines import bgs_interval, least_factor_by_search, U_predicate
from src.trackC.mechanisms import random_shift_probe, truncated_FN_test

RAW = "experiments/raw/trackC_mechanisms.jsonl"


def contains_multiple(A, L, p):
    """[A, A+L) contains a multiple of p."""
    return (-A) % p < L


# ------------------------------------------------------------------ E1
def E1(out):
    L, t = 2048, 40
    for seed in range(12):
        inst = I.balanced_semiprime(32, seed)
        N, (p, q) = inst.N, inst.labels["factors"]
        rng = random.Random(f"E1|{seed}")
        for want in (True, False):
            while True:
                A = rng.randrange(1, N - L)
                if (contains_multiple(A, L, p) or contains_multiple(A, L, q)) == want:
                    break
            probes = random_shift_probe(N, L, t, seed * 2 + int(want))
            hits = [pr for pr in probes if pr["gcd"] > 1]
            null = 1 - (1 - L / p) * (1 - L / q)
            rec = {"exp": "E1", "seed": seed, "bits": 32, "L": L, "t": t, "target_I": want,
                   "hits": len(hits), "hit_rate": len(hits) / t, "null_rate": null,
                   "hits_explained_by_own_interval": sum(
                       1 for pr in probes if pr["gcd"] > 1 and (contains_multiple(pr["s"], L, p) or contains_multiple(pr["s"], L, q)))}
            append_jsonl(out, rec)
            print(rec)


# ------------------------------------------------------------------ E2
def E2(out):
    rng = random.Random("E2")
    primes = []
    while len(primes) < 8:
        x = rng.randrange(500, 4000)
        if gmpy2.is_prime(x) and x not in primes:
            primes.append(x)
    for p in primes:
        for L in (5, 50, p // 4):
            for e in (1, 2, 3, 5, 8, 13, p - 2, p - 1, p):
                # sliding window over all A mod p:  S(A) = sum_{A<=j<A+L} j^e mod p
                S = sum(pow(j, e, p) for j in range(0, L)) % p
                zin = zout = 0
                for A in range(p):
                    if S == 0:
                        if contains_multiple(A, L, p):
                            zin += 1
                        else:
                            zout += 1
                    S = (S - pow(A, e, p) + pow(A + L, e, p)) % p
                rec = {"exp": "E2", "p": p, "L": L, "e": e, "n_in": L, "n_out": p - L,
                       "zeros_in": zin, "zeros_out": zout,
                       "rate_in": zin / L, "rate_out": zout / (p - L)}
                append_jsonl(out, rec)
                print(rec)


# ------------------------------------------------------------------ E3
def E3(out):
    inst = I.balanced_semiprime(64, 3)
    N, p = inst.N, inst.labels["factors"][0]
    for e in range(6, 17):
        K = 1 << e
        m1, m2 = Metrics(), Metrics()
        t0 = time.perf_counter()
        z = truncated_FN_test(N, K, m1)
        t1 = time.perf_counter()
        u = U_predicate(N, K, "bgs", m2)
        t2 = time.perf_counter()
        rec = {"exp": "E3", "bits": 64, "K": K, "Z": z["Z"], "U": u["U"], "truth": K <= p,
               "fn_q_mul": m1.get("q_mul") + m1.get("q_sqr"), "fn_dimsum": m1.get("q_mul_dimsum"),
               "fn_model": (m1.get("q_mul") + m1.get("q_sqr")) * 2 * K * max(1, (2 * K - 1).bit_length()),
               "bgs_model": m2.get("pmul_model"), "fn_wall": t1 - t0, "bgs_wall": t2 - t1}
        append_jsonl(out, rec)
        print(rec)
    # correctness on small composites, including Carmichael and prime powers
    bad = 0
    for N in list(range(4, 600)) + [561, 1105, 1729, 41041]:
        if gmpy2.is_prime(N):
            continue
        lp = next(d for d in range(2, N + 1) if N % d == 0)
        for K in range(2, min(N, lp + 4)):
            if truncated_FN_test(N, K)["Z"] != (K <= lp):
                bad += 1
    rec = {"exp": "E3_correct", "mismatches": bad}
    append_jsonl(out, rec)
    print(rec)


# ------------------------------------------------------------------ E4
def lrs_seq(c, init, n, p):
    u = list(init)
    R = len(c)
    while len(u) < n:
        u.append(sum(c[i] * u[-R + i] for i in range(R)) % p)
    return u


def E4(out):
    for p, Rmax in ((5, 3), (7, 3), (11, 2)):
        for R in range(1, Rmax + 1):
            n = 2 * p + R + 2
            found_prefix = found_interval = 0
            checked = 0
            for c in itertools.product(range(p), repeat=R):
                for init in itertools.product(range(p), repeat=R):
                    u = lrs_seq(c, init, n, p)          # u[K] for K = 0..n-1
                    checked += 1
                    # prefix pattern on K=1..p+R: nonzero iff K <= p
                    if all((u[K] != 0) == (K <= p) for K in range(1, p + R + 1)):
                        found_prefix += 1
                    # interval pattern (L = R): u_A = 0 iff [A, A+L) contains multiple of p, A=1..2p-1
                    L = R
                    if all((u[A] == 0) == contains_multiple(A, L, p) for A in range(1, 2 * p)):
                        found_interval += 1
            rec = {"exp": "E4", "p": p, "R": R, "sequences": checked,
                   "prefix_pattern_found": found_prefix, "interval_pattern_L_eq_R_found": found_interval}
            append_jsonl(out, rec)
            print(rec)
    # tightness: order L+1 suffices for the interval pattern (u_A = prod_{i<L}(A+i))
    p, L = 11, 3
    u = [math.prod(A + i for i in range(L)) % p for A in range(0, 2 * p + 1)]
    ok = all((u[A] == 0) == contains_multiple(A, L, p) for A in range(1, 2 * p))
    rec = {"exp": "E4_tight", "p": p, "L": L, "order": L + 1, "pattern_ok": ok}
    append_jsonl(out, rec)
    print(rec)


# ------------------------------------------------------------------ E5
def E5(out):
    worst = []
    for p in [x for x in range(5, 140) if gmpy2.is_prime(x)]:
        divs = [g for g in range(1, p) if (p - 1) % g == 0]
        for L in range(4, p - 1):
            Z = {(-i) % p for i in range(L)}            # A with [A,A+L) containing a multiple
            best = None
            for g in divs:
                e = g  # the fibre structure depends only on gcd(e, p-1) = g
                imgs = {}
                for y in range(p):
                    imgs.setdefault(pow(y, e, p), set()).add(y)
                for u in range(p):
                    T = {(z + u) % p for z in Z}        # y = A + u
                    # T must be a union of fibres
                    vals = {pow(y, e, p) for y in T}
                    if all(imgs[v] <= T for v in vals):
                        D = len(vals)
                        if best is None or D < best[0]:
                            best = (D, g, u)
            ratio = best[0] / ((L - 1) / 2)
            worst.append((ratio, p, L, best))
            if best[0] < (L - 1) / 2:
                print("VIOLATION", p, L, best)
    worst.sort()
    rec = {"exp": "E5", "pairs_checked": len(worst), "min_ratio_D_over_(L-1)/2": worst[0][0],
           "argmin": worst[0][1:3], "best_at_argmin": worst[0][3],
           "violations": sum(1 for w in worst if w[0] < 1)}
    append_jsonl(out, rec)
    print(rec)


# ------------------------------------------------------------------ E6
def E6(out):
    worst = 0.0
    arg = None
    for p in range(5, 20000):
        if not gmpy2.is_prime(p):
            continue
        rs = [r for r in range(2, p) if (p - 1) % r == 0 and gmpy2.is_prime(r)]
        nmax = 0
        for r in rs:
            n = next(k for k in range(2, p) if pow(k, (p - 1) // r, p) != 1)
            nmax = max(nmax, n)
        ratio = (nmax - 1) / (math.sqrt(p) + 1)  # max L with [1,L] in a proper subgroup
        if ratio > worst:
            worst, arg = ratio, (p, nmax)
    rec = {"exp": "E6", "p_max": 20000, "max_ratio_Lmax_over_sqrtp+1": worst, "argmax": arg}
    append_jsonl(out, rec)
    print(rec)


# ------------------------------------------------------------------ E7
def E7(out):
    cases = ([I.balanced_semiprime(b, s) for b in (32, 48, 64) for s in range(2)]
             + [I.unbalanced_semiprime(64, 20, 0), I.prime_power(48, 2, 0), I.prime_power(48, 3, 0),
                I.multi_prime(64, 3, 0), I.multi_prime(64, 4, 0), I.carmichael_small(5),
                I.prime_instance(48, 0)])
    for inst in cases:
        m = Metrics()
        t0 = time.perf_counter()
        r = least_factor_by_search(inst.N, "bgs", m)
        dt = time.perf_counter() - t0
        truth = inst.labels["factors"][0] if len(inst.labels["factors"]) > 1 else None
        rec = {"exp": "E7", "family": inst.family, "params": inst.params, "seed": inst.seed,
               "bits": inst.N.bit_length(), "found": r["p"], "correct": r["p"] == truth,
               "calls": r["calls"], "pmul_model": r["counts"].get("pmul_model"),
               "rn_mul": r["counts"].get("rn_mul"), "gcd": r["counts"].get("gcd"), "wall": dt,
               "p_over_N^(1/4)": (truth / inst.N ** 0.25) if truth else None}
        append_jsonl(out, rec)
        print(rec)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", required=True)
    ap.add_argument("--out", default=RAW)
    a = ap.parse_args()
    globals()[a.exp](a.out)



# ------------------------------------------------------------------ E8
def E8(out):
    """Two-layer power maps: P(x) = Q(phi(x+u)), phi(y) = (y^e1 + w)^e2.

    Minimal number of distinct phi-values on an interval T of length L that is a
    union of phi-fibres (= minimal number of distinct roots of Q, a lower bound
    for deg Q).  Compare with one layer (e2 = 1, w = 0) bound ceil((L-1)/2).
    Exhaustive over e1, e2 in [1, p-1], w in F_p and all p cyclic intervals.
    """
    import os
    variant = os.environ.get("E8_VARIANT", "all")   # all | no_zero_exp | gcd_le2
    primes = (7, 11, 13, 17, 19, 23) if variant == "all" else (11, 13, 17, 19, 23, 29, 31)

    def ok_exp(e, p):
        if variant == "all":
            return True
        if variant == "no_zero_exp":
            return e % (p - 1) != 0
        return gcd(e, p - 1) <= 2

    for p in primes:
        best = {L: (L, None) for L in range(2, p - 1)}
        exps = [e for e in range(1, p) if ok_exp(e, p)]
        for e1 in exps:
            for e2 in exps:
                for w in range(p):
                    phi = [pow((pow(y, e1, p) + w) % p, e2, p) for y in range(p)]
                    fib = {}
                    for y, v in enumerate(phi):
                        fib.setdefault(v, []).append(y)
                    for L in range(2, p - 1):
                        if best[L][0] <= 1:
                            continue
                        for a in range(p):
                            T = {(a + i) % p for i in range(L)}
                            vals = {phi[y] for y in T}
                            if len(vals) >= best[L][0]:
                                continue
                            if all(set(fib[v]) <= T for v in vals):
                                best[L] = (len(vals), (e1, e2, w, a))
        for L, (D, arg) in best.items():
            rec = {"exp": "E8", "variant": variant, "p": p, "L": L, "min_distinct_roots_2layer": D,
                   "one_layer_bound": -(-(L - 1) // 2), "argmin(e1,e2,w,a)": arg}
            append_jsonl(out, rec)
            print(rec, flush=True)


if __name__ == "__main__":
    main()
