"""HARNESS: numerical checks of the Track A propositions (notes/trackA_analysis.md).

  check_A1_exhaustive  Gap lemma / Theorem A1: over F_2 and F_3, the minimum
                       linear complexity (= minimum rank of an algebra B with
                       x in B, lambda: B -> F realising u_i = lambda(x^i)) of a
                       sequence of length L that vanishes off a mask and is
                       nonzero somewhere after position g, exhaustively.
  check_A2_order       Codex r1 Sec. 4.2 / our A2: Berlekamp-Massey order of
                       t -> [X^r]((1+X)^t mod X^M - c) over F_p equals M.
  check_C_props        Props C1-C3 on balanced semiprimes (explicit Lucas support):
                       C1 all classes hit for M < p; C2 <= 3 terms per section
                       for M > p; C3 first level with a pure class vs
                       min(p-1, q-p).
Usage: .venv/bin/python -m src.trackA.theory_checks --out experiments/raw/trackA_theory.jsonl
"""
from __future__ import annotations

import argparse
import itertools
import random
from collections import Counter
from typing import List

from ..common import instances as I
from ..common.metrics import append_jsonl
from .diagnostic import support_explicit


def berlekamp_massey(s: List[int], p: int) -> int:
    """Linear complexity of s over F_p."""
    C, B = [1], [1]
    L, m, b = 0, 1, 1
    for n in range(len(s)):
        d = s[n] % p
        for i in range(1, L + 1):
            d = (d + C[i] * s[n - i]) % p
        if d == 0:
            m += 1
            continue
        coef = d * pow(b, -1, p) % p
        T = C[:]
        if len(C) < len(B) + m:
            C += [0] * (len(B) + m - len(C))
        for i, bi in enumerate(B):
            C[i + m] = (C[i + m] - coef * bi) % p
        if 2 * L <= n:
            L, B, b, m = n + 1 - L, T, d, 1
        else:
            m += 1
    return L


def min_lc_masked(mask: List[int], p: int, after: int) -> int:
    """min linear complexity over u in F_p^L with supp(u) in mask and u_i != 0
    for some i >= after (exhaustive over values on mask positions)."""
    pos = [i for i, b in enumerate(mask) if b]
    L = len(mask)
    best = L + 1
    for vals in itertools.product(range(p), repeat=len(pos)):
        if not any(v and i >= after for v, i in zip(vals, pos)):
            continue
        u = [0] * L
        for v, i in zip(vals, pos):
            u[i] = v
        best = min(best, berlekamp_massey(u, p))
    return best


def check_A1_exhaustive(out):
    L = 20
    cases = []
    for M in (2, 4, 8):
        for r in range(M):
            cases.append((f"class r={r} mod {M}", [int(i % M == r) for i in range(L)], M, M))
    for j in (1, 2):  # (j=3 needs L >= 24 to have a one after the gap)
        cases.append((f"bit {j} = 1", [(i >> j) & 1 for i in range(L)], 1 << (j + 1), (1 << j) + 1))
    for a, b in ((5, 8), (9, 12), (0, 6)):
        # interval [a,b): forward bound a+1 (nontrivial at index >= a) ; Fitting bound b' over a field
        cases.append((f"interval [{a},{b})", [int(a <= i < b) for i in range(L)], a, a + 1))
    for name, mask, after, predicted in cases:
        for p in (2, 3):
            if p == 3 and sum(mask) > 10:
                continue
            got = min_lc_masked(mask, p, after)
            rec = {"experiment": "trackA_theory", "check": "A1_exhaustive", "mask": name,
                   "p": p, "L": L, "nontrivial_after": after, "predicted_lower": predicted,
                   "min_linear_complexity": got, "ok": got >= predicted}
            append_jsonl(out, rec)
            print(rec, flush=True)


def check_A2_order(out):
    rng = random.Random(1)
    for p in (101, 257, 1009):
        for M in (2, 4, 8, 16, 32):
            for _ in range(3):
                r = rng.randrange(M)
                c = rng.randrange(2, p)
                # s_t = [X^r] ((1+X)^t mod X^M - c), t = 0 .. 3M
                v = [1] + [0] * (M - 1)
                seq = []
                for t in range(3 * M + 4):
                    seq.append(v[r])
                    w = [(v[i] + v[i - 1]) % p for i in range(M)]
                    w[0] = (v[0] + c * v[M - 1]) % p
                    v = w
                order = berlekamp_massey(seq, p)
                expect = M if (c - (-1) ** M) % p else M - 1
                rec = {"experiment": "trackA_theory", "check": "A2_order", "p": p, "M": M,
                       "r": r, "c": c, "bm_order": order, "expected_M": M,
                       "ok": order in (M, expect)}
                append_jsonl(out, rec)
    print("A2 done", flush=True)


def check_C_props(out, bits_list=(16, 20, 24, 28, 32, 36, 40), seeds=3):
    for bits in bits_list:
        for seed in range(seeds):
            inst = I.balanced_semiprime(bits, seed)
            N = inst.N
            p, q = inst.factors
            sup = support_explicit(N, [p, q])
            assert len(sup) == 2 * q - p - 1
            first_pure = None
            rows = []
            k = 0
            while (1 << k) <= 2 * N:
                M = 1 << k
                hitp, hitq = Counter(), Counter()
                terms = Counter()
                for i in sup:
                    terms[i % M] += 1
                    if i % p == 0:
                        hitp[i % M] += 1
                    else:
                        hitq[i % M] += 1
                classes = len(terms)
                pure = sum(1 for c in terms if (c in hitp) != (c in hitq))
                maxt = max(terms.values())
                if pure and first_pure is None:
                    first_pure = M
                rows.append({"M": M, "classes_hit": classes, "all_hit": classes == M,
                             "max_terms": maxt, "pure_classes": pure})
                if M < p:
                    assert classes == M, ("C1", N, M)
                if M > p:
                    assert maxt <= 3, ("C2", N, M, maxt)
                k += 1
            rec = {"experiment": "trackA_theory", "check": "C_props", "N": N, "bits": bits,
                   "seed": seed, "p": p, "q": q, "q_minus_p": q - p,
                   "first_pure_M": first_pure, "min_pm1_qmp": min(p - 1, q - p), "rows": rows}
            append_jsonl(out, rec)
            print(bits, seed, "p", p, "q-p", q - p, "first pure M", first_pure, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="experiments/raw/trackA_theory.jsonl")
    a = ap.parse_args()
    check_A2_order(a.out)
    check_C_props(a.out)
    check_A1_exhaustive(a.out)


if __name__ == "__main__":
    main()
