"""C-2R mathematical attack: exact checks of the lemmas in notes/repair_C2R_results.md §3.

DIAGNOSTIC (p known).  Functions return dicts; `main` writes JSONL.

  check_character_mechanism   Prop. C2R-2: E = (p+1)/2, p = 3 mod 4  =>  A = H n T, and the
                              boundary form of exactness.
  check_pair_identity         Lemma C2R-1(b): exact  =>  4 r*(c) = |A| + S_A(c) + [c in A] - 2[c/2 in A].
  check_square_exponents      Thm C2R-3: E = 2 (all p) and E = 2 + (p-1)/2 (p = 1 mod 8) have
                              r*(c) <= tau(c0) resp. 2 tau(|c0|) for EVERY c.
  twisted_scaling             twisted identities (p-1) | k(E-1): max_c r*(c) grows linearly in m
                              (no exactness); exact folds among them.
  character_scaling           best exact d / m for the character exponent at larger p.
  payoff_demo                 E with (p-1) | 2E-2 for p | N  =>  gcd(a^(2E-2) - 1, N) reveals p.

Run:  .venv/bin/python -m src.repair.c2r_theory --out experiments/raw/repair_C2R_theory.jsonl
"""
from __future__ import annotations

import argparse
import json
import random
import time
from concurrent.futures import ProcessPoolExecutor
from math import gcd, isqrt
from typing import Dict, List

import numpy as np

from src.repair.c2r_search import criterion_fold, is_prime, params, squares_mask


def legendre(a: int, p: int) -> int:
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def tau(n: int) -> int:
    n = abs(n)
    t, d = 0, 1
    while d * d <= n:
        if n % d == 0:
            t += 1 if d * d == n else 2
        d += 1
    return t


def pair_counts(A: np.ndarray, p: int):
    iu = np.triu_indices(len(A), k=1)
    s = (A[iu[0]] + A[iu[1]]) % p
    return np.unique(s, return_counts=True)


def signed(x: int, p: int) -> int:
    x %= p
    return x - p if x > p // 2 else x


# ---------------------------------------------------------------------------
def check_character_mechanism(pmax: int = 3000) -> Dict:
    n_p = bad = 0
    for p in range(7, pmax, 4):
        if not is_prime(p):
            continue
        n_p += 1
        m, _ = params(p)
        E = (p + 1) // 2
        assert gcd(E, p - 1) == 2
        H = squares_mask(p)
        A = {pow(j, E, p) for j in range(m + 1)}
        HT = {t % p for t in range(-m, m + 1) if H[t % p]}
        bad += A != HT
        for c in range(-2 * m, 2 * m + 1):
            ex, d = criterion_fold(p, E, c % p, m, H)
            if c >= 0:
                outside = [a for a in range(-m, c - m) if H[a % p]]
            else:
                outside = [a for a in range(c + m + 1, m + 1) if H[a % p]]
            ex2 = all(legendre(c - a, p) == -1 for a in outside)
            lo, hi = max(-m, c - m), min(m, c + m)
            pr = sum(1 for a in range(lo, hi + 1) if 2 * a < c and H[a % p] and H[(c - a) % p])
            bad += (ex != ex2) or (d != m + 1 - pr)
    return {"check": "character_mechanism", "primes": n_p, "violations": bad}


def check_pair_identity(primes=(1031, 1543, 4099, 4129, 6151), trials: int = 400, seed: int = 0) -> Dict:
    rng = random.Random(seed)
    n = bad = n_exact = 0
    for p in primes:
        m, _ = params(p)
        H = squares_mask(p)
        Es = [E for E in range(2, p - 1, 2) if gcd(E, p - 1) == 2]
        for _ in range(trials):
            E = rng.choice(Es)
            A = np.array([pow(j, E, p) for j in range(m + 1)], dtype=np.int64)
            Aset = set(A.tolist())
            cs, cnt = pair_counts(A, p)
            for c in list(cs[:50]) + [rng.randrange(p)]:
                c = int(c)
                ex, d = criterion_fold(p, E, c, m, H)
                if not ex:
                    continue
                n_exact += 1
                r = (m + 1) - d
                S = sum(legendre(c - a, p) for a in Aset)
                half = (c * pow(2, -1, p)) % p
                n += 1
                bad += 4 * r != len(Aset) + S + (c in Aset) - 2 * (half in Aset)
    return {"check": "pair_identity_under_exactness", "exact_cases": n_exact, "violations": bad}


def check_square_exponents(pmax: int = 6000) -> Dict:
    n2 = n2c = bad = 0
    worst = 0
    for p in range(11, pmax):
        if not is_prime(p):
            continue
        m, _ = params(p)
        cases = [(2, 1)]
        if p % 8 == 1:
            E = 2 + (p - 1) // 2
            assert gcd(E, p - 1) == 2
            cases.append((E, 2))
        for E, factor in cases:
            A = np.array([pow(j, E, p) for j in range(m + 1)], dtype=np.int64)
            assert len(set(A.tolist())) == m + 1
            cs, cnt = pair_counts(A, p)
            for c, r in zip(cs.tolist(), cnt.tolist()):
                c0 = signed(c, p) if E != 2 else c % p
                bound = factor * tau(c0) if c0 else 0
                bad += r > bound
                worst = max(worst, r)
            if E == 2:
                n2 += 1
            else:
                n2c += 1
    return {"check": "square_exponents_pair_bound", "primes_E2": n2,
            "primes_E2_plus_half": n2c, "violations": bad, "max_pairs_seen": worst}


def _twisted(p: int, k: int):
    """Smallest E = 1 + s(p-1)/k (s coprime to k, s < k) with gcd(E, p-1) = 2, else None."""
    if (p - 1) % k:
        return None
    for s in range(1, k):
        if gcd(s, k) != 1:
            continue
        E = (1 + s * (p - 1) // k) % (p - 1)
        if gcd(E, p - 1) == 2:
            return E
    return None


def twisted_row(args) -> Dict:
    p, k = args
    t0 = time.process_time()
    E = _twisted(p, k)
    if E is None:
        return {"p": p, "k": k, "E": None}
    m, L = params(p)
    H = squares_mask(p)
    A = np.array([pow(j, E, p) for j in range(m + 1)], dtype=np.int64)
    inA = np.zeros(p, dtype=bool)
    inA[A] = True
    cs, cnt = pair_counts(A, p)
    top = np.argsort(-cnt)[:2000]
    B = (cs[top][:, None] - A[None, :]) % p
    exact = ~(H[B] & ~inA[B]).any(axis=1)
    ex_pairs = int(cnt[top][exact].max()) if exact.any() else 0
    # fixed translations c = zeta * c', zeta in the image of theta = y^(E-1), 1 <= c' <= 10:
    # heuristic mean count (2m - c')/(2 k^2) ~ m/k^2 (same-copy pairs u + w = c').
    cdict = dict(zip(cs.tolist(), cnt.tolist()))
    zetas = sorted({pow(j, E - 1, p) for j in range(1, m + 1)})
    fixed = [cdict.get(z * c1 % p, 0) for z in zetas for c1 in range(1, 11)]
    return {"p": p, "k": k, "E": E, "m": m, "L": L, "n_zeta": len(zetas),
            "mean_pairs_fixed_small_c": sum(fixed) / len(fixed),
            "max_pairs_any_c": int(cnt.max()),
            "ratio_pairs_over_m": float(cnt.max()) / m, "heuristic_m_over_k2": m / k ** 2,
            "max_pairs_exact_c_top2000": ex_pairs, "best_exact_d": (m + 1 - ex_pairs),
            "cpu_seconds": time.process_time() - t0}


def char_row(p: int) -> Dict:
    m, L = params(p)
    E = (p + 1) // 2
    H = squares_mask(p)
    best, bc, nex = m + 1, None, 0
    for c in range(1, 2 * m + 1):
        for cc in (c, -c):
            ex, d = criterion_fold(p, E, cc % p, m, H)
            if ex:
                nex += 1
                if d < best:
                    best, bc = d, cc
    return {"p": p, "k": 2, "E": E, "m": m, "L": L, "best_exact_d": best, "c": bc,
            "d_over_m": best / m, "passes_3L8": 8 * best <= 3 * L, "n_exact_c": nex}


def scaling_primes(lo: int, hi: int, per_decade: int, cond) -> List[int]:
    out = []
    x = lo
    while x < hi:
        q = x
        while not (is_prime(q) and cond(q)):
            q += 1
        out.append(q)
        x = int(x * 10 ** (1 / per_decade)) + 1
    return out


def payoff_demo(seed: int = 1) -> Dict:
    rng = random.Random(seed)
    ok = n = 0
    for _ in range(50):
        while True:
            p = rng.randrange(10 ** 5, 10 ** 6) | 3
            if is_prime(p) and p % 4 == 3:
                break
        while True:
            q = rng.randrange(10 ** 5, 10 ** 6) | 1
            if is_prime(q) and q != p:
                break
        N = p * q
        E = (p + 1) // 2
        g = gcd(pow(2, 2 * E - 2, N) - 1, N)
        n += 1
        ok += g == p or (1 < g < N)
    return {"check": "payoff_character_exponent_factors", "instances": n, "split": ok}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="experiments/raw/repair_C2R_theory.jsonl")
    ap.add_argument("--procs", type=int, default=32)
    a = ap.parse_args(argv)
    rows = [check_character_mechanism(), check_pair_identity(), check_square_exponents(), payoff_demo()]
    tw = []
    for k, cond in [(4, lambda q: q % 4 == 1), (6, lambda q: q % 6 == 1), (8, lambda q: q % 8 == 1)]:
        tw += [(q, k) for q in scaling_primes(2 * 10 ** 4, 4 * 10 ** 6, 4, cond)]
    chp = scaling_primes(2 * 10 ** 4, 4 * 10 ** 6, 4, lambda q: q % 4 == 3)
    with ProcessPoolExecutor(max_workers=a.procs) as pool:
        rows += [dict(r, check="twisted_scaling") for r in pool.map(twisted_row, tw)]
        rows += [dict(r, check="character_scaling") for r in pool.map(char_row, chp)]
    with open(a.out, "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
            print(json.dumps(r))


if __name__ == "__main__":
    main()
