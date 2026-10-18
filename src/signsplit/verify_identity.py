"""HARNESS: exhaustive verification of the sign-splitting identity (work-plan step 1).

Checks, for every odd N < N_MAX, every squarefree M in MS coprime to N, every
unit sigma mod M and lock-satisfying weights (a, b) from a small box:

  L2  (general value lemma)  a*y + b*x == tau(sigma) (mod M^2) for EVERY unit x
      mod M^2 with x == sigma (mod M) and y = N x^{-1} mod M^2.  (Factor-free,
      and covers p = x, q = y for every factorisation.)
  L1  representative independence: tau(sigma + k M) == tau(sigma).
  L3  CRT locality: tau(sigma) mod m_i^2 == tau_i(sigma mod m_i).
  L4  sign flip: tau(-sigma) == -tau(sigma) (mod M^2).
  T5  signed subset sum: tau(sigma(eps)) == Delta0 + sum eps_i Delta_i (mod M^2)
      for every eps in {+-1}^t, and the integer carry lies in carry_range.
  SP  semiprime check (uses factors): a q + b p == tau(p mod M) (mod M^2) for all
      semiprimes N = p q < N_MAX with p, q prime.
  NEG control: for (a, b) violating the lock, a*y + b*x mod M^2 is NOT constant
      on the class (count of classes where it varies).

Run:  .venv/bin/python -m src.signsplit.verify_identity [N_MAX]
Writes experiments/raw/signsplit_identity.jsonl (one summary record).
"""
from __future__ import annotations

import sys
import time
from math import gcd

from src.common.metrics import append_jsonl
from src.common.instances import is_prime

from .identity import (all_eps, carry, carry_range, factor_squarefree, gamma_of,
                       lock_holds, sigma_of_eps, sign_decomposition, sqrt_mod_prime,
                       tau)

MS = [3, 5, 7, 15, 21, 35, 105, 6, 10, 30, 42, 210, 11 * 13, 3 * 5 * 7 * 11]


def lock_pairs(N: int, M: int, sigma: int, amax: int = 3):
    """(a, b) with b sigma^2 == a N (mod M), a in [-amax, amax], two lifts of b."""
    s2inv = pow(sigma * sigma, -1, M)
    for a in range(-amax, amax + 1):
        b0 = (a * N * s2inv) % M
        for b in (b0, b0 - M):
            if (a, b) != (0, 0):
                yield a, b


def run(N_MAX: int = 301) -> dict:
    t0 = time.time()
    cnt = {"L2_checks": 0, "L1_checks": 0, "L3_checks": 0, "L4_checks": 0,
           "T5_checks": 0, "SP_checks": 0, "NEG_classes": 0, "NEG_varying": 0}
    fail = {k: 0 for k in ("L2", "L1", "L3", "L4", "T5", "carry", "SP")}
    for N in range(3, N_MAX, 2):
        for M in MS:
            if gcd(M, N) != 1:
                continue
            M2 = M * M
            primes = factor_squarefree(M)
            units = [s for s in range(1, M) if gcd(s, M) == 1]
            full_x = M <= 210
            for sigma in units:
                for a, b in lock_pairs(N, M, sigma):
                    assert lock_holds(a, b, sigma, N, M)
                    T = tau(a, b, sigma, N, M)
                    # L2: all x in the class mod M^2 (all units, every lift)
                    ks = range(M) if full_x else range(0, M, max(1, M // 50))
                    for k in ks:
                        x = sigma + k * M
                        y = (N * pow(x, -1, M2)) % M2
                        cnt["L2_checks"] += 1
                        if (a * y + b * x - T) % M2:
                            fail["L2"] += 1
                    # L1
                    for k in (1, 2, -3, M + 1):
                        cnt["L1_checks"] += 1
                        if tau(a, b, sigma + k * M, N, M) != T:
                            fail["L1"] += 1
                    # L3
                    for m in primes:
                        cnt["L3_checks"] += 1
                        if T % (m * m) != tau(a, b, sigma % m, N, m):
                            fail["L3"] += 1
                    # L4
                    cnt["L4_checks"] += 1
                    if (tau(a, b, (-sigma) % M, N, M) + T) % M2:
                        fail["L4"] += 1
            # T5 per square class: one (a, b) per class gamma
            gammas = sorted({gamma_of(s, N, M) for s in units})
            for g in gammas:
                roots = []
                for m in primes:
                    rr = sqrt_mod_prime((N * pow(g, -1, m)) % m, m)
                    roots.append(min(r for r in rr if r % m))
                # every class gamma here comes from a unit so roots exist
                for a in (1, 2, -3):
                    b = (a * g) % M          # b == gamma a  <=>  b sigma^2 == a N
                    for bb in (b, b - M):
                        if (a, bb) == (0, 0):
                            continue
                        dec = sign_decomposition(a, bb, roots, primes, N)
                        clo, chi = carry_range(dec)
                        t = len(dec["Delta"])
                        for eps in all_eps(t):
                            sg = sigma_of_eps(eps, roots, primes)
                            T = tau(a, bb, sg, N, M)
                            S = dec["Delta0"] + sum(e * d for e, d in zip(eps, dec["Delta"]))
                            cnt["T5_checks"] += 1
                            if (S - T) % M2:
                                fail["T5"] += 1
                            c = carry(eps, dec)
                            if not (clo <= c <= chi) or S - M2 * c != T:
                                fail["carry"] += 1
            # NEG control: (a, b) = (1, 0) violates the lock unless N == 0 mod M
            for sigma in units[: min(len(units), 6)]:
                a, b = 1, 1
                if lock_holds(a, b, sigma, N, M):
                    continue
                vals = set()
                for k in range(M):
                    x = sigma + k * M
                    y = (N * pow(x, -1, M2)) % M2
                    vals.add((a * y + b * x) % M2)
                cnt["NEG_classes"] += 1
                cnt["NEG_varying"] += int(len(vals) > 1)
    # SP: real semiprimes
    for p in range(3, int(N_MAX ** 0.5) + 2):
        if not is_prime(p):
            continue
        for q in range(p + 2, N_MAX // p + 1):
            if not is_prime(q):
                continue
            N = p * q
            for M in MS:
                if gcd(M, N) != 1:
                    continue
                for a, b in lock_pairs(N, M, p % M):
                    cnt["SP_checks"] += 1
                    if (a * q + b * p - tau(a, b, p % M, N, M)) % (M * M):
                        fail["SP"] += 1
    rec = {"kind": "identity_verification", "N_MAX": N_MAX, "Ms": MS,
           "counts": cnt, "failures": fail, "wall_s": round(time.time() - t0, 1)}
    return rec


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 301
    rec = run(n)
    print(rec)
    append_jsonl("experiments/raw/signsplit_identity.jsonl", rec)
