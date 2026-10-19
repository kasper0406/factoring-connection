"""E1 (notes/polyimpl_prereg.md): residual mass of cyclotomic group-order exposure.

DIAGNOSTIC code. It factors the cyclotomic values Phi_d(p) of sampled primes p, never a
composite N = pq. No operational algorithm imports it.

A prime p is "exposed by d" iff P+(Phi_d(p)) <= B2, the most generous criterion: some
B2-smooth exponent annihilates the whole torus T_d(F_p) of order Phi_d(p).
"""
from __future__ import annotations

import json
import math
import random
import sys
from fractions import Fraction

import flint
import gmpy2

D_SET = (1, 2, 3, 4, 6)
SIZES = (20, 26, 32)
SAMPLES = 20000


def phi_d(d: int, p: int) -> int:
    """Cyclotomic polynomial Phi_d evaluated at p, for d in {1,2,3,4,6}."""
    if d == 1:
        return p - 1
    if d == 2:
        return p + 1
    if d == 3:
        return p * p + p + 1
    if d == 4:
        return p * p + 1
    if d == 6:
        return p * p - p + 1
    raise ValueError(d)


def largest_prime_factor(m: int) -> int:
    if m <= 1:
        return 1
    return max(int(f) for f, _ in flint.fmpz(m).factor())


def b2_bound(k: int) -> int:
    """B2 = X^2 with X = 2^{k/3} (N ~ p^2 ~ 2^{2k}, X = N^{1/6}); floor of 2^{2k/3}."""
    # exact floor of 2^{2k/3}: largest integer b with b^3 <= 2^{2k}
    target = 1 << (2 * k)
    b = int(round(2 ** (2 * k / 3)))
    while b ** 3 > target:
        b -= 1
    while (b + 1) ** 3 <= target:
        b += 1
    return b


def exposure_profile(p: int, b2: int) -> dict:
    return {d: largest_prime_factor(phi_d(d, p)) <= b2 for d in D_SET}


def sample_primes(k: int, count: int) -> list[int]:
    rng = random.Random(1019 * k)
    lo, hi = 1 << (k - 1), (1 << k) - 1
    seen: set[int] = set()
    out: list[int] = []
    while len(out) < count:
        x = rng.randint(lo, hi)
        if x in seen:
            continue
        if gmpy2.is_prime(x):
            seen.add(x)
            out.append(x)
    return out


def run_size(k: int, count: int = SAMPLES) -> dict:
    b2 = b2_bound(k)
    primes = sample_primes(k, count)
    recs = []
    for p in primes:
        prof = exposure_profile(p, b2)
        recs.append({"k": k, "p": p, **{f"exp{d}": prof[d] for d in D_SET}})
    resid = [not any(r[f"exp{d}"] for d in D_SET) for r in recs]
    resid1 = [not r["exp1"] for r in recs]
    resid12 = [not (r["exp1"] or r["exp2"]) for r in recs]
    rng = random.Random(7 * k)
    idx = list(range(len(primes)))
    rng.shuffle(idx)
    pairs = []
    for i in range(0, len(idx) - 1, 2):
        a, b = idx[i], idx[i + 1]
        p, q = sorted((primes[a], primes[b]))
        if q <= 2 * p:
            pairs.append((resid[a], resid[b]))  # "both residual" is symmetric
    r2 = sum(1 for x, y in pairs if x and y) / len(pairs)
    n = len(recs)

    def se(f):
        return math.sqrt(f * (1 - f) / n)

    r1 = sum(resid) / n
    summ = {
        "k": k, "B2": b2, "n_primes": n, "n_pairs": len(pairs),
        "r1_D1": sum(resid1) / n, "r1_D12": sum(resid12) / n, "r1": r1,
        "se_r1": se(r1), "r2": r2, "se_r2": math.sqrt(r2 * (1 - r2) / len(pairs)),
        "per_d_unexposed": {d: sum(1 for r in recs if not r[f"exp{d}"]) / n for d in D_SET},
    }
    return {"summary": summ, "records": recs}


def model_values() -> dict:
    """Independent Dickman-model predictions (prereg P1-P4)."""
    rho3 = 0.04860838  # Dickman rho(3)
    u1 = math.log(1.5)  # Pr[P+(n) > n^{2/3}]
    return {"r1_D1": u1, "r1_D12": u1 ** 2, "r1": u1 ** 2 * (1 - rho3) ** 3,
            "r2": (u1 ** 2 * (1 - rho3) ** 3) ** 2}


def run(out_raw="experiments/raw/polyimpl_e1.jsonl",
        out_sum="experiments/analysis/polyimpl_e1_summary.json"):
    summaries = []
    with open(out_raw, "w") as fh:
        for k in SIZES:
            res = run_size(k)
            for r in res["records"]:
                fh.write(json.dumps(r) + "\n")
            summaries.append(res["summary"])
            print(json.dumps(res["summary"]), flush=True)
    s = {s["k"]: s for s in summaries}
    checks = {
        "P1": all(0.33 <= x["r1_D1"] <= 0.48 for x in summaries),
        "P2": all(0.12 <= x["r1_D12"] <= 0.22 for x in summaries),
        "P3": all(0.09 <= x["r1"] <= 0.20 for x in summaries),
        "P4": all(0.008 <= x["r2"] <= 0.045 for x in summaries),
        "P5": 0.7 <= s[32]["r1"] / s[20]["r1"] <= 1.3,
    }
    verdict = ("CB+ experimentally supported; Rank-4 refuted for cyclotomic portfolios"
               if checks["P3"] and checks["P4"] and checks["P5"]
               else ("residual shrinks; no CB+ support"
                     if s[32]["r1"] < 0.05 or s[32]["r1"] / s[20]["r1"] < 0.7
                     else "mixed; see numbers"))
    with open(out_sum, "w") as fh:
        json.dump({"summaries": summaries, "model": model_values(), "checks": checks,
                   "verdict": verdict}, fh, indent=2)
    print(json.dumps({"checks": checks, "verdict": verdict}))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        run()
    else:
        print(__doc__)
