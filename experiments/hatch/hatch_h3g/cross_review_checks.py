#!/usr/bin/env python3
"""
Exact verification and counterexample checks for Round 2 hatch cross-review (G reviewing C-H1, C-H2 + self-audit of G-H1, G-H3).
Read-only against Round 1 outputs; writes nothing outside stdout.
"""
import json
import math
from pathlib import Path

ROOT = Path("/Volumes/git/factoring-connection")


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def primes_upto(B: int):
    return [q for q in range(2, B + 1) if is_prime(q)]


def factor_upto(n: int, primes):
    rem = abs(n)
    fac = {}
    for p in primes:
        if p * p > rem:
            break
        if rem % p == 0:
            cnt = 0
            while rem % p == 0:
                rem //= p
                cnt += 1
            fac[p] = cnt
    if rem > 1 and rem <= primes[-1] and is_prime(rem):
        fac[rem] = fac.get(rem, 0) + 1
        rem = 1
    return fac, rem


def check_ch1_m2_missed_witnesses():
    print("=== CHECK 1: C-H1 M2 missed witnesses (limit = 4*K*K bug) ===")
    runs_path = ROOT / "codex/scratch/hatch_h1c/runs.jsonl"
    K = 64
    limit_codex = 4 * K * K  # 16384
    low_codex = 2 * K * K    # 8192

    total_base = 0
    total_rec = 0
    missed_details = []

    with open(runs_path) as f:
        for line in f:
            row = json.loads(line)
            bits = row["bits"]
            rep = row["rep"]
            N = row["N"]
            a, b, c = row["m2"]["form"]
            B = row["m2"]["B"]
            plist = primes_upto(B)
            # Good primes in C-H1 experiment.py: q > 2, q not dividing a or N, (-N/q) == 1
            good = set()
            for q in plist:
                if q == 2 or a % q == 0 or N % q == 0:
                    continue
                if pow(-N % q, (q - 1) // 2, q) == 1:
                    good.add(q)

            base_hits = []
            rec_hits = []
            for u in range(1, K + 1):
                for v in range(1, K + 1):
                    if math.gcd(u, v) != 1:
                        continue
                    val = a * u * u + b * u * v + c * v * v
                    fac, rem = factor_upto(val, plist)
                    if rem != 1:
                        continue
                    base_hits.append((u, v, val, fac))
                    # Compute good-prime part M_good and check if ANY divisor of M_good lies in (2K^2, 4K^2]
                    good_pows = [(q, e) for q, e in sorted(fac.items()) if q in good]
                    # Note: in C-H1 forced_moduli, if q divides val to power e0, choosing q uses pe = q^1, q^2, ...
                    # Wait: ANY divisor d of M_good is reachable in C-H1's DFS iff we multiply prime powers pe <= q^e_max
                    # without exceeding limit = 4*K*K, stopping as soon as d > 2*K*K!
                    divs = {1}
                    for q, e in good_pows:
                        new_divs = set(divs)
                        for k in range(1, e + 1):
                            qk = q ** k
                            for d0 in list(divs):
                                new_divs.add(d0 * qk)
                        divs = new_divs
                    in_codex_band = sorted(d for d in divs if low_codex < d <= limit_codex)
                    in_lemma_w_band = sorted(d for d in divs if low_codex < d <= low_codex * B)
                    M_good = math.prod(q ** e for q, e in good_pows)
                    if in_codex_band:
                        rec_hits.append((u, v))
                    else:
                        missed_details.append({
                            "bits": bits,
                            "rep": rep,
                            "u": u,
                            "v": v,
                            "Q": val,
                            "B": B,
                            "fac": fac,
                            "M_good": M_good,
                            "has_div_in_2K2_4K2": bool(in_codex_band),
                            "min_div_gt_2K2": min((d for d in divs if d > low_codex), default=None),
                            "max_div_le_2K2": max((d for d in divs if d <= low_codex), default=1),
                            "has_div_in_2K2_2K2B": bool(in_lemma_w_band),
                        })
            total_base += len(base_hits)
            total_rec += len(rec_hits)
            assert len(base_hits) == row["m2"]["base"]["hits"], (bits, rep, len(base_hits), row["m2"]["base"]["hits"])
            assert len(rec_hits) == row["m2"]["recovery"]["hits"], (bits, rep, len(rec_hits), row["m2"]["recovery"]["hits"])

    print(f"Total base hits = {total_base}, recovery hits = {total_rec}, missed = {len(missed_details)}")
    for item in missed_details:
        if item["bits"] >= 64 or len(missed_details) <= 15:
            print("  MISSED:", item)


def check_ch2_grid_entropy_and_det():
    print("\n=== CHECK 2: C-H2 Mahler / grid-entropy floor & determinant check ===")
    res_path = ROOT / "codex/scratch/hatch_h2c/results.jsonl"
    aud_path = ROOT / "codex/scratch/hatch_h2c/audit.jsonl"
    with open(res_path) as fr, open(aud_path) as fa:
        for lr, la in zip(fr, fa):
            rr = json.loads(lr)
            ra = json.loads(la)
            bits = rr["bits"]
            n = rr["n"]
            N = rr["N"]
            best_F = ra["best_norms"][0] * ra["best_norms"][1]
            ratio = math.log(best_F) / math.log(N)
            # Theoretical base Mahler norm of (1, T, T^3): 18/n * log N
            mahler_ratio = 18.0 / n
            # Bulk sublevel heuristic best of M = 29701 samples in D = 6n embeddings, r = 9 real coords:
            # Wait: one scalar norm |F| is product of 2 sides; extreme-value of M samples of sum of 2 sides
            print(f"  bits={bits} n={n} seed={rr['seed']}: best log|F|/logN = {ratio:.4f} vs base 18/n = {mahler_ratio:.3f}")


def check_gh3_dag_pool_ramified_prime_powers():
    print("\n=== CHECK 3: G-H3 Theorem DAG-POOL ramified prime-power check ===")
    # Test W(X, Y) = X^3 * (X + Y) + p^3 * Y^4 for prime p = 101 (D = 4)
    # Count primitive-containing index-p^3 sublattices Lambda_t = {(a,b): a == p*t*b (mod p^3)} with p^3 | W(Lambda_t)
    p = 17
    Q = p ** 3
    # W(X, Y) = X^3 (X + Y) + Q * Y^4
    def W(x, y):
        return (x ** 3) * (x + y) + Q * (y ** 4)

    valid_t = []
    for t in range(p * p):
        # Basis of Lambda_t: M_t(u, v) = (Q*u + p*t*v, v), det M_t = Q
        # Check all coefficients of W(Q*u + p*t*v, v) mod Q:
        ok = True
        for u in range(5):
            for v in range(5):
                if W(Q * u + p * t * v, v) % Q != 0:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            valid_t.append(t)
    print(f"  At p={p}, Q=p^3={Q}, W(X,Y)=X^3(X+Y)+p^3 Y^4: {len(valid_t)} == p^2={p*p} index-p^3 sublattices satisfy Q | W(Lambda_t) (vs D^omega(Q) = 4^1 = 4 for squarefree Q)!")
    # Also check: what power of p divides disc(W)?
    # For quartic a X^4 + b X^3 Y + c X^2 Y^2 + d X Y^3 + e Y^4 with a=1, b=1, c=0, d=0, e=Q:
    a, b, c, d, e = 1, 1, 0, 0, Q
    Delta = (
        256 * (a**3) * (e**3)
        - 192 * (a**2) * b * d * (e**2)
        - 128 * (a**2) * (c**2) * (e**2)
        + 144 * (a**2) * c * (d**2) * e
        - 27 * (a**2) * (d**4)
        + 144 * a * (b**2) * c * (e**2)
        - 6 * a * (b**2) * (d**2) * e
        - 80 * a * b * (c**2) * d * e
        + 18 * a * b * c * (d**3)
        + 16 * a * (c**4) * e
        - 4 * a * (c**3) * (d**2)
        - 27 * (b**4) * (e**2)
        + 18 * (b**3) * c * d * e
        - 4 * (b**3) * (d**3)
        - 4 * (b**2) * (c**3) * e
        + (b**2) * (c**2) * (d**2)
    )
    vp = 0
    tmp = abs(Delta)
    while tmp % p == 0:
        tmp //= p
        vp += 1
    print(f"  disc(W) = {Delta}, v_p(disc(W)) = {vp} (i.e. p^6 | disc(W), so |disc(W)| >= p^6 automatically unless p | N!)")


if __name__ == "__main__":
    check_ch1_m2_missed_witnesses()
    check_ch2_grid_entropy_and_det()
    check_gh3_dag_pool_ramified_prime_powers()
