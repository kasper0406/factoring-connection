"""C-2R: exhaustive search for exact two-layer interval folds (agent `repair`).

DIAGNOSTIC code: the prime p is known.  Independent implementation (not the audit's code).

Model (notes/repair_C2R_prereg.md).  p an odd prime, m = isqrt(p)//2, L = 2m+1,
Z_L = {0, -1, ..., -(L-1)} (mod p), T = Z_L + m = [-m, m].  A two-layer fold is

    phi(x) = ((x + m)^E + v)^2,   Q(Y) = prod_{w in phi(Z_L)} (Y - w),   d = |phi(Z_L)|.

It is EXACT if {x in F_p : Q(phi(x)) = 0} = Z_L.  Restricted to gcd(E, p-1) = 2 with
A = {j^E : 0 <= j <= m}, H = {0} u squares, c = -2v, exactness is the reflection
closure   a in A, c - a in H  ==>  c - a in A,   and  d = |A| - r*(c),  where
r*(c) = #{unordered {a, b} subset A : a != b, a + b = c}.  `brute_fold` checks both
statements by direct evaluation on all of F_p (used in tests for small p).

Arms (minimum exact d):  1 = public exponents 2 <= E <= 64;  2 = all exponents;
3 = all exponents with (2E - 2) % (p-1) != 0 (excludes the character exponent (p+1)/2).
Success per instance and arm: 8 d <= 3 L.

Run:  .venv/bin/python -m src.repair.c2r_search --out experiments/raw/repair_C2R.jsonl
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from math import gcd, isqrt
from typing import Dict, List, Optional, Tuple

import numpy as np

PUBLIC_EMAX = 64
CPU_CAP_DEFAULT = 600.0  # seconds per instance; timeouts are inconclusive


# ---------------------------------------------------------------------------
# instances
# ---------------------------------------------------------------------------
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


def first_primes_from(a: int, k: int = 4) -> List[int]:
    out = []
    while len(out) < k:
        if is_prime(a):
            out.append(a)
        a += 1
    return out


def instance_primes() -> List[int]:
    """First four primes >= 2^j and >= 3*2^(j-1), j = 8, 10, 12, 14 (audit section 4)."""
    out = []
    for mult in (1, 3):
        for j in (8, 10, 12, 14):
            start = (1 << j) if mult == 1 else 3 * (1 << (j - 1))
            out.extend(first_primes_from(start))
    return out


def params(p: int) -> Tuple[int, int]:
    m = isqrt(p) // 2
    return m, 2 * m + 1


def squares_mask(p: int) -> np.ndarray:
    H = np.zeros(p, dtype=bool)
    j = np.arange(p, dtype=np.int64)
    H[(j * j) % p] = True  # includes 0
    return H


def character_exponent(p: int, E: int) -> bool:
    return (2 * E - 2) % (p - 1) == 0


def twist_order(p: int, E: int, kmax: int = 12) -> int:
    """Least k <= kmax with (p-1) | k(E-1), i.e. y^E = chi(y) * y for a character chi of
    order dividing k ('twisted identity'); 0 if none."""
    for k in range(1, kmax + 1):
        if (k * (E - 1)) % (p - 1) == 0:
            return k
    return 0


# ---------------------------------------------------------------------------
# brute force (small p): direct evaluation of the fold on all of F_p
# ---------------------------------------------------------------------------
def brute_fold(p: int, E: int, v: int, m: Optional[int] = None) -> Tuple[bool, int]:
    """Return (exact, d) for phi(x) = ((x+m)^E + v)^2 by evaluating on every x in F_p."""
    if m is None:
        m, _ = params(p)
    L = 2 * m + 1
    vals = [pow((pow((x + m) % p, E, p) + v) % p, 2, p) for x in range(p)]
    ZL = {(-i) % p for i in range(L)}
    R = {vals[x] for x in ZL}
    zero = {x for x in range(p) if vals[x] in R}
    return zero == ZL, len(R)


def criterion_fold(p: int, E: int, c: int, m: Optional[int] = None,
                   H: Optional[np.ndarray] = None) -> Tuple[bool, int]:
    """(exact, d) from the reflection-closure criterion; requires gcd(E, p-1) == 2."""
    assert gcd(E, p - 1) == 2
    if m is None:
        m, _ = params(p)
    if H is None:
        H = squares_mask(p)
    A = {pow(j, E, p) for j in range(m + 1)}
    exact = all((not H[(c - a) % p]) or ((c - a) % p in A) for a in A)
    pairs = sum(1 for a in A if (c - a) % p in A and a < (c - a) % p)
    return exact, len(A) - pairs


# ---------------------------------------------------------------------------
# the exhaustive search
# ---------------------------------------------------------------------------
def search_prime(p: int, cpu_cap: float = CPU_CAP_DEFAULT) -> Dict:
    t0 = time.process_time()
    m, L = params(p)
    H = squares_mask(p)
    base = m + 1
    best = [base, base, base]
    wit: List[Optional[Tuple[int, int]]] = [None, None, None]
    n_E = n_pairsums = n_exact_c = 0
    # diagnostics on non-character exponents
    nc_max_pairs_any = (0, None)      # max r*(c) over all c (exact or not)
    nc_max_pairs_exact = (0, None)    # max r*(c) over exact c
    char_info = None
    iu = np.triu_indices(m + 1, k=1)
    j = np.arange(m + 1, dtype=np.int64)
    for E in range(2, p - 1, 2):
        if gcd(E, p - 1) != 2:
            continue
        if time.process_time() - t0 > cpu_cap:
            return {"p": p, "m": m, "L": L, "status": "inconclusive: CPU cap",
                    "cpu_seconds": time.process_time() - t0}
        n_E += 1
        A = np.array([pow(int(x), E, p) for x in j], dtype=np.int64)
        assert len(np.unique(A)) == m + 1  # fibres of y^E are {+-y}; [0, m] is a transversal
        inA = np.zeros(p, dtype=bool)
        inA[A] = True
        sums = (A[iu[0]] + A[iu[1]]) % p
        cs, cnt = np.unique(sums, return_counts=True)  # cnt = r*(c): distinct a, b
        n_pairsums += len(cs)
        B = (cs[:, None] - A[None, :]) % p
        exact = ~(H[B] & ~inA[B]).any(axis=1)
        n_exact_c += int(exact.sum())
        is_char = character_exponent(p, E)
        if not is_char:
            k = int(cnt.argmax())
            if cnt[k] > nc_max_pairs_any[0]:
                nc_max_pairs_any = (int(cnt[k]), (E, int(cs[k]), twist_order(p, E), bool(exact[k])))
        if exact.any():
            ce = cnt[exact]
            k = int(ce.argmax())
            dmin = base - int(ce[k])
            c = int(cs[exact][k])
            if not is_char and int(ce[k]) > nc_max_pairs_exact[0]:
                nc_max_pairs_exact = (int(ce[k]), (E, c, twist_order(p, E)))
            if is_char:
                char_info = {"E": E, "d": dmin, "c": c, "n_exact_c": int(exact.sum())}
            enabled = [E <= PUBLIC_EMAX, True, not is_char]
            for arm in range(3):
                if enabled[arm] and dmin < best[arm]:
                    best[arm] = dmin
                    wit[arm] = (E, (-c * pow(2, -1, p)) % p)
        elif is_char:
            char_info = {"E": E, "d": None, "c": None, "n_exact_c": 0}
    return {
        "p": p, "m": m, "L": L, "p_mod_4": p % 4, "status": "ok",
        "best_d": best, "witness_E_v": wit,
        "passes": [8 * d <= 3 * L for d in best],
        "threshold_3L_over_8": 3 * L / 8, "one_layer_d": base,
        "n_exponents": n_E, "n_candidate_c": n_pairsums, "n_exact_c": n_exact_c,
        "character_exponent": char_info,
        "noncharacter_max_pairs_any_c": {"pairs": nc_max_pairs_any[0],
                                         "E_c_twist_exact": nc_max_pairs_any[1]},
        "noncharacter_max_pairs_exact_c": {"pairs": nc_max_pairs_exact[0],
                                           "E_c_twist": nc_max_pairs_exact[1]},
        "cpu_seconds": time.process_time() - t0,
        "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }


def verify_certificate(p: int = 1031, E: int = 516, v: int = -7) -> Dict:
    """Exhaustive check of the audit's certificate: zero set of Q(phi(x)) over all x in F_p."""
    m, L = params(p)
    R_claim = {1, 4, 9, 16, 25, 36, 49, 64, 81, 196, 441}
    phi = [pow((pow((x + m) % p, E, p) + v) % p, 2, p) for x in range(p)]
    ZL = {(-i) % p for i in range(L)}
    R = {phi[x] for x in ZL}
    Qzero = {x for x in range(p) if phi[x] in R_claim}
    return {"p": p, "m": m, "L": L, "E": E, "v": v, "E_is_(p+1)/2": E == (p + 1) // 2,
            "roots_equal_claim": R == R_claim, "deg_Q": len(R_claim),
            "zero_set_is_Z_L": Qzero == ZL, "n_inputs_checked": p}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="experiments/raw/repair_C2R.jsonl")
    ap.add_argument("--procs", type=int, default=32)
    ap.add_argument("--cpu-cap", type=float, default=CPU_CAP_DEFAULT)
    ap.add_argument("--primes", type=int, nargs="*")
    a = ap.parse_args(argv)
    primes = a.primes or instance_primes()
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.procs) as pool, open(a.out, "a") as fh:
        futs = [pool.submit(search_prime, p, a.cpu_cap) for p in primes]
        for f in futs:
            row = f.result()
            row["experiment"] = "C2R_main"
            fh.write(json.dumps(row) + "\n")
            fh.flush()
            print(json.dumps({k: row.get(k) for k in ("p", "L", "best_d", "passes", "cpu_seconds")}),
                  flush=True)
    print(f"wall {time.time() - t0:.1f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
