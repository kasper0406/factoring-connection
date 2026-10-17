"""Reproductions of the r3 audit's Track A counterexamples (agent `repair`).  DIAGNOSTIC.

Each function returns a dict of exact results; `main` writes them to JSONL.
Run:  .venv/bin/python -m src.repair.audit_trackA --out experiments/raw/repair_audit_trackA.jsonl
"""
from __future__ import annotations

import argparse
import itertools
import json
from math import comb, gcd, isqrt
from typing import Dict, List

from src.repair.c2r_search import is_prime


# ---------------------------------------------------------------------------
# small helpers (own implementations)
# ---------------------------------------------------------------------------
def coeffs_FN(N: int) -> List[int]:
    """a_i = [X^i] F_N mod N, i = 0..N (a_0 = a_N = 0)."""
    return [0] + [comb(N, i) % N for i in range(1, N)] + [0]


def section(N: int, M: int, r: int) -> Dict[int, int]:
    """S_{r,M}(Y) as {exponent j: coefficient a_{r+jM}} (nonzero terms only)."""
    a = coeffs_FN(N)
    return {(i - r) // M: a[i] for i in range(r, N + 1, M) if a[i]}


def eval_section(S: Dict[int, int], c: int, N: int) -> int:
    return sum(v * pow(c, e, N) for e, v in S.items()) % N


def polymulmod_trunc(f: List[int], g: List[int], K: int, N: int) -> List[int]:
    out = [0] * K
    for i, x in enumerate(f[:K]):
        if x:
            for j, y in enumerate(g[: K - i]):
                out[i + j] = (out[i + j] + x * y) % N
    return out


def berlekamp_massey(seq: List[int], q: int) -> int:
    """Linear complexity over the prime field F_q."""
    s = [x % q for x in seq]
    C, B, L, m, b = [1], [1], 0, 1, 1
    for n in range(len(s)):
        d = s[n]
        for i in range(1, L + 1):
            if i < len(C):
                d = (d + C[i] * s[n - i]) % q
        if d == 0:
            m += 1
            continue
        coef = d * pow(b, -1, q) % q
        T = C[:]
        C = C + [0] * max(0, len(B) + m - len(C))
        for i, bi in enumerate(B):
            C[i + m] = (C[i + m] - coef * bi) % q
        if 2 * L <= n:
            L, B, b, m = n + 1 - L, T, d, 1
        else:
            m += 1
    return L


def balanced_semiprimes(nmax: int, odd: bool = True):
    ps = [x for x in range(2, isqrt(nmax) + 2) if is_prime(x)]
    out = []
    for p in ps:
        for q in range(p + 1, 2 * p):
            if is_prime(q) and p * q <= nmax and (not odd or p > 2):
                out.append((p * q, p, q))
    return sorted(out)


# ---------------------------------------------------------------------------
# A-sel: child selector 2^{k+1} claim; corrected bound D >= min(M, i+1)
# ---------------------------------------------------------------------------
def selector_counterexample() -> Dict:
    N, M, r, K = 35, 16, 5, 6
    a = coeffs_FN(N)
    # observable: coefficient 5 of (.) in R_35[X]/(X^6); L(X^i) = [i == 5]
    L_mono = [1 if i == 5 else 0 for i in range(N + 1)]
    isolates = all(L_mono[i] == 0 for i in range(N + 1) if i % M != r)
    trunc = [1]
    base = [1, 1] + [0] * (K - 2)
    e = N
    while e:  # (1+X)^35 mod (35, X^6) by repeated squaring
        if e & 1:
            trunc = polymulmod_trunc(trunc + [0] * K, base, K, N)
        base = polymulmod_trunc(base, base, K, N)
        e >>= 1
    value = trunc[5]
    return {"item": "A-selector", "N": N, "M": M, "r": r, "a_5": a[5], "observable_rank": K,
            "isolates_class_r_on_[0,N]": isolates, "value_on_F_N": value,
            "rank_below_M": K < M, "support_in_class": [i for i in range(N + 1) if a[i] and i % M == r]}


def corrected_bound_exhaustive(q: int = 2, n: int = 14, Ms=(2, 3, 4, 5, 8)) -> Dict:
    """Every class-supported sequence over F_q has LC >= min(M, i+1) for each nonzero index i,
    and the bound is attained (u = delta_i for i < M; the class indicator for i >= M)."""
    viol = checked = 0
    attained = {}
    for M in Ms:
        for r in range(M):
            pos = [t for t in range(n) if t % M == r]
            for vals in itertools.product(range(q), repeat=len(pos)):
                if not any(vals):
                    continue
                u = [0] * n
                for t, v in zip(pos, vals):
                    u[t] = v
                lc = berlekamp_massey(u, q)
                imax = max(t for t in range(n) if u[t])
                need = min(M, imax + 1)
                checked += 1
                viol += lc < need
    # tightness: delta_r (i = r < M) has LC r+1; the full class indicator has LC M
    tight = True
    for M in Ms:
        for r in range(M):
            delta = [1 if t == r else 0 for t in range(n)]
            ind = [1 if t % M == r else 0 for t in range(n)]
            last = max(t for t in range(n) if ind[t])
            tight &= berlekamp_massey(delta, q) == r + 1 and berlekamp_massey(ind, q) == min(M, last + 1)
    return {"item": "A-corrected-bound", "q": q, "n": n, "sequences": checked,
            "violations": viol, "bound_attained_(delta_r: r+1; class indicator: M)": tight}


def balanced_scoped_bound(nmax: int = 20000) -> Dict:
    """Balanced N = pq (p < q < 2p): every support index is >= p, so an isolating observable
    that is nonzero on F_N has rank >= min(M, p+1).  Also: classes M = 2^k never have M = p."""
    n = bad = 0
    for N, p, q in balanced_semiprimes(nmax):
        supp = [p * j for j in list(range(1, q - p + 1)) + list(range(p, q))] + [q * j for j in range(1, p)]
        n += 1
        bad += min(supp) != p
        bad += any(1 << k == p for k in range(1, 40))
    return {"item": "A-balanced-scoped", "instances": n, "violations": bad}


# ---------------------------------------------------------------------------
def coinflip_example() -> Dict:
    p = 101
    ev = lambda coeffs: sum(c * pow(2, i, p) for i, c in enumerate(coeffs)) % p
    v1, vX = ev([1]), ev([0, 1])
    return {"item": "A-3", "field": p, "f(2) for f=1": v1, "f(2) for f=X": vX,
            "distinguishes": v1 != vX, "isolates_class_0_mod_2": vX == 0,
            "isolates_class_1_mod_2": v1 == 0,
            "budget_selector_is_coin_by_construction": "src/trackA/selectors.py BudgetSelector "
            "returns rng.getrandbits(1) when 2M > D"}


def c2_scope() -> Dict:
    S = section(35, 7, 0)
    gen_viol = gen_n = cop_bad = 0
    for N, p, q in balanced_semiprimes(3000):
        for M in range(p + 1, 2 * q + 1):
            for r in range(M):
                t = len([i for i in range(r, N + 1, M) if comb(N, i) % N and 0 < i < N])
                if gcd(M, N) == 1:
                    gen_n += 1
                    cop_bad += t > 3
                elif t > 3:
                    gen_viol += 1
    return {"item": "A-C2", "N": 35, "M": 7, "r": 0, "section": S, "terms": len(S),
            "coprime_M>p_cases": gen_n, "coprime_M_violations": cop_bad,
            "noncoprime_M_with_>3_terms": gen_viol}


def c3_gap() -> Dict:
    S = section(35, 8, 1)
    return {"item": "A-C3", "N": 35, "M": 8, "r": 1, "section": S,
            "value_at_0": eval_section(S, 0, 35), "gcd_at_0": gcd(eval_section(S, 0, 35), 35),
            "values_all_c": sorted({gcd(eval_section(S, c, 35), 35) for c in range(35)})}


def c4_repair(nmax_eval: int = 1200, nmax_struct: int = 200000) -> Dict:
    S6 = section(6, 2, 0)
    zero_fn_6 = all(eval_section(S6, c, 6) == 0 for c in range(6))
    printed_ineq_fails = []
    for N, p, q in balanced_semiprimes(5000):
        if (q - 1) / 2 + 1 > p - 1:  # printed bound with M = 2
            printed_ineq_fails.append((p, q))
    # structural check of the repaired argument: exponents distinct mod (l-1) on each side
    struct_bad = struct_n = 0
    for N, p, q in balanced_semiprimes(nmax_struct):
        for k in range(1, N.bit_length()):
            M = 1 << k
            span_p = (q - 2) // M
            span_q = (p - 2) // M
            struct_n += 1
            struct_bad += not (span_p <= p - 2 and span_q <= q - 2)
    # direct exhaustive evaluation (componentwise functions on F_l, combined by CRT)
    eval_bad = eval_n = 0
    for N, p, q in balanced_semiprimes(nmax_eval):
        a = coeffs_FN(N)
        for k in range(1, N.bit_length() + 1):
            M = 1 << k
            for r in range(min(M, N + 1)):
                terms = [((i - r) // M, a[i]) for i in range(r, N + 1, M) if a[i]]
                if not terms:
                    continue
                eval_n += 1
                ok = False
                for l in (p, q):
                    if any(v % l for _, v in terms):
                        if any(sum(v * pow(c, e, l) for e, v in terms) % l for c in range(l)):
                            ok = True
                eval_bad += not ok
    return {"item": "A-C4", "N=6 section": S6, "N=6 zero function": zero_fn_6,
            "printed_inequality_fails_(p,q)_first": printed_ineq_fails[:5],
            "printed_inequality_fail_count_N<=5000": len(printed_ineq_fails),
            "span_bound_checks": struct_n, "span_bound_violations": struct_bad,
            "nonzero_sections_evaluated": eval_n, "scalar_false_zeros": eval_bad}


def s04_F51(nmax: int = 1500) -> Dict:
    S = section(51, 4, 0)
    vanishes = all(eval_section(S, c, 51) == 0 for c in range(51))
    # census of scalar traps M = 2^k >= 2 for odd squarefree N <= nmax
    traps = []
    for N in range(15, nmax + 1, 2):
        fac = [l for l in range(3, N + 1) if N % l == 0 and is_prime(l)]
        if len(fac) < 2 or any(N % (l * l) == 0 for l in fac):
            continue
        a = coeffs_FN(N)
        for k in range(1, N.bit_length()):
            M = 1 << k
            for r in range(M):
                terms = [((i - r) // M, a[i]) for i in range(r, N + 1, M) if a[i]]
                if not terms:
                    continue
                if all(all(sum(v * pow(c, e, l) for e, v in terms) % l == 0 for c in range(l)) for l in fac):
                    traps.append((N, M, r))
    balanced = {N for N, _, _ in balanced_semiprimes(nmax)}
    return {"item": "A-S04-F51", "section": S, "zero_function_mod_51": vanishes,
            "trap_count_odd_squarefree_N<=%d" % nmax: len(traps),
            "trap_moduli": sorted({t[0] for t in traps})[:30],
            "traps_on_balanced_semiprimes": sum(t[0] in balanced for t in traps)}


def mask_trailing_block() -> Dict:
    p = 5
    u = [(1 - pow(2, t, p)) % p for t in range(8)]
    return {"item": "A-mask", "field": p, "u_t=1-2^t, t=0..7": u,
            "linear_complexity": berlekamp_massey(u, p),
            "zero_at_4": u[4] == 0, "nonzero_at_3": u[3] != 0, "u_5": u[5],
            "printed_prefix_bound(last_nontrivial+1)": 4}


def circuit_sparse_coordinates(N: int = 35, M: int = 8, c: int = 3) -> Dict:
    """Gates of the X^N chain (monomials X^(2^j)) have ONE nonzero coordinate mod (X^M - c);
    gates (1+X)^(2^j) can also have vanishing coordinates for specific N, c."""
    mono_coords = []
    for j in range(1, 12):
        e = 1 << j
        mono_coords.append(sum(1 for t in range(M) if t == e % M and pow(c, e // M, N)))
    zero_coords = []
    for j in range(1, 12):
        e = 1 << j
        coords = [0] * M
        for i in range(e + 1):
            coords[i % M] = (coords[i % M] + comb(e, i) * pow(c, i // M, N)) % N
        zero_coords.append(sum(1 for x in coords if x == 0))
    return {"item": "A-circuit", "N": N, "M": M, "c": c,
            "nonzero_coords_of_X^(2^j)_j=1..11": mono_coords,
            "zero_coords_of_(1+X)^(2^j)_j=1..11": zero_coords}


def dimension_vs_cost(N: int = 1009 * 1013, M: int = 1 << 12, r: int = 2213, c: int = 7) -> Dict:
    """Rank-M isolating observable lambda = coefficient r in R_N[X]/(X^M - c) applied to f = X^N:
    value c^(N // M) [N mod M == r], computable with O(log N) multiplications."""
    fast = pow(c, N // M, N) if N % M == r else 0
    # dense check on a smaller instance
    Ns, Ms, rs = 1009, 16, 1009 % 16
    dense = [0] * Ms
    dense[1 % Ms] = 1  # X
    acc = [1] + [0] * (Ms - 1)
    e = Ns
    base = dense
    def mul(f, g):
        out = [0] * Ms
        for i, x in enumerate(f):
            if x:
                for j, y in enumerate(g):
                    if y:
                        k = i + j
                        out[k % Ms] = (out[k % Ms] + x * y * pow(c, k // Ms, Ns)) % Ns
        return out
    while e:
        if e & 1:
            acc = mul(acc, base)
        base = mul(base, base)
        e >>= 1
    return {"item": "A-dim-vs-cost", "N": N, "M": M, "r": r, "N_mod_M": N % M, "fast_value": fast,
            "fast_multiplications": 2 * N.bit_length(),
            "small_check_dense_equals_fast": acc[rs] == pow(c, Ns // Ms, Ns)}


def fitting_gap_exhaustive(q: int = 3, D: int = 3, n: int = 14) -> Dict:
    """All sequences with characteristic polynomial x^e g(x) (deg D, g(0) != 0) over F_q:
    D consecutive zeros anywhere => u_t = 0 for all t >= e, including blocks that start before e."""
    viol = early = checked = 0
    for coeffs in itertools.product(range(q), repeat=D):  # x^D + c_{D-1} x^{D-1} + ... + c_0
        e = 0
        while e < D and coeffs[e] == 0:
            e += 1
        for init in itertools.product(range(q), repeat=D):
            u = list(init)
            while len(u) < n:
                t = len(u)
                u.append((-sum(coeffs[i] * u[t - D + i] for i in range(D))) % q)
            for t0 in range(n - D + 1):
                if all(x == 0 for x in u[t0:t0 + D]):
                    checked += 1
                    early += t0 < e
                    viol += any(u[t] for t in range(e, n))
                    break
    return {"item": "A-LemmaG-Fitting", "q": q, "D": D, "n": n, "sequences_with_zero_block": checked,
            "blocks_starting_before_e": early, "violations": viol}


ALL = [selector_counterexample, corrected_bound_exhaustive, balanced_scoped_bound, coinflip_example,
       c2_scope, c3_gap, c4_repair, s04_F51, mask_trailing_block, circuit_sparse_coordinates,
       dimension_vs_cost, fitting_gap_exhaustive]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="experiments/raw/repair_audit_trackA.jsonl")
    a = ap.parse_args(argv)
    with open(a.out, "w") as fh:
        for f in ALL:
            r = f()
            fh.write(json.dumps(r, default=str) + "\n")
            print(json.dumps(r, default=str))


if __name__ == "__main__":
    main()
