#!/usr/bin/env python3
"""
Pre-registered toy experiments for Hatch 1:
  E1: Exact index & LLL decoding comparison: Z[Y] (Boneh/HG) vs Int(Z) (Polya/Mahler, containing all Buium p-derivations) vs Linear Buium delta_2 jet lattice.
  E2: Local q-derivation redundancy at w=0 and CRT-gluing height blowup across coprime primes.
  E3: Galois-folded O_K CRT decoding beyond Johnson over K = Q(zeta_5) vs smooth-norm split-prime failure (A_orbit = 0).
  E4: Toy M2 quadratic-form smooth-search wall (Z[k], Int(Z), and linear delta_2).
"""

from __future__ import annotations
import json
import math
import random
import statistics
import time
from fractions import Fraction
from pathlib import Path

OUT_DIR = Path(__file__).resolve().parent


# ============================================================================
# 1. Exact Arithmetic Helpers: Primes, Binomials, Buium p-Derivations, Modular HNF, Exact Rational LLL
# ============================================================================

def primes_up_to(n: int) -> list[int]:
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            sieve[p * p : n + 1 : p] = b"\x00" * (((n - p * p) // p) + 1)
    return [i for i, v in enumerate(sieve) if v]


def v_p(n: int, p: int) -> int:
    if n == 0:
        return 10**9
    c = 0
    n = abs(n)
    while n % p == 0:
        n //= p
        c += 1
    return c


def delta_p(x: int, p: int) -> int:
    """Buium's p-derivation (Fermat quotient): delta_p(x) = (x - x^p) // p."""
    num = x - pow(x, p)
    assert num % p == 0
    return num // p


def delta_p_iter(x: int, p: int, order: int) -> int:
    val = x
    for _ in range(order):
        val = delta_p(val, p)
    return val


def delta_p_mod(x: int, p: int, order: int, mod: int) -> int:
    """Evaluate delta_p^{(order)}(x) mod `mod` when gcd(p, mod) == 1."""
    inv_p = pow(p, -1, mod)
    val = x % mod
    for _ in range(order):
        val = ((val - pow(val, p, mod)) * inv_p) % mod
    return val


def binom_int(n: int, k: int) -> int:
    """Exact integer-valued polynomial binom(n, k) = n(n-1)...(n-k+1)/k! for any n in Z, k >= 0."""
    if k < 0:
        return 0
    if k == 0:
        return 1
    num = 1
    for i in range(k):
        num *= (n - i)
    return num // math.factorial(k)


def ext_gcd(a: int, b: int) -> tuple[int, int, int]:
    """Return (g, x, y) with a*x + b*y = g >= 0."""
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r < 0:
        return -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def add_congruence_to_basis(basis: list[list[int]], coeffs: list[int], mod: int) -> tuple[list[list[int]], int]:
    """
    Given a Z-basis `basis` (d rows in Z^d) of a lattice L, intersect L with
    { x in Z^d : sum_j coeffs[j] * x[j] == 0 (mod `mod`) }.
    Returns (new_basis, index_factor).
    """
    d = len(basis)
    vals = [sum(row[j] * coeffs[j] for j in range(d)) % mod for row in basis]
    # Reduce vals using row operations over Z so only at most one nonzero residue remains
    B = [row[:] for row in basis]
    piv = -1
    for i in range(d):
        if vals[i] == 0:
            continue
        if piv == -1:
            piv = i
        else:
            a, b = vals[piv], vals[i]
            g, s, t = ext_gcd(a, b)
            u, w = -b // g, a // g
            # New piv row = s*B[piv] + t*B[i] (residue g)
            # New i row   = u*B[piv] + w*B[i] (residue 0)
            row_piv = [s * B[piv][k] + t * B[i][k] for k in range(d)]
            row_i = [u * B[piv][k] + w * B[i][k] for k in range(d)]
            B[piv], B[i] = row_piv, row_i
            vals[piv], vals[i] = g % mod, 0
    if piv == -1 or vals[piv] == 0:
        return B, 1
    g_final = math.gcd(vals[piv], mod)
    idx = mod // g_final
    B[piv] = [x * idx for x in B[piv]]
    return B, idx


def hnf_reduce_rows(basis: list[list[int]]) -> list[list[int]]:
    """Row-Hermite Normal Form of full-rank d x d integer basis to keep entries compact before LLL."""
    d = len(basis)
    B = [row[:] for row in basis]
    for col in range(d):
        # Find pivot in rows >= col
        piv = -1
        for r in range(col, d):
            if B[r][col] != 0:
                if piv == -1:
                    piv = r
                else:
                    a, b = B[piv][col], B[r][col]
                    g, s, t = ext_gcd(a, b)
                    u, w = -b // g, a // g
                    rp = [s * B[piv][k] + t * B[r][k] for k in range(d)]
                    rr = [u * B[piv][k] + w * B[r][k] for k in range(d)]
                    B[piv], B[r] = rp, rr
        assert piv != -1, "Basis not full rank"
        if piv != col:
            B[col], B[piv] = B[piv], B[col]
        if B[col][col] < 0:
            B[col] = [-x for x in B[col]]
        diag = B[col][col]
        for r in range(col):
            q = B[r][col] // diag
            if q != 0:
                B[r] = [B[r][k] - q * B[col][k] for k in range(d)]
    return B


def exact_lll(basis_int: list[list[int]]) -> list[list[int]]:
    """Exact integral LLL (Cohen Algorithm 2.6.7, delta = 3/4) — zero Fractions, exact integer arithmetic."""
    n = len(basis_int)
    m = len(basis_int[0])
    B = [row[:] for row in basis_int]

    def dot(u: list[int], v: list[int]) -> int:
        return sum(a * b for a, b in zip(u, v))

    d = [0] * (n + 1)
    d[0] = 1
    lam = [[0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1):
            u = dot(B[i], B[j])
            for l in range(j):
                u = (d[l + 1] * u - lam[i][l] * lam[j][l]) // d[l]
            if j < i:
                lam[i][j] = u
            else:
                d[i + 1] = u

    def redi(k_idx: int, l_idx: int) -> None:
        if 2 * abs(lam[k_idx][l_idx]) > d[l_idx + 1]:
            q = (2 * lam[k_idx][l_idx] + d[l_idx + 1]) // (2 * d[l_idx + 1])
            B[k_idx] = [B[k_idx][t] - q * B[l_idx][t] for t in range(m)]
            lam[k_idx][l_idx] -= q * d[l_idx + 1]
            for i_idx in range(l_idx):
                lam[k_idx][i_idx] -= q * lam[l_idx][i_idx]

    def swapi(k_idx: int) -> None:
        B[k_idx], B[k_idx - 1] = B[k_idx - 1], B[k_idx]
        for j_idx in range(k_idx - 1):
            lam[k_idx][j_idx], lam[k_idx - 1][j_idx] = lam[k_idx - 1][j_idx], lam[k_idx][j_idx]
        L_val = lam[k_idx][k_idx - 1]
        B_new = (d[k_idx - 1] * d[k_idx + 1] + L_val * L_val) // d[k_idx]
        for i_idx in range(k_idx + 1, n):
            t_val = lam[i_idx][k_idx]
            lam[i_idx][k_idx] = (d[k_idx + 1] * lam[i_idx][k_idx - 1] - L_val * t_val) // d[k_idx]
            lam[i_idx][k_idx - 1] = (B_new * t_val + L_val * lam[i_idx][k_idx]) // d[k_idx + 1]
        d[k_idx] = B_new

    k = 1
    while k < n:
        redi(k, k - 1)
        if 4 * (d[k + 1] * d[k - 1] + lam[k][k - 1] * lam[k][k - 1]) >= 3 * d[k] * d[k]:
            for l in range(k - 2, -1, -1):
                redi(k, l)
            k += 1
        else:
            swapi(k)
            k = max(k - 1, 1)
    return B


def crt_pair(a1: int, m1: int, a2: int, m2: int) -> tuple[int, int]:
    g, s, _ = ext_gcd(m1, m2)
    assert (a2 - a1) % g == 0
    lcm = (m1 // g) * m2
    x = (a1 + m1 * (((a2 - a1) // g) * s % (m2 // g))) % lcm
    return x, lcm


def crt_list(residues: list[int], moduli: list[int]) -> tuple[int, int]:
    r, m = 0, 1
    for a, mod in zip(residues, moduli):
        r, m = crt_pair(r, m, a % mod, mod)
    return r, m


# ============================================================================
# 2. Experiment E1: Int(Z) vs Z[Y] vs Linear Buium delta_2
# ============================================================================

def build_intZ_J_lattice(P: int, w: int, r: int, L: int) -> tuple[list[list[int]], int]:
    """J_{P, r, L}: F(Y) = sum_{k=0}^L b_k binom(Y, k) with P^r | F(w + P*u) for u=0..L."""
    d = L + 1
    B = [[1 if i == j else 0 for j in range(d)] for i in range(d)]
    total_idx = 1
    mod = P**r
    for u in range(L + 1):
        pt = w + P * u
        coeffs = [binom_int(pt, k) % mod for k in range(d)]
        B, idx = add_congruence_to_basis(B, coeffs, mod)
        B = hnf_reduce_rows(B)
        total_idx *= idx
    return B, total_idx


def build_intZ_Istar_lattice(prime_powers: list[tuple[int, int]], w: int, r: int, L: int) -> tuple[list[list[int]], int]:
    """I^*_{P, r, L}: for each q^e || P and 1 <= f <= e, q^(f*r) | F((w mod q^f) + q^f * u) for u=0..L."""
    d = L + 1
    B = [[1 if i == j else 0 for j in range(d)] for i in range(d)]
    total_idx = 1
    for q, e in prime_powers:
        for f in range(1, e + 1):
            qf = q**f
            qfr = q ** (f * r)
            wf = w % qf
            for u in range(L + 1):
                pt = wf + qf * u
                coeffs = [binom_int(pt, k) % qfr for k in range(d)]
                B, idx = add_congruence_to_basis(B, coeffs, qfr)
                total_idx *= idx
            B = hnf_reduce_rows(B)
    return B, total_idx


def build_boneh_ZY_lattice(P: int, w: int, r: int, L: int) -> list[list[int]]:
    """Standard Boneh/Howgrave-Graham basis in Z[Y]_{<= L}: P^{max(0, r-k)} (Y - w)^k."""
    d = L + 1
    rows = []
    for k in range(d):
        mult = P ** max(0, r - k)
        # Expand mult * (Y - w)^k = sum_{j=0}^k mult * C(k, j) * (-w)^(k-j) * Y^j
        row = [0] * d
        for j in range(k + 1):
            row[j] = mult * math.comb(k, j) * ((-w) ** (k - j))
        rows.append(row)
    return hnf_reduce_rows(rows)


def eval_intZ_poly(b: list[int], y: int) -> int:
    return sum(b[k] * binom_int(y, k) for k in range(len(b)))


def eval_ZY_poly(c: list[int], y: int) -> int:
    return sum(c[k] * (y**k) for k in range(len(c)))


def run_E1() -> dict:
    rng = random.Random(20261120)
    ps_10 = primes_up_to(29)  # 10 primes
    P_10 = math.prod(ps_10)
    pp_10 = [(q, 1) for q in ps_10]

    # E1a: Exact index checks
    index_rows = []
    all_p1a = True
    for r, L in [(1, 2), (1, 4), (2, 4), (2, 6), (3, 6)]:
        w = rng.randrange(P_10)
        _, idx_J = build_intZ_J_lattice(P_10, w, r, L)
        _, idx_Istar = build_intZ_Istar_lattice(pp_10, w, r, L)
        copp_idx = P_10 ** (r * (r + 1) // 2)
        holds_chain = (idx_Istar >= idx_J >= copp_idx)
        r1_exact = True
        if r == 1:
            lucas_pred = math.prod(q ** ((L // q) + 1) for q in ps_10)
            r1_exact = (idx_Istar == lucas_pred)
        all_p1a = all_p1a and holds_chain and r1_exact
        index_rows.append({
            "r": r,
            "L": L,
            "ln_idx_Coppersmith": round(math.log(copp_idx), 4),
            "ln_idx_J": round(math.log(idx_J), 4),
            "ln_idx_Istar": round(math.log(idx_Istar), 4),
            "holds_chain": holds_chain,
            "r1_lucas_exact": r1_exact,
        })

    # E1b: Empirical LLL decoding radius across 30 planted trials (P = prod_{q<=23} q, H = 2^12, L=4, r=2)
    ps_9 = primes_up_to(23)  # 9 primes: 2,3,5,7,11,13,17,19,23; P = 223,092,870
    P_9 = math.prod(ps_9)
    pp_9 = [(q, 1) for q in ps_9]
    H_1b = 2**12
    L_1b, r_1b = 4, 2
    W_9 = math.log(P_9)
    b_1b = math.log(H_1b)
    johnson_9 = math.sqrt(W_9 * b_1b)

    # Generate all 512 divisor subsets sorted by ln(d)
    divisors = [1]
    for q in ps_9:
        divisors += [d * q for d in divisors]
    divisors = sorted(set(divisors))
    # Pick candidate divisors spanning [0.5 * johnson, 1.35 * johnson]
    cand_divs = [d for d in divisors if 0.55 * johnson_9 <= math.log(d) <= 1.30 * johnson_9]
    # Subsample ~20 evenly spaced in log(d)
    step = max(1, len(cand_divs) // 20)
    cand_divs = cand_divs[::step]

    min_ZY, min_Istar, min_J = [], [], []
    scales_ZY = [H_1b**k for k in range(L_1b + 1)]
    scales_intZ = [max(1, binom_int(H_1b, k)) for k in range(L_1b + 1)]

    for trial in range(30):
        j0 = rng.randint(1, H_1b - 1)
        got_zy = got_is = got_j = float("inf")
        for d in cand_divs:
            # Plant X such that gcd(X + j0, P_9) == d
            # Choose cofactor u coprime to P_9
            while True:
                u = rng.randrange(1, P_9)
                if math.gcd(u, P_9 // d) == 1:
                    break
            X_mod = (d * u - j0) % P_9
            w = (-X_mod) % P_9
            lnd = math.log(d)

            if got_zy == float("inf"):
                B_zy = build_boneh_ZY_lattice(P_9, w, r_1b, L_1b)
                B_zy_sc = [[row[k] * scales_ZY[k] for k in range(L_1b + 1)] for row in B_zy]
                red = exact_lll(B_zy_sc)
                c_vec = [red[0][k] // scales_ZY[k] for k in range(L_1b + 1)]
                if any(c_vec) and eval_ZY_poly(c_vec, j0) == 0:
                    got_zy = lnd

            if got_is == float("inf"):
                B_is, _ = build_intZ_Istar_lattice(pp_9, w, r_1b, L_1b)
                B_is_sc = [[row[k] * scales_intZ[k] for k in range(L_1b + 1)] for row in B_is]
                red = exact_lll(B_is_sc)
                b_vec = [red[0][k] // scales_intZ[k] for k in range(L_1b + 1)]
                if any(b_vec) and eval_intZ_poly(b_vec, j0) == 0:
                    got_is = lnd

            if got_j == float("inf"):
                # Note: J_{P, r, L} only enforces P^r | F(w + P*u), NOT d^r | F(w + d*u)!
                # Still test if its shortest vector happens to vanish at j0!
                B_j, _ = build_intZ_J_lattice(P_9, w, r_1b, L_1b)
                B_j_sc = [[row[k] * scales_intZ[k] for k in range(L_1b + 1)] for row in B_j]
                red = exact_lll(B_j_sc)
                b_vec = [red[0][k] // scales_intZ[k] for k in range(L_1b + 1)]
                if any(b_vec) and eval_intZ_poly(b_vec, j0) == 0:
                    got_j = lnd

        min_ZY.append(got_zy)
        min_Istar.append(got_is)
        min_J.append(got_j)

    mean_ZY = statistics.mean(x for x in min_ZY if x < float("inf"))
    mean_Istar = statistics.mean(x for x in min_Istar if x < float("inf"))
    finite_J = [x for x in min_J if x < float("inf")]
    mean_J = statistics.mean(finite_J) if finite_J else W_9
    p1b_holds = (
        min(min_Istar) >= 0.85 * min(min_ZY)
        and min(min_J) >= 0.85 * min(min_ZY)
        and min(min_Istar) > 0.70 * johnson_9
        and min(min_J) > 0.70 * johnson_9
    )

    # E1c: Linear Buium delta_2 vs r=1 Coppersmith (s=3, 4; H=2^10; P_odd = prod_{3<=q<=29} q)
    ps_odd = [q for q in primes_up_to(29) if q >= 3]
    P_odd = math.prod(ps_odd)
    H_1c = 2**10
    divs_odd = [1]
    for q in ps_odd:
        divs_odd += [d * q for d in divs_odd]
    divs_odd = sorted(set(divs_odd))
    step_odd = max(1, len(divs_odd) // 32)
    cand_odd = divs_odd[::step_odd]

    lin_comp = []
    p1c_holds = True
    for s in (2, 3, 4):
        buium_scales = [1] + [max(1, abs(delta_p_iter(H_1c, 2, u))) for u in range(s)]
        copp_scales = [H_1c**u for u in range(s + 1)]
        b_mins, c_mins = [], []
        for trial in range(15):
            j0 = rng.randint(3, H_1c - 1)
            gb = gc = math.log(P_odd)
            for d in cand_odd:
                while True:
                    u = rng.randrange(1, P_odd)
                    if math.gcd(u, P_odd // d) == 1:
                        break
                w = (j0 - d * u) % P_odd  # w == j0 mod d
                lnd = math.log(d)
                if gb == math.log(P_odd):
                    d_dim = s + 1
                    B0 = [[1 if i == j else 0 for j in range(d_dim)] for i in range(d_dim)]
                    coeffs = [1] + [delta_p_mod(w, 2, u_idx, P_odd) for u_idx in range(s)]
                    B_lat, _ = add_congruence_to_basis(B0, coeffs, P_odd)
                    B_sc = [[row[k] * buium_scales[k] for k in range(d_dim)] for row in B_lat]
                    red = exact_lll(B_sc)
                    c_vec = [red[0][k] // buium_scales[k] for k in range(d_dim)]
                    val = c_vec[0] + sum(c_vec[u_idx + 1] * delta_p_iter(j0, 2, u_idx) for u_idx in range(s))
                    if any(c_vec) and val == 0:
                        gb = lnd
                if gc == math.log(P_odd):
                    d_dim = s + 1
                    B0 = [[1 if i == j else 0 for j in range(d_dim)] for i in range(d_dim)]
                    coeffs = [pow(w, u_idx, P_odd) for u_idx in range(s + 1)]
                    B_lat, _ = add_congruence_to_basis(B0, coeffs, P_odd)
                    B_sc = [[row[k] * copp_scales[k] for k in range(d_dim)] for row in B_lat]
                    red = exact_lll(B_sc)
                    c_vec = [red[0][k] // copp_scales[k] for k in range(d_dim)]
                    val = sum(c_vec[u_idx] * (j0**u_idx) for u_idx in range(s + 1))
                    if any(c_vec) and val == 0:
                        gc = lnd
            b_mins.append(gb)
            c_mins.append(gc)
        mean_b = statistics.mean(b_mins)
        mean_c = statistics.mean(c_mins)
        if s >= 3 and not (mean_b > mean_c):
            p1c_holds = False
        lin_comp.append({
            "s": s,
            "theory_b_coeff_Buium": round((2**s - 1) / (s + 1), 4),
            "theory_b_coeff_Copp": round(s / 2, 4),
            "mean_min_lnd_Buium_delta2": round(mean_b, 4),
            "mean_min_lnd_Copp_r1": round(mean_c, 4),
            "ln_P_odd": round(math.log(P_odd), 4),
        })

    return {
        "P1a_holds": all_p1a,
        "index_table": index_rows,
        "P1b_holds": p1b_holds,
        "E1b_summary": {
            "W_9": round(W_9, 4),
            "b_1b": round(b_1b, 4),
            "johnson_sqrt_bW": round(johnson_9, 4),
            "min_decoded_ZY": round(min(min_ZY), 4),
            "mean_decoded_ZY": round(mean_ZY, 4),
            "min_decoded_IntZ_Istar": round(min(min_Istar), 4),
            "mean_decoded_IntZ_Istar": round(mean_Istar, 4),
            "min_decoded_IntZ_J": round(min(min_J), 4),
            "finite_hits_IntZ_J_out_of_30": len(finite_J),
            "mean_decoded_IntZ_J_when_finite": round(mean_J, 4),
        },
        "P1c_holds": p1c_holds,
        "E1c_linear_buium_vs_coppersmith": lin_comp,
    }


# ============================================================================
# 3. Experiment E2: Local q-Derivation w=0 Redundancy & CRT-Gluing Height Blowup
# ============================================================================

def run_E2() -> dict:
    rng = random.Random(20261121)
    val_checks = 0
    val_exact = 0
    for q in (2, 3, 5, 7):
        for m in (1, 2, 3, 4):
            for _ in range(60):
                u_cof = rng.randint(1, 50)
                while u_cof % q == 0:
                    u_cof += 1
                z = (q**m) * u_cof
                for ord_u in range(m + 1):
                    val_checks += 1
                    val = delta_p_iter(z, q, ord_u)
                    if v_p(val, q) == m - ord_u:
                        val_exact += 1

    ps = primes_up_to(19)  # 2,3,5,7,11,13,17,19 -> M = 9,699,690
    M = math.prod(ps)
    H = 256
    lifts_abs = []
    for j in range(2, H + 1):
        res = [delta_p(j, q) % q for q in ps]
        y, mod = crt_list(res, ps)
        assert mod == M
        if y > M // 2:
            y -= M
        lifts_abs.append(abs(y))

    min_lift = min(lifts_abs)
    med_ratio = statistics.median(lifts_abs) / M
    p2_holds = (val_exact == val_checks) and (min_lift > H) and (0.15 <= med_ratio <= 0.35)
    return {
        "P2_holds": p2_holds,
        "valuation_checks": val_checks,
        "valuation_exact_matches": val_exact,
        "M_modulus": M,
        "H_window": H,
        "min_abs_Y1_for_j_ge_2": min_lift,
        "median_abs_Y1_over_M": round(med_ratio, 4),
    }


# ============================================================================
# 4. Experiment E3: Galois-Folded O_K Over K = Q(zeta_5) vs Smooth-Norm Split-Prime Failure
# ============================================================================

# Arithmetic in Z[zeta_5] with basis (1, z, z^2, z^3), z^4 = -1 - z - z^2 - z^3
def z5_add(a: list[int], b: list[int]) -> list[int]:
    return [a[i] + b[i] for i in range(4)]


def z5_mul(a: list[int], b: list[int]) -> list[int]:
    raw = [0] * 7
    for i in range(4):
        for j in range(4):
            raw[i + j] += a[i] * b[j]
    # Reduce deg 6, 5, 4 using z^4 = -1 - z - z^2 - z^3, z^5 = 1, z^6 = z
    res = raw[:4]
    # raw[6] * z^6 = raw[6] * z
    res[1] += raw[6]
    # raw[5] * z^5 = raw[5] * 1
    res[0] += raw[5]
    # raw[4] * z^4 = -raw[4] * (1 + z + z^2 + z^3)
    for i in range(4):
        res[i] -= raw[4]
    return res


def z5_sigma(a: list[int], power: int = 1) -> list[int]:
    """Galois generator sigma(z) = z^2 on Q(zeta_5)."""
    cur = a[:]
    for _ in range(power % 4):
        # 1 -> 1; z -> z^2; z^2 -> z^4 = (-1,-1,-1,-1); z^3 -> z^6 = z
        c0, c1, c2, c3 = cur
        cur = [c0 - c2, c3 - c2, c1 - c2, -c2]
    return cur


def z5_mul_matrix(a: list[int]) -> list[list[int]]:
    """4x4 matrix M (column-action: (a*x)_row = sum_col M[row][col] * x[col])."""
    cols = []
    for col in range(4):
        basis_vec = [1 if i == col else 0 for i in range(4)]
        cols.append(z5_mul(a, basis_vec))
    return [[cols[col][row] for col in range(4)] for row in range(4)]


def rat_rank_and_det(mat: list[list[int]]) -> tuple[int, Fraction]:
    n = len(mat)
    A = [[Fraction(x) for x in row] for row in mat]
    rank = 0
    det = Fraction(1)
    for col in range(n):
        piv = -1
        for r in range(rank, n):
            if A[r][col] != 0:
                piv = r
                break
        if piv == -1:
            det = Fraction(0)
            continue
        if piv != rank:
            A[rank], A[piv] = A[piv], A[rank]
            det = -det
        diag = A[rank][col]
        det *= diag
        A[rank] = [x / diag for x in A[rank]]
        for r in range(n):
            if r != rank and A[r][col] != 0:
                f = A[r][col]
                A[r] = [A[r][k] - f * A[rank][k] for k in range(n)]
        rank += 1
    if rank < n:
        det = Fraction(0)
    return rank, det


def z5_norm(a: list[int]) -> int:
    _, d = rat_rank_and_det(z5_mul_matrix(a))
    assert d.denominator == 1
    return int(d.numerator)


def is_B_smooth(n: int, primes: list[int]) -> bool:
    n = abs(n)
    if n == 0:
        return False
    for p in primes:
        while n % p == 0:
            n //= p
        if n == 1:
            return True
    return n == 1


def run_E3() -> dict:
    rng = random.Random(20261122)
    # 11 rational primes >= 103
    all_ps = primes_up_to(200)
    qs = [q for q in all_ps if q >= 103][:11]  # 103..157
    P = math.prod(qs)
    W = 4.0 * sum(math.log(q) for q in qs)
    B_c = 60000
    b = 4.0 * math.log(2 * B_c + 1)
    R_rate = b / W
    johnson_frac = math.sqrt(R_rate)
    S_scale = 4 * B_c  # 240,000

    s3_success = 0
    s3_unique_from_basis = 0
    s1_success = 0
    max_ker_dim = 0
    measured_agr_fracs = []

    for trial in range(30):
        # Random message alpha in [-B_c, B_c]^4 not in Q
        while True:
            alpha = [rng.randint(-B_c, B_c) for _ in range(4)]
            if len({alpha[1], alpha[2], alpha[3]}) > 1:
                break
        # Pick 5 agreeing primes with A/W < johnson_frac
        while True:
            agree_qs = rng.sample(qs, 5)
            agr_frac = (4.0 * sum(math.log(q) for q in agree_qs)) / W
            if agr_frac < johnson_frac - 0.003:
                break
        measured_agr_fracs.append(agr_frac)

        # Build received word w in (Z/P Z)^4: equals alpha mod q for q in agree_qs, random mod q otherwise
        w_coords = []
        residues_per_coord = [[] for _ in range(4)]
        for q in qs:
            if q in agree_qs:
                for i in range(4):
                    residues_per_coord[i].append(alpha[i] % q)
            else:
                for i in range(4):
                    residues_per_coord[i].append(rng.randrange(q))
        for i in range(4):
            val, _ = crt_list(residues_per_coord[i], qs)
            w_coords.append(val)

        # --- s = 3 Galois-folded lattice (rank 16 over Z): (c0, c1, c2, c3) in O_K^4 ---
        # Congruence: c0 + c1*w + c2*sigma(w) + c3*sigma^2(w) == 0 (mod P O_K) -> 4 integer congruences mod P!
        w0 = w_coords
        w1 = z5_sigma(w_coords, 1)
        w2 = z5_sigma(w_coords, 2)
        M0 = z5_mul_matrix(w0)
        M1 = z5_mul_matrix(w1)
        M2 = z5_mul_matrix(w2)

        d16 = 16
        B16 = [[1 if i == j else 0 for j in range(d16)] for i in range(d16)]
        for row_idx in range(4):
            coeffs = [0] * 16
            coeffs[row_idx] = 1
            for col in range(4):
                coeffs[4 + col] = M0[row_idx][col] % P
                coeffs[8 + col] = M1[row_idx][col] % P
                coeffs[12 + col] = M2[row_idx][col] % P
            B16, _ = add_congruence_to_basis(B16, coeffs, P)
        B16 = hnf_reduce_rows(B16)
        scales16 = [1] * 4 + [S_scale] * 12
        B16_sc = [[row[k] * scales16[k] for k in range(16)] for row in B16]
        red16 = exact_lll(B16_sc)

        # Check shortest vector
        c_vec = [red16[0][k] // scales16[k] for k in range(16)]
        c0, c1, c2, c3 = c_vec[0:4], c_vec[4:8], c_vec[8:12], c_vec[12:16]
        a0, a1, a2 = alpha, z5_sigma(alpha, 1), z5_sigma(alpha, 2)
        E_val = z5_add(z5_add(c0, z5_mul(c1, a0)), z5_add(z5_mul(c2, a1), z5_mul(c3, a2)))

        # Compute Q-linear matrix of L(x) = c1*x + c2*sigma(x) + c3*sigma^2(x)
        L_cols = []
        for col in range(4):
            e_vec = [1 if i == col else 0 for i in range(4)]
            img = z5_add(
                z5_mul(c1, e_vec),
                z5_add(z5_mul(c2, z5_sigma(e_vec, 1)), z5_mul(c3, z5_sigma(e_vec, 2))),
            )
            L_cols.append(img)
        L_mat = [[L_cols[col][row] for col in range(4)] for row in range(4)]
        rk, _ = rat_rank_and_det(L_mat)
        ker_dim = 4 - rk
        if E_val == [0, 0, 0, 0] and (any(c1) or any(c2) or any(c3)) and ker_dim <= 2:
            s3_success += 1
            max_ker_dim = max(max_ker_dim, ker_dim)

        # Also check if intersecting all basis vectors with small norm uniquely determines alpha (rank 4)
        stacked_rows = []
        for b_idx in range(4):
            cv = [red16[b_idx][k] // scales16[k] for k in range(16)]
            cc0, cc1, cc2, cc3 = cv[0:4], cv[4:8], cv[8:12], cv[12:16]
            ev = z5_add(z5_add(cc0, z5_mul(cc1, a0)), z5_add(z5_mul(cc2, a1), z5_mul(cc3, a2)))
            if ev == [0, 0, 0, 0]:
                cols_b = []
                for col in range(4):
                    e_vec = [1 if i == col else 0 for i in range(4)]
                    img = z5_add(
                        z5_mul(cc1, e_vec),
                        z5_add(z5_mul(cc2, z5_sigma(e_vec, 1)), z5_mul(cc3, z5_sigma(e_vec, 2))),
                    )
                    cols_b.append(img)
                for row in range(4):
                    stacked_rows.append([cols_b[col][row] for col in range(4)])
        if len(stacked_rows) >= 4:
            # Compute rank of stacked_rows (m_rows x 4)
            A_st = [[Fraction(x) for x in r_st] for r_st in stacked_rows]
            rk_st = 0
            for col in range(4):
                piv = next((r_i for r_i in range(rk_st, len(A_st)) if A_st[r_i][col] != 0), -1)
                if piv != -1:
                    A_st[rk_st], A_st[piv] = A_st[piv], A_st[rk_st]
                    dg = A_st[rk_st][col]
                    A_st[rk_st] = [x / dg for x in A_st[rk_st]]
                    for r_i in range(len(A_st)):
                        if r_i != rk_st and A_st[r_i][col] != 0:
                            f = A_st[r_i][col]
                            A_st[r_i] = [A_st[r_i][k] - f * A_st[rk_st][k] for k in range(4)]
                    rk_st += 1
            if rk_st == 4:
                s3_unique_from_basis += 1

        # --- s = 1 baseline (rank 8 over Z): c0 + c1*w == 0 (mod P O_K) ---
        B8 = [[1 if i == j else 0 for j in range(8)] for i in range(8)]
        for row_idx in range(4):
            coeffs = [0] * 8
            coeffs[row_idx] = 1
            for col in range(4):
                coeffs[4 + col] = M0[row_idx][col] % P
            B8, _ = add_congruence_to_basis(B8, coeffs, P)
        B8 = hnf_reduce_rows(B8)
        scales8 = [1] * 4 + [S_scale] * 4
        B8_sc = [[row[k] * scales8[k] for k in range(8)] for row in B8]
        red8 = exact_lll(B8_sc)
        c8 = [red8[0][k] // scales8[k] for k in range(8)]
        E8 = z5_add(c8[0:4], z5_mul(c8[4:8], alpha))
        if E8 == [0, 0, 0, 0] and any(c8[4:8]):
            s1_success += 1

    # --- E3b: Primitive 163-smooth elements in [-6, 6]^4 \ {0} ---
    ps_163 = primes_up_to(163)
    prim_smooth_count = 0
    nonzero_orbit_count = 0
    for a0 in range(-3, 4):  # Wait: check [-6, 6]^4 = 13^4 = 28,561 elements! Fast? Let's check speed!
        pass
    # Note: 28,561 integer 4x4 determinants in pure Python Fraction might take ~8s; optimize norm via complex embeddings rounded to int!
    # Exact algebraic norm in Q(zeta_5): prod_{j=1..4} sum_{k=0..3} a_k exp(2*pi*i*j*k/5)!
    import cmath
    roots = [cmath.exp(2j * math.pi * j / 5.0) for j in range(1, 5)]
    for a0 in range(-6, 7):
        for a1 in range(-6, 7):
            for a2 in range(-6, 7):
                for a3 in range(-6, 7):
                    if a0 == 0 and a1 == 0 and a2 == 0 and a3 == 0:
                        continue
                    g = math.gcd(math.gcd(abs(a0), abs(a1)), math.gcd(abs(a2), abs(a3)))
                    if g != 1:
                        continue
                    prod_c = 1.0 + 0.0j
                    for zc in roots:
                        prod_c *= a0 + a1 * zc + a2 * (zc * zc) + a3 * (zc * zc * zc)
                    nrm = round(prod_c.real)
                    if is_B_smooth(nrm, ps_163):
                        prim_smooth_count += 1
                        # q O_K | alpha iff q | a0, a1, a2, a3 (since 1, z, z^2, z^3 is an integral basis of O_K!)
                        if g > 1:
                            nonzero_orbit_count += 1

    p3a_holds = (s3_success >= 24) and (s1_success == 0)
    p3b_holds = (prim_smooth_count > 0) and (nonzero_orbit_count == 0)
    return {
        "P3a_holds": p3a_holds,
        "P3b_holds": p3b_holds,
        "E3a_galois_folded": {
            "field": "Q(zeta_5), m=4, disc=125",
            "N_primes": 11,
            "W_nats": round(W, 4),
            "b_nats": round(b, 4),
            "rate_R": round(R_rate, 4),
            "johnson_fraction_sqrt_R": round(johnson_frac, 4),
            "max_trial_agreement_fraction": round(max(measured_agr_fracs), 4),
            "mean_trial_agreement_fraction": round(statistics.mean(measured_agr_fracs), 4),
            "s3_exact_annihilator_hits_out_of_30": s3_success,
            "s3_unique_point_from_basis_out_of_30": s3_unique_from_basis,
            "s3_max_kernel_dim": max_ker_dim,
            "s1_hits_out_of_30": s1_success,
        },
        "E3b_smooth_norm_orbit": {
            "box": "[-6, 6]^4 primitive",
            "primitive_163_smooth_elements": prim_smooth_count,
            "elements_with_A_orbit_gt_0": nonzero_orbit_count,
        },
    }


# ============================================================================
# 5. Experiment E4: Toy M2 Smooth-Search Wall on Discriminant -4N Quadratic Norms
# ============================================================================

def run_E4() -> dict:
    rng = random.Random(20261123)
    ps_47 = primes_up_to(47)  # 15 primes <= 47
    H = 512
    # Generate 28-bit semiprimes N = p*q where k0^2 + N is 47-smooth for some k0 in [1, H]
    small_ps = [p for p in primes_up_to(32000) if p > 12000]
    instances = []
    while len(instances) < 40:
        p = rng.choice(small_ps)
        q = rng.choice(small_ps)
        if p == q:
            continue
        N = p * q
        if N.bit_length() != 28:
            continue
        # Find if [1, H] has a 47-smooth value of k^2 + N
        for k0 in range(1, H + 1):
            val = k0 * k0 + N
            if is_B_smooth(val, ps_47):
                instances.append((N, k0, val))
                break

    # Full factor-base prime-power modulus P = prod_{q <= 47} q^{e_q} with q^{e_q} <= H^2 + 2^28
    X_max = H * H + (1 << 28)
    pp = []
    P_full = 1
    for q in ps_47:
        e = int(math.log(X_max) / math.log(q))
        pp.append((q, e))
        P_full *= q**e
    P_odd = math.prod(q**e for q, e in pp if q != 2)

    L, r = 4, 2
    hits_ZY = 0
    hits_IntZ = 0
    hits_Buium_lin = 0

    scales_k_ZY = [H**i for i in range(L + 1)]
    scales_k_IntZ = [max(1, binom_int(H, i)) for i in range(L + 1)]
    Q_max = X_max
    scales_Buium_Q = [1, Q_max, max(1, abs(delta_p(Q_max, 2)))]

    for N, k0, val0 in instances:
        # 1. Z[k] lattice of degree <= 4 vanishing mod P_full^2 on Q(k) = k^2 + N == 0 (mod P_full):
        #    Enforce P_full^r | F(k) whenever k^2 + N == 0 (mod P_full)? Note: N may not be -square mod all q!
        #    Use value-variable Z = Q(k, 1) in [N, N + H^2], shift Z' = Z - N = k^2 in [1, H^2]!
        #    Wait: testing directly on k in [1, H] where Q(k) = k^2 + N:
        #    Build lattice of F(k) in Z[k]_{<=4} (resp. Int(Z)_{<=4}) such that for each q^e || P_full and root r_q of r_q^2 + N == 0 mod q^f (if any),
        #    even simpler and MORE generous to the decoder: evaluate on Y = k^2 in [1, H^2] where Y + N = val0 | P_full!
        #    Since val0 | P_full, Y0 = k0^2 satisfies Y0 == -N (mod val0) with FULL agreement ln(val0) ~ 19.3 nats!
        H_sq = H * H
        w = (-N) % P_full
        sc_zy = [H_sq**i for i in range(L + 1)]
        sc_iz = [max(1, binom_int(H_sq, i)) for i in range(L + 1)]

        B_zy = build_boneh_ZY_lattice(P_full, w, r, L)
        red_zy = exact_lll([[row[i] * sc_zy[i] for i in range(L + 1)] for row in B_zy])
        c_zy = [red_zy[0][i] // sc_zy[i] for i in range(L + 1)]
        if any(c_zy) and eval_ZY_poly(c_zy, k0 * k0) == 0:
            hits_ZY += 1

        B_iz, _ = build_intZ_J_lattice(P_full, w, r, L)
        red_iz = exact_lll([[row[i] * sc_iz[i] for i in range(L + 1)] for row in B_iz])
        b_iz = [red_iz[0][i] // sc_iz[i] for i in range(L + 1)]
        if any(b_iz) and eval_intZ_poly(b_iz, k0 * k0) == 0:
            hits_IntZ += 1

        # 3. Linear Buium jet on Y = k^2: c0 + c1*Y + c2*delta_2(Y) == 0 (mod P_odd)
        w_odd = (-N) % P_odd
        sc_b = [1, H_sq, max(1, abs(delta_p(H_sq, 2)))]
        B_b0 = [[1 if i == j else 0 for j in range(3)] for i in range(3)]
        coeffs_b = [1, w_odd, delta_p_mod(w_odd, 2, 1, P_odd)]
        B_b, _ = add_congruence_to_basis(B_b0, coeffs_b, P_odd)
        red_b = exact_lll([[row[i] * sc_b[i] for i in range(3)] for row in B_b])
        c_b = [red_b[0][i] // sc_b[i] for i in range(3)]
        if any(c_b) and (c_b[0] + c_b[1] * (k0 * k0) + c_b[2] * delta_p(k0 * k0, 2) == 0):
            hits_Buium_lin += 1

    p4_holds = (hits_ZY == 0) and (hits_IntZ == 0) and (hits_Buium_lin == 0)
    return {
        "P4_holds": p4_holds,
        "instances": 40,
        "ln_P_full": round(math.log(P_full), 4),
        "ln_H_sq": round(math.log(H * H), 4),
        "johnson_threshold_nats": round(math.sqrt(math.log(P_full) * math.log(H * H)), 4),
        "mean_smooth_witness_ln_val": round(statistics.mean(math.log(v) for _, _, v in instances), 4),
        "hits_Boneh_ZY": hits_ZY,
        "hits_IntZ_Polya": hits_IntZ,
        "hits_Linear_Buium_delta2": hits_Buium_lin,
    }


def main() -> None:
    t0 = time.time()
    r1 = run_E1()
    r2 = run_E2()
    r3 = run_E3()
    r4 = run_E4()
    elapsed = time.time() - t0

    summary = {
        "elapsed_seconds": round(elapsed, 3),
        "scorecard": {
            "P1a_exact_IntZ_index": r1["P1a_holds"],
            "P1b_IntZ_vs_Boneh_LLL_radius": r1["P1b_holds"],
            "P1c_linear_Buium_degradation": r1["P1c_holds"],
            "P2_local_q_redundancy_and_CRT_blowup": r2["P2_holds"],
            "P3a_Galois_folded_OK_beyond_Johnson": r3["P3a_holds"],
            "P3b_primitive_smooth_norm_zero_orbit": r3["P3b_holds"],
            "P4_toy_M2_smooth_search_wall": r4["P4_holds"],
        },
        "E1": r1,
        "E2": r2,
        "E3": r3,
        "E4": r4,
    }
    out_path = OUT_DIR / "results.json"
    out_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
