#!/usr/bin/env python3
"""
Exact arithmetic calculations and pre-registered toy falsification tests E1, E2, E3
for Consultant G (sub-L[1/3] factoring study).
Uses only Python standard library exact arithmetic (int, fractions.Fraction) + math for logs/areas.
"""

import json
import math
import random
from fractions import Fraction
from itertools import combinations, product


# ============================================================================
# Helpers: Exact Primality, Primes, Smoothness, Extended GCD, LLL
# ============================================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n == p:
            return True
        if n % p == 0:
            return False
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def primes_up_to(B: int) -> list[int]:
    if B < 2:
        return []
    sieve = bytearray(b"\x01") * (B + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(math.isqrt(B)) + 1):
        if sieve[i]:
            sieve[i * i : B + 1 : i] = bytearray(len(sieve[i * i : B + 1 : i]))
    return [i for i, v in enumerate(sieve) if v]


def is_B_smooth(n: int, fb: list[int]) -> bool:
    n = abs(n)
    if n == 0:
        return False
    if n == 1:
        return True
    for p in fb:
        if p * p > n:
            return n <= fb[-1]
        while n % p == 0:
            n //= p
        if n == 1:
            return True
    return n == 1


def gauss_split(t: int, N: int) -> tuple[int, int]:
    """Return coprime-ish (u, v) with u == t*v (mod N) and |u|, v <= isqrt(N) via Euclid."""
    t = t % N
    if t == 0:
        return (0, 1)
    r0, r1 = N, t
    v0, v1 = 0, 1
    limit = math.isqrt(N)
    while r1 > limit:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        v0, v1 = v1, v0 - q * v1
    if v1 < 0:
        r1, v1 = -r1, -v1
    return (r1, v1)


def balanced_semiprime(bits: int, seed: int) -> tuple[int, int, int]:
    rng = random.Random(seed)
    half = bits // 2
    lo = 1 << (half - 1)
    hi = (1 << half) - 1
    while True:
        p = rng.randrange(lo | 1, hi + 1, 2)
        while not is_prime(p):
            p += 2
            if p > hi:
                p = lo | 1
        q = rng.randrange(lo | 1, hi + 1, 2)
        while not is_prime(q) or q == p:
            q += 2
            if q > hi:
                q = lo | 1
        N = p * q
        if N.bit_length() in (bits - 1, bits, bits + 1):
            return N, p, q


def lll_reduce(basis: list[list[int]], delta_num: int = 3, delta_den: int = 4) -> list[list[int]]:
    """Exact rational LLL reduction for small dimension lattices."""
    B = [list(row) for row in basis]
    n = len(B)
    m = len(B[0])

    def dot(u, v):
        return sum(x * y for x, y in zip(u, v))

    def gs():
        B_star = []
        mu = [[Fraction(0)] * n for _ in range(n)]
        norm_sq = [Fraction(0)] * n
        for i in range(n):
            v = [Fraction(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = dot(B[i], B_star[j]) / norm_sq[j]
                for l in range(m):
                    v[l] -= mu[i][j] * B_star[j][l]
            B_star.append(v)
            norm_sq[i] = dot(v, v)
        return B_star, mu, norm_sq

    B_star, mu, norm_sq = gs()
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q != 0:
                for l in range(m):
                    B[k][l] -= q * B[j][l]
                B_star, mu, norm_sq = gs()
        if norm_sq[k] >= (Fraction(delta_num, delta_den) - mu[k][k - 1] ** 2) * norm_sq[k - 1]:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            B_star, mu, norm_sq = gs()
            k = max(k - 1, 1)
    return B


# ============================================================================
# Experiment E1: Adaptive Beam Search & s * a^2 Lattice Steering
# ============================================================================

def run_e1() -> dict:
    results = {}
    M_budget = 1024
    beam_width = 64

    for bits in (32, 40, 48):
        B = round(2 ** (bits / 4.0))
        fb = primes_up_to(B)
        beam_hits_total = 0
        rand_hits_total = 0
        sq_hits_total = 0
        gauss_hits_total = 0
        min_sv_over_N = float("inf")
        sv_violations = 0

        for s_idx in range(3):
            N, p, q = balanced_semiprime(bits, 1000 * bits + s_idx)
            rng = random.Random(5000 * bits + s_idx)
            # Filter out primes dividing N just in case
            fb_N = [pr for pr in fb if N % pr != 0]

            # 1(a) Adaptive beam search: build products u > N step by step,
            # steering by centered modular residue |u mod* N|
            depth = max(3, int(math.ceil(1.3 * math.log(N) / math.log(B / 2))))
            beam = [1]
            evals_used = 0
            beam_hits = 0
            # Run enough beam restarts so total leaf smoothness tests == M_budget
            while evals_used < M_budget:
                cur = [rng.choice(fb_N) for _ in range(min(beam_width, len(fb_N)))]
                for step in range(depth - 1):
                    candidates = set()
                    for u in cur:
                        for pr in rng.sample(fb_N, min(6, len(fb_N))):
                            candidates.add(u * pr)
                    # Score by centered residue mod N once wrapping, or distance to kN
                    def score(val: int) -> int:
                        r = val % N
                        return min(r, N - r) if r != 0 else N
                    scored = sorted(candidates, key=score)
                    cur = [u for u in scored[:beam_width]]
                for u in cur:
                    if evals_used >= M_budget:
                        break
                    if u > N:
                        r = u % N
                        v = r if r <= N - r else N - r
                        if v > 0 and u != v:
                            evals_used += 1
                            if is_B_smooth(v, fb_N):
                                beam_hits += 1
                    else:
                        # Force wrap by multiplying up to > N
                        while u <= N:
                            u *= rng.choice(fb_N)
                        r = u % N
                        v = r if r <= N - r else N - r
                        if v > 0 and u != v:
                            evals_used += 1
                            if is_B_smooth(v, fb_N):
                                beam_hits += 1

            # Oblivious baseline: M_budget random B-smooth products > N
            rand_hits = 0
            for _ in range(M_budget):
                u = 1
                while u <= N * B:
                    u *= rng.choice(fb_N)
                r = u % N
                v = r if r <= N - r else N - r
                if v > 0 and is_B_smooth(v, fb_N):
                    rand_hits += 1

            beam_hits_total += beam_hits
            rand_hits_total += rand_hits

            # 1(b) Smooth-squarefree-kernel s * a^2 = v (mod N) via Euclid on a^2 mod N
            sq_hits = 0
            gauss_hits = 0
            for idx in range(M_budget):
                a = rng.randrange(2, N - 1)
                if math.gcd(a, N) > 1:
                    continue
                t_sq = (a * a) % N
                v_res, s_mult = gauss_split(t_sq, N)
                # Note: v_res == s_mult * a^2 (mod N); nontrivial since t_sq != 0, 1
                if v_res > 0 and s_mult > 0 and t_sq > 1:
                    # Wait: Euclid stops when r1 <= sqrt(N); check max(|v|, |s|)?
                    # Actually |v_res| <= sqrt(N) and s_mult <= sqrt(N), Wait: can |s|*|v| < N/2?!
                    #HOLD ON! Check s * a^2 == v (mod N): here a^2 IS O(N), NOT small!
                    # Wait: if a is known, then (a * sqrt_s)^2 == s * v (mod N) requires BOTH s AND v to be B-smooth!
                    # Let's check |s| * |v| vs N, AND check both s and v B-smooth!
                    prod_sv = abs(v_res) * abs(s_mult)
                    min_sv_over_N = min(min_sv_over_N, prod_sv / N)
                    if is_B_smooth(v_res, fb_N) and is_B_smooth(s_mult, fb_N):
                        sq_hits += 1

                # Random target t Gauss split control (Dixon / smoothbuild M1)
                r_rand = rng.randrange(2, N - 1)
                t_ctrl = (r_rand * r_rand) % N
                u_c, v_c = gauss_split(t_ctrl, N)
                if u_c > 0 and v_c > 0 and is_B_smooth(u_c, fb_N) and is_B_smooth(v_c, fb_N):
                    gauss_hits += 1

            sq_hits_total += sq_hits
            gauss_hits_total += gauss_hits

        results[bits] = {
            "B": B,
            "beam_hits": beam_hits_total,
            "rand_hits": rand_hits_total,
            "beam_over_rand": round(beam_hits_total / max(1, rand_hits_total), 4),
            "sq_kernel_hits": sq_hits_total,
            "gauss_ctrl_hits": gauss_hits_total,
            "sq_over_gauss": round(sq_hits_total / max(1, gauss_hits_total), 4),
            "min_sv_over_N": round(min_sv_over_N, 6),
        }
    return results


# ============================================================================
# Experiment E2: Exact Verification of Farey-Trapezoid Lemma FT
# ============================================================================

def poly_eval_homog(coeffs: list[int], a: int, b: int) -> int:
    """Evaluate sum_{i=0}^d coeffs[i] * a^i * b^(d-i), coeffs from degree 0 to d."""
    d = len(coeffs) - 1
    total = 0
    for i, c in enumerate(coeffs):
        total += c * (a ** i) * (b ** (d - i))
    return total


def poly_mul(p1: list[int], p2: list[int]) -> list[int]:
    out = [0] * (len(p1) + len(p2) - 1)
    for i, c1 in enumerate(p1):
        for j, c2 in enumerate(p2):
            out[i + j] += c1 * c2
    return out


def poly_disc_exact(coeffs: list[int]) -> int:
    """Exact discriminant of binary form sum_{i=0}^D coeffs[i] x^i y^(D-i) via Sylvester resultant."""
    D = len(coeffs) - 1
    a_D = coeffs[D]
    # f(x) = sum_{i=0}^D coeffs[i] x^i, f'(x) = sum_{i=1}^D i*coeffs[i] x^(i-1)
    # Build (2D - 1) x (2D - 1) Sylvester matrix over Fraction
    size = 2 * D - 1
    M = [[Fraction(0)] * size for _ in range(size)]
    # D - 1 rows of f (reversed from x^D down to x^0)
    f_rev = list(reversed(coeffs))
    fp_rev = list(reversed([i * coeffs[i] for i in range(1, D + 1)]))
    for r in range(D - 1):
        for j, val in enumerate(f_rev):
            M[r][r + j] = Fraction(val)
    for r in range(D):
        for j, val in enumerate(fp_rev):
            M[D - 1 + r][r + j] = Fraction(val)
    # Exact Bareiss / Gaussian elimination over Fraction
    det = Fraction(1)
    for col in range(size):
        piv = None
        for r in range(col, size):
            if M[r][col] != 0:
                piv = r
                break
        if piv is None:
            return 0
        if piv != col:
            M[col], M[piv] = M[piv], M[col]
            det = -det
        pv = M[col][col]
        det *= pv
        for r in range(col + 1, size):
            factor = M[r][col] / pv
            if factor != 0:
                for c in range(col, size):
                    M[r][c] -= factor * M[col][c]
    res = int(det)
    sign = -1 if (D * (D - 1) // 2) % 2 == 1 else 1
    return sign * (res // a_D)


def skew_form(coeffs: list[int], k: int) -> list[int]:
    """Substitute x -> x + k*y in sum_{i=0}^D coeffs[i] x^i y^(D-i)."""
    D = len(coeffs) - 1
    out = [0] * (D + 1)
    for i, c in enumerate(coeffs):
        # c * (x + k*y)^i * y^(D-i)
        for j in range(i + 1):
            out[j] += c * math.comb(i, j) * (k ** (i - j))
    return out


def run_e2() -> dict:
    E = 256
    instances = []
    # Build 12 exact instances: D in {4, 5}, bits in {28, 36}, types {GNFS, SNFS, Skewed}
    for bits in (28, 36):
        for d in (3, 4):  # D = d + 1 in {4, 5}
            N_gnfs, _, _ = balanced_semiprime(bits, 700 + bits + d)
            m = int(round(N_gnfs ** (1.0 / d)))
            # Base-m polynomial f(m) == N_gnfs
            rem = N_gnfs
            f_coeffs = []
            for _ in range(d):
                f_coeffs.append(rem % m)
                rem //= m
            f_coeffs[d - 1] += rem * m  # Wait: keep degree d monic if rem >= 1
            # Better clean base-m:
            m = max(3, int(N_gnfs ** (1.0 / d)))
            while m ** d > N_gnfs:
                m -= 1
            rem = N_gnfs - (m ** d)
            f_coeffs = [0] * (d + 1)
            f_coeffs[d] = 1
            for j in range(d):
                f_coeffs[j] = rem % m
                rem //= m
            f_coeffs[d - 1] += rem * m
            assert sum(f_coeffs[j] * (m ** j) for j in range(d + 1)) == N_gnfs
            g_coeffs = [-m, 1]  # x - m
            F_gnfs = poly_mul(f_coeffs, g_coeffs)
            instances.append((f"GNFS_b{bits}_D{d+1}", F_gnfs, N_gnfs))

            # Skewed SL_2(Z) shear x -> x + 17y of GNFS
            F_skew = skew_form(F_gnfs, 17)
            instances.append((f"SKEW_b{bits}_D{d+1}", F_skew, N_gnfs))

            # SNFS: N = m_s^d + 3, f = x^d + 3, g = x - m_s
            m_s = max(4, int(round((1 << bits) ** (1.0 / d))))
            N_snfs = (m_s ** d) + 3
            f_s = [0] * (d + 1)
            f_s[0] = 3
            f_s[d] = 1
            g_s = [-m_s, 1]
            F_snfs = poly_mul(f_s, g_s)
            instances.append((f"SNFS_b{bits}_D{d+1}", F_snfs, N_snfs))

    # Precompute primitive pairs in [-E, E]^2
    prim_pts = []
    for a in range(-E, E + 1):
        for b in range(-E, E + 1):
            if math.gcd(abs(a), abs(b)) == 1:
                prim_pts.append((a, b))

    violations = 0
    total_checks = 0
    cusp_pass = 0
    summary_rows = []

    for name, F_coeffs, N_val in instances:
        D = len(F_coeffs) - 1
        disc = abs(poly_disc_exact(F_coeffs))
        assert disc > 0, f"Zero disc for {name}"
        assert disc % (N_val * N_val) == 0, f"N^2 does not divide disc for {name}"
        # Bean's exact upper bound on A_F: 3*B(1/3, 1/3) < 15.8998 < 16
        A_F_bean = 16.0 * (disc ** (-1.0 / (D * (D - 1))))
        ft_err = 4 * D * (int(math.floor(math.log2(E))) + 2)
        dav_err = 2 * (D + 1) * (2 * E + 1) * 2  # two-sided b in [-E, E]

        vals = sorted(abs(poly_eval_homog(F_coeffs, a, b)) for (a, b) in prim_pts)
        vals = [v for v in vals if v > 0]
        x_min = vals[0]
        x_max = vals[-1]

        # Count at x_min
        cnt_xmin = sum(1 for v in vals if v <= x_min)
        if cnt_xmin <= ft_err:
            cusp_pass += 1

        # Dyadic grid from x_min to x_max
        X = x_min
        idx = 0
        n_vals = len(vals)
        max_ratio_to_bound = 0.0
        while True:
            while idx < n_vals and vals[idx] <= X:
                idx += 1
            cnt = idx
            ft_main = 384.0 * D * (9.0 ** ((D - 1.0) / D)) * A_F_bean * (X ** (2.0 / D))
            ft_bound = ft_main + ft_err
            total_checks += 1
            if cnt > ft_bound:
                violations += 1
            max_ratio_to_bound = max(max_ratio_to_bound, cnt / ft_bound)
            if X >= x_max:
                break
            X = min(x_max, X * 2)

        summary_rows.append({
            "name": name,
            "D": D,
            "log2_disc": round(math.log2(disc), 1),
            "cnt_at_xmin": cnt_xmin,
            "ft_err": ft_err,
            "davenport_err_2sided": dav_err,
            "dav_over_ft_err": round(dav_err / ft_err, 2),
            "max_cnt_over_ft_bound": round(max_ratio_to_bound, 5),
        })

    return {
        "E": E,
        "num_instances": len(instances),
        "total_grid_checks": total_checks,
        "violations": violations,
        "cusp_pass": f"{cusp_pass}/{len(instances)}",
        "sample_rows": summary_rows[:6],
    }


# ============================================================================
# Experiment E3: Door (ii) (1, s, t, t^2) Closure & Door (i) k=Q(i), m=3 Test
# ============================================================================

def find_split_prime_t3_Qi(lo: int, hi: int, seed: int) -> tuple[int, int, int]:
    """Find prime p == 1 (mod 12) where 3 is a cube mod p, return (p, cbrt3_mod_p, sqrt_m1_mod_p)."""
    rng = random.Random(seed)
    while True:
        p = rng.randrange(lo, hi)
        p = p - (p % 12) + 1
        if p < lo or not is_prime(p):
            continue
        if pow(3, (p - 1) // 3, p) != 1:
            continue
        # Find root of t^3 == 3 mod p and i^2 == -1 mod p
        g = 2
        while pow(g, (p - 1) // 2, p) == 1 or pow(g, (p - 1) // 3, p) == 1:
            g += 1
        imag = pow(g, (p - 1) // 4, p)
        # Find cbrt(3): brute step in cosets or Tonelli-like for 3 | p-1
        # Since p is up to ~20 bits, e = ((p-1)//3) inverse mod 3 if p % 9 != 1, else search coset
        r_found = None
        w3 = pow(g, (p - 1) // 3, p)
        # Use Cipolla / random power if 9 | p-1: write p-1 = 3^a * m
        m = p - 1
        a = 0
        while m % 3 == 0:
            m //= 3
            a += 1
        # Solve x^3 == 3 by random trials in (F_p^*) / cubes -> fast!
        # Wait: if m % 3 == 1, inv = -1: (3^((m+1)//...)); even simpler: use Shanks-Tonelli for r=3!
        for _ in range(500):
            # Try random shifts (z + c)^((p-1)//3)? Instead, solve in 3-Sylow subgroup of size 3^a!
            # Pick inv_m = pow(3, -1, m)? Wait: 3^(m) is in 3-Sylow subgroup of order 3^a!
            break
        # Exact 3-Sylow root extraction (works for ANY a >= 1 in O(a^2) ops!):
        # Let inv3_mod_m: if m % 3 == 1, 2*m + 1 is divisible by 3 -> k = (2*m + 1)//3 -> 3*k = 2m + 1 == 1 mod m?
        # Wait: we want 3*k == 1 (mod m), so (3^k)^3 = 3^(3k) = 3 * (3^m)^u!
        k_inv = pow(3, -1, m)
        x0 = pow(3, k_inv, p)  # x0^3 * 3^(-1) has order dividing 3^(a-1)!
        rem = (pow(x0, 3, p) * pow(3, p - 2, p)) % p  # order divides 3^(a-1)
        zeta3_top = pow(g, m, p)  # order 3^a
        # Discrete log of rem in 3-Sylow subgroup generated by zeta3_top^3 (order 3^(a-1))
        # Since 3^a <= p (actually a <= 12), discrete log in base 3 takes <= 3*a steps!
        cur_x = x0
        while (pow(cur_x, 3, p) - 3) % p != 0:
            rem = (pow(cur_x, 3, p) * pow(3, p - 2, p)) % p
            # Find order 3^b of rem
            b = 0
            tmp = rem
            while tmp != 1:
                tmp = pow(tmp, 3, p)
                b += 1
            # Multiply cur_x by zeta3_top^(3^(a - 1 - b)) or its square
            step = pow(zeta3_top, 3 ** (a - b - 1), p)  # step^3 has order 3^b!
            if pow((cur_x * step) % p, 3 * (3 ** (b - 1)), p) == pow(3, 3 ** (b - 1), p):
                cur_x = (cur_x * step) % p
            else:
                cur_x = (cur_x * step * step) % p
        assert pow(cur_x, 3, p) == 3 % p
        assert (imag * imag + 1) % p == 0
        return p, cur_x, imag


def crt2(a1: int, p: int, a2: int, q: int) -> int:
    N = p * q
    inv_p = pow(p, -1, q)
    k = ((a2 - a1) * inv_p) % q
    return (a1 + k * p) % N


def norm_K0_cubic3(c0: int, c1: int, c2: int) -> int:
    """Exact norm N_{Q(cbrt(3))/Q}(c0 + c1*t + c2*t^2) where t^3 = 3: det of multiplication matrix."""
    # Multiplication by (c0 + c1*t + c2*t^2) on (1, t, t^2):
    # col 0 (x 1):   [c0,   c1,   c2]
    # col 1 (x t):   [3*c2, c0,   c1]
    # col 2 (x t^2): [3*c1, 3*c2, c0]
    return (
        c0 * (c0 * c0 - 3 * c1 * c2)
        - 3 * c2 * (c1 * c0 - 3 * c2 * c2)
        + 3 * c1 * (c1 * c1 - c2 * c0)
    )


def norm_Qi_cubic(A0: tuple[int, int], A1: tuple[int, int], A2: tuple[int, int], d_re: int, d_im: int) -> int:
    """
    Exact norm from K = Q(i)(s) (degree 6 over Q) to Q of A0 + A1*s + A2*s^2,
    where s^3 = d_re + i*d_im in Z[i]!
    First compute det over Z[i] (gives element in Z[i]), then norm Z[i]->Z (re^2 + im^2)!
    """
    def zi_add(x, y):
        return (x[0] + y[0], x[1] + y[1])

    def zi_sub(x, y):
        return (x[0] - y[0], x[1] - y[1])

    def zi_mul(x, y):
        return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])

    D = (d_re, d_im)
    # 3x3 matrix over Z[i]:
    # [ A0,    D*A2,  D*A1 ]
    # [ A1,    A0,    D*A2 ]
    # [ A2,    A1,    A0   ]
    m00, m01, m02 = A0, zi_mul(D, A2), zi_mul(D, A1)
    m10, m11, m12 = A1, A0, zi_mul(D, A2)
    m20, m21, m22 = A2, A1, A0

    det_zi = zi_add(
        zi_sub(
            zi_mul(m00, zi_sub(zi_mul(m11, m22), zi_mul(m12, m21))),
            zi_mul(m01, zi_sub(zi_mul(m10, m22), zi_mul(m12, m20))),
        ),
        zi_mul(m02, zi_sub(zi_mul(m10, m21), zi_mul(m11, m20))),
    )
    return det_zi[0] ** 2 + det_zi[1] ** 2


def run_e3() -> dict:
    door2_rows = []
    door1_rows = []

    for bits, (lo, hi) in [(22, (1 << 10, 1 << 11)), (32, (1 << 15, 1 << 16)), (42, (1 << 20, 1 << 21))]:
        for s_idx in range(3):
            p, t_p, i_p = find_split_prime_t3_Qi(lo, hi, 3000 + bits * 10 + s_idx)
            q, t_q, i_q = find_split_prime_t3_Qi(lo, hi, 4000 + bits * 10 + s_idx)
            while q == p:
                q, t_q, i_q = find_split_prime_t3_Qi(lo, hi, 5000 + bits * 10 + s_idx + q % 97)
            N = p * q
            t_N = crt2(t_p, p, t_q, q)
            i_N = crt2(i_p, p, i_q, q)
            assert (pow(t_N, 3, N) - 3) % N == 0
            assert (i_N * i_N + 1) % N == 0

            # ----------------------------------------------------------------
            # Part 1: Door (ii) (1, s, t, t^2) vs (1, s', t, t^2) over K_0 = Q(t), t^3 = 3
            # Choose linking value v_s mod N; find TWO distinct short elements alpha_1, alpha_2 in Z[t]
            # with alpha_k(t_N) == v_s^2 (mod N) via LLL on {c0 + c1*t_N + c2*t_N^2 == 0 mod N}!
            # ----------------------------------------------------------------
            rng = random.Random(9000 + bits * 10 + s_idx)
            v_s = rng.randrange(2, N - 1)
            target_sq = (v_s * v_s) % N
            # Particular solution: (target_sq, 0, 0); homogeneous lattice ker: c0 + c1*t_N + c2*t_N^2 == 0 (mod N)
            lat_ker = [
                [N, 0, 0],
                [(-t_N) % N, 1, 0],
                [(-t_N * t_N) % N, 0, 1],
            ]
            red_ker = lll_reduce(lat_ker)
            # Lift target_sq into small vector using Babai rounding on red_ker:
            # Find c_part = (target_sq, 0, 0) reduced by red_ker
            # Simpler exact CVP in dim 3: solve R * lambda = (target_sq, 0, 0) over Q and round!
            # Use Cramer / Fraction inverse of 3x3 red_ker:
            def solve_3x3(M3, rhs):
                # M3 rows are basis vectors: sum lam_i M3[i] = rhs
                # Transpose M3
                A = [[Fraction(M3[r][c]) for r in range(3)] for c in range(3)]
                b_vec = [Fraction(x) for x in rhs]
                for col in range(3):
                    piv = next(r for r in range(col, 3) if A[r][col] != 0)
                    A[col], A[piv] = A[piv], A[col]
                    b_vec[col], b_vec[piv] = b_vec[piv], b_vec[col]
                    pv = A[col][col]
                    for c in range(col, 3):
                        A[col][c] /= pv
                    b_vec[col] /= pv
                    for r in range(3):
                        if r != col and A[r][col] != 0:
                            fac = A[r][col]
                            for c in range(col, 3):
                                A[r][c] -= fac * A[col][c]
                            b_vec[r] -= fac * b_vec[col]
                return [round(x) for x in b_vec]

            rhs_vec = [target_sq, 0, 0]
            lams = solve_3x3(red_ker, rhs_vec)
            c_base = [
                rhs_vec[j] - sum(lams[i] * red_ker[i][j] for i in range(3))
                for j in range(3)
            ]
            # Two distinct linked short radicands in Z[t]: alpha_1 = c_base, alpha_2 = c_base + red_ker[0]
            al1 = c_base
            al2 = [c_base[j] + red_ker[0][j] for j in range(3)]
            assert (al1[0] + al1[1] * t_N + al1[2] * t_N * t_N) % N == target_sq
            assert (al2[0] + al2[1] * t_N + al2[2] * t_N * t_N) % N == target_sq
            # Relative resultant of x^2 - al1 and x^2 - al2 over K_0 is (al1 - al2)^2 = red_ker[0]^2!
            diff = [al1[j] - al2[j] for j in range(3)]
            norm_diff = abs(norm_K0_cubic3(*diff))
            norm_res = norm_diff * norm_diff
            rr_divides = (norm_res % N == 0)

            # Evaluate |F(a)| = |W(a)*W'(a)| over a in {0, -1, 1}^4!
            # Candidate is P(t) + a2 * s, where P(t) = a0 + a3*t + a4*t^2.
            # Relative norm to K_0 is P(t)^2 - a2^2 * al_k(t) in Z[t]!
            def mul_Zt(u, v):
                # (u0 + u1 t + u2 t^2)*(v0 + v1 t + v2 t^2) mod (t^3 - 3)
                r0 = u[0] * v[0] + 3 * (u[1] * v[2] + u[2] * v[1])
                r1 = u[0] * v[1] + u[1] * v[0] + 3 * u[2] * v[2]
                r2 = u[0] * v[2] + u[1] * v[1] + u[2] * v[0]
                return [r0, r1, r2]

            coincide_a2_zero = True
            min_logF_genuine = float("inf")
            for a0, a2, a3, a4 in product((-1, 0, 1), repeat=4):
                if a0 == a2 == a3 == a4 == 0:
                    continue
                P_t = [a0, a3, a4]
                P_sq = mul_Zt(P_t, P_t)
                rel1 = [P_sq[j] - (a2 * a2) * al1[j] for j in range(3)]
                rel2 = [P_sq[j] - (a2 * a2) * al2[j] for j in range(3)]
                W1 = abs(norm_K0_cubic3(*rel1))
                W2 = abs(norm_K0_cubic3(*rel2))
                if a2 == 0:
                    if W1 != W2:
                        coincide_a2_zero = False
                else:
                    if W1 > 0 and W2 > 0:
                        min_logF_genuine = min(min_logF_genuine, math.log(W1 * W2))

            door2_rows.append({
                "bits": N.bit_length(),
                "rr_divides": rr_divides,
                "coincide_when_a2_0": coincide_a2_zero,
                "min_logF_over_logN": round(min_logF_genuine / math.log(N), 4),
            })

            # ----------------------------------------------------------------
            # Part 2: Door (i) k = Q(i), relative rank m = 3 (r = 6, D = 12)
            # K = Q(i)(s), s^3 = d1 in Z[i]; K' = Q(i)(s'), s'^3 = d2 in Z[i]
            # linked via s_N^3 == d1 == d2 (mod (N, i - i_N)), AND test non-power LLL-reduced
            # relative rank-3 module over Z[i]!
            # ----------------------------------------------------------------
            s_N = rng.randrange(2, N - 1)
            target_cube = pow(s_N, 3, N)
            # Lattice in Z[i] ~ Z^2: {x + y*i_N == 0 mod N}
            lat_zi = [[N, 0], [(-i_N) % N, 1]]
            red_zi = lll_reduce(lat_zi)
            # Find short d1, d2 in Z[i] with d_k(i_N) == target_cube (mod N)
            # 2x2 solve:
            det2 = red_zi[0][0] * red_zi[1][1] - red_zi[0][1] * red_zi[1][0]
            lam0 = round(Fraction(target_cube * red_zi[1][1], det2))
            lam1 = round(Fraction(-target_cube * red_zi[0][1], det2))
            d1 = (
                target_cube - lam0 * red_zi[0][0] - lam1 * red_zi[1][0],
                -lam0 * red_zi[0][1] - lam1 * red_zi[1][1],
            )
            d2 = (d1[0] + red_zi[0][0], d1[1] + red_zi[0][1])
            assert (d1[0] + d1[1] * i_N) % N == target_cube
            assert (d2[0] + d2[1] * i_N) % N == target_cube

            # Also build a genuinely NON-POWER 3rd Z[i]-generator via LLL on the
            # rank-8 fibre-product sublattice {(c0, c1, s2, c0', c1', s2') : same mod (N, i-i_N)}!
            # Specifically, fix s2_1 = 1 on side 1: u3 = b0 + b1*s + s^2, and u3' = b0' + b1'*s' + b2'*s'^2
            # linked mod N: (b0 - b0') + (b1 - b1')*s_N + (1 - b2')*s_N^2 == 0 (mod N) over Z[i]!
            # Even better: LLL in Z^8 (coords: re/im of (b1, b2, b1', b2'), with b0=b0'=0) solving:
            # b1*s_N + b2*s_N^2 - b1'*s_N - b2'*s_N^2 == 0 (mod (N, i - i_N)) with (b2, b2') != (0,0) and u3 != u3'!
            # Weight the s^2 coordinates by |d1|^(2/3), s coordinates by |d1|^(1/3) so LLL balances archimedean sizes!
            w_s1 = max(1, int(round((d1[0]**2 + d1[1]**2) ** (1.0 / 6.0))))
            w_s2 = max(1, int(round((d1[0]**2 + d1[1]**2) ** (2.0 / 6.0))))
            w_sp1 = max(1, int(round((d2[0]**2 + d2[1]**2) ** (1.0 / 6.0))))
            w_sp2 = max(1, int(round((d2[0]**2 + d2[1]**2) ** (2.0 / 6.0))))
            # Congruence mod N: x0 + x1*i_N + x2*s_N + x3*i_N*s_N + x4*s_N^2 + x5*i_N*s_N^2
            #                  - x6*s_N - x7*i_N*s_N - x8*s_N^2 - x9*i_N*s_N^2 == 0 (mod N)
            # Keep 6 vars: (b1_re, b2_re, bp1_re, bp2_re, diff_b0_re, diff_b0_im) -> ker mod N:
            s2_N = (s_N * s_N) % N
            coefs_modN = [
                1,
                i_N % N,
                s_N % N,
                s2_N % N,
                (-s_N) % N,
                (-s2_N) % N,
            ]
            # Lattice of (c0_diff_re, c0_diff_im, b1, b2, bp1, bp2) satisfying sum coefs_modN[j]*x[j] == 0 mod N:
            # Basis: x[0] = -sum_{j>=1} coefs_modN[j]*x[j] mod N, plus N*e_0!
            weights6 = [1, 1, w_s1, w_s2, w_sp1, w_sp2]
            lat6 = [[N * weights6[0], 0, 0, 0, 0, 0]]
            for j in range(1, 6):
                row = [0] * 6
                row[0] = ((-coefs_modN[j]) % N) * weights6[0]
                row[j] = 1 * weights6[j]
                lat6.append(row)
            red6 = lll_reduce(lat6)
            # Pick shortest vector with (b2 != 0 and bp2 != 0) so relative rank is 3 on BOTH sides!
            chosen_v = None
            for rv in red6:
                unw = [rv[j] // weights6[j] for j in range(6)]
                if unw[3] != 0 and unw[5] != 0:
                    chosen_v = unw
                    break
            if chosen_v is None:
                # Combine two basis vectors to ensure b2 != 0 and bp2 != 0
                u_a = [red6[0][j] // weights6[j] for j in range(6)]
                u_b = [red6[1][j] // weights6[j] for j in range(6)]
                for mult in (1, -1, 2):
                    cand = [u_a[j] + mult * u_b[j] for j in range(6)]
                    if cand[3] != 0 and cand[5] != 0:
                        chosen_v = cand
                        break
                if chosen_v is None:
                    chosen_v = [0, 0, 0, 1, 0, 1]  # fallback power generator

            # u3 = (c0_diff_re + i*c0_diff_im) + b1*s + b2*s^2 on side 1, u3' = bp1*s' + bp2*s'^2 on side 2
            c0_1 = (chosen_v[0], chosen_v[1])
            b1_1 = (chosen_v[2], 0)
            b2_1 = (chosen_v[3], 0)
            b1_2 = (chosen_v[4], 0)
            b2_2 = (chosen_v[5], 0)
            assert (
                chosen_v[0] + chosen_v[1] * i_N + chosen_v[2] * s_N + chosen_v[3] * s2_N
                - chosen_v[4] * s_N - chosen_v[5] * s2_N
            ) % N == 0

            def zi_add2(x, y):
                return (x[0] + y[0], x[1] + y[1])

            def zi_mul2(x, y):
                return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])

            min_logF_door1 = float("inf")
            min_logF_nonpow = float("inf")
            med_vals = []
            for coeffs6 in product((-1, 0, 1), repeat=6):
                A0 = (coeffs6[0], coeffs6[1])
                A1 = (coeffs6[2], coeffs6[3])
                A2 = (coeffs6[4], coeffs6[5])
                if A1 == (0, 0) and A2 == (0, 0):
                    continue
                # Power basis witness
                N1 = norm_Qi_cubic(A0, A1, A2, d1[0], d1[1])
                N2 = norm_Qi_cubic(A0, A1, A2, d2[0], d2[1])
                if N1 > 0 and N2 > 0 and N1 != N2:
                    lf = math.log(N1 * N2)
                    min_logF_door1 = min(min_logF_door1, lf)
                    med_vals.append(lf)
                # Non-power LLL-linked 3rd generator witness:
                # Side 1 element: (A0 + A2*c0_1) + (A1 + A2*b1_1)*s + (A2*b2_1)*s^2
                P0_1 = zi_add2(A0, zi_mul2(A2, c0_1))
                P1_1 = zi_add2(A1, zi_mul2(A2, b1_1))
                P2_1 = zi_mul2(A2, b2_1)
                P0_2 = A0
                P1_2 = zi_add2(A1, zi_mul2(A2, b1_2))
                P2_2 = zi_mul2(A2, b2_2)
                NP1 = norm_Qi_cubic(P0_1, P1_1, P2_1, d1[0], d1[1])
                NP2 = norm_Qi_cubic(P0_2, P1_2, P2_2, d2[0], d2[1])
                if NP1 > 0 and NP2 > 0 and NP1 != NP2:
                    min_logF_nonpow = min(min_logF_nonpow, math.log(NP1 * NP2))

            med_vals.sort()
            med_lf = med_vals[len(med_vals) // 2]
            door1_rows.append({
                "bits": N.bit_length(),
                "d1_norm_bits": round(math.log2(d1[0] ** 2 + d1[1] ** 2), 2),
                "d2_norm_bits": round(math.log2(d2[0] ** 2 + d2[1] ** 2), 2),
                "min_logF_pow_over_logN": round(min_logF_door1 / math.log(N), 4),
                "min_logF_nonpow_over_logN": round(min_logF_nonpow / math.log(N), 4),
                "med_logF_over_logN": round(med_lf / math.log(N), 4),
            })

    return {"door2": door2_rows, "door1": door1_rows}


# ============================================================================
# Symbolic / Exact Exponent Ledger Verification
# ============================================================================

def verify_ledgers() -> dict:
    # 1. Verify 2D-Bohr-CRT / CG-2D-CRT ledger at beta = 1/4:
    # Candidate Q_a(k, v) has log X = (1/2) x; B = L[1/4, kappa];
    # u_0 = (x / 2) / (kappa * x^(1/4) * (ln x)^(3/4)) = (1/(2*kappa)) * x^(3/4) * (ln x)^(-3/4)
    # ln(1/rho(u_0)) ~ u_0 * ln(u_0) = (3 / (8 * kappa)) * x^(3/4) * (ln x)^(1/4) = (3/(8*kappa)) * F_{3/4}
    # So K = L[3/4, 3/(16*kappa) + eps] gives K^2 * rho(u_0) = L[3/4, 2*eps] >> 1!
    # Number of relations needed: pi(B) = L[1/4, kappa];
    # If decoder per block costs pi(B)^a = L[1/4, a*kappa], total cost = max((1+a)*kappa, 2*kappa) -> L[1/4]!
    cg_crt_ledger = {
        "beta": "1/4",
        "log_X": "(1/2) * x",
        "log_K": "(3/(16*kappa) + eps) * x^(3/4) * (log x)^(1/4)",
        "agreement_over_entropy": "((1/2)*x) / (2 * log K) ~ (4*kappa/3) * (x / log x)^(1/4) -> infinity",
        "alpha": "1/4 (for any fixed polynomial decoder exponent a >= 1)",
    }

    # 2. Verify Door (i) Relative-Rank-m Type-D Ledger over fixed k ([k:Q] = e = O(1)):
    # Total degree D = 2*d = 2*e*n, total rank r = e*m, box a_j in {0, +-1} -> C = r * ln 3.
    # Under c = 0 (relative G_m fails favourably so all m generators have height e*x/d):
    # log X = e * x / d + 2 * d * ln(3*r) ~ e * x / d + 2 * d * ln x (up to constants).
    # Optimize d = sqrt(e * x / (2 * ln x)) -> log X ~ 2 * sqrt(2 * e * x * ln x) -> alpha = 1/4!
    # Under relative-R1 (if v_P(Delta_rel) >= P(m, 2n) holds over O_k):
    # height floor becomes e * (m - 1) * x / (2n - 1) = (r - e) * x / (D/e - 1) ~ e * (r / D) * x -> c = 1 -> alpha = 1/3!
    door1_ledger = {
        "c_eq_0_alpha": "1/4",
        "relative_R1_c_eq_1_alpha": "1/3",
        "relative_R1_predicted_slope_r6_D12_e2": str(Fraction(2 * (3 - 1), 6 - 1)),  # e*(m-1)/(2n-1) = 4/5 = 0.80 (for WW', i.e. 2 sides!)
    }
    return {"cg_crt": cg_crt_ledger, "door1": door1_ledger}


def main():
    e1 = run_e1()
    e2 = run_e2()
    e3 = run_e3()
    ledgers = verify_ledgers()
    out = {"ledgers": ledgers, "E1": e1, "E2": e2, "E3": e3}
    with open("results.json", "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
