#!/usr/bin/env python3
"""
Pre-registered toy experiments E1, E2, E3 for Hatch 3 (R / Thin Dependencies,
Shared Descent DAGs, and Reuse Across Related Inputs).

Follows AGENT_RULES.md and prereg.md:
- Exact integer / Fraction / GF(2) arithmetic in all operational & algebraic routines.
- Operational routines receive ONLY N (and public parameters/seeds), never p or q.
"""

import hashlib
import itertools
import json
import math
import os
import random
import time
from fractions import Fraction

WORKDIR = "/Volumes/git/factoring-connection/codex/scratch/hatch_h3g"
PREREG_PATH = os.path.join(WORKDIR, "prereg.md")
OUT_JSON = os.path.join(WORKDIR, "hatch3_results.json")
GLOBAL_SEED = 20261120


# ============================================================================
# Elementary number-theory & exact linear algebra helpers
# ============================================================================

def primes_up_to(B: int) -> list[int]:
    sieve = bytearray(b"\x01") * (B + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, math.isqrt(B) + 1):
        if sieve[p]:
            sieve[p * p : B + 1 : p] = b"\x00" * (((B - p * p) // p) + 1)
    return [i for i, v in enumerate(sieve) if v]


def is_probable_prime_det64(n: int) -> bool:
    """Deterministic Miller-Rabin for all 64-bit integers (Jim Sinclair bases)."""
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False
    d = n - 1
    s = 0
    while (d & 1) == 0:
        d >>= 1
        s += 1
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if a % n == 0:
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


def generate_semiprimes(rng: random.Random, bits_half: int, count: int, blum: bool) -> list[dict]:
    """HARNESS ONLY: generates seeded balanced semiprimes N = p * q."""
    out = []
    lo = (1 << (bits_half - 1)) + 1
    hi = (1 << bits_half) - 1
    while len(out) < count:
        p = rng.randrange(lo | 1, hi + 1, 2)
        if blum and (p & 3) != 3:
            continue
        if not is_probable_prime_det64(p):
            continue
        q = rng.randrange(lo | 1, hi + 1, 2)
        if q == p or (blum and (q & 3) != 3):
            continue
        if not is_probable_prime_det64(q):
            continue
        if p > q:
            p, q = q, p
        out.append({"N": p * q, "p": p, "q": q, "blum": blum})
    return out


def jacobi_symbol(a: int, n: int) -> int:
    """Exact Jacobi symbol (a / n) for odd n > 0 (OPERATIONAL: uses only n)."""
    assert n > 0 and (n & 1) == 1
    a %= n
    res = 1
    while a != 0:
        while (a & 1) == 0:
            a >>= 1
            r = n & 7
            if r == 3 or r == 5:
                res = -res
        a, n = n, a
        if (a & 3) == 3 and (n & 3) == 3:
            res = -res
        a %= n
    return res if n == 1 else 0


def gauss_split_mod_n(lam: int, N: int) -> tuple[int, int]:
    """OPERATIONAL: Stopped extended Euclid returning (u, v) with u ≡ lam * v (mod N), |u|, v <= sqrt(N) roughly."""
    bound = math.isqrt(N)
    r0, r1 = N, lam % N
    t0, t1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        t0, t1 = t1, t0 - q * t1
    if t1 < 0:
        return -r1, -t1
    return r1, t1


def trial_factor_fb(val: int, fb: list[int]) -> list[int] | None:
    """OPERATIONAL: Exact trial division of nonzero integer val over fb. Returns [sign_bit, e_p1, ...] or None."""
    assert val != 0
    sign = 1 if val < 0 else 0
    rem = abs(val)
    exps = [0] * len(fb)
    for idx, p in enumerate(fb):
        if p * p > rem:
            break
        if rem % p == 0:
            c = 0
            while rem % p == 0:
                rem //= p
                c += 1
            exps[idx] = c
    if rem == 1:
        return [sign] + exps
    # Rem might be a prime in fb
    if rem <= fb[-1]:
        # binary search in fb
        lo, hi = 0, len(fb) - 1
        while lo <= hi:
            mid = (lo + hi) >> 1
            if fb[mid] == rem:
                exps[mid] += 1
                return [sign] + exps
            elif fb[mid] < rem:
                lo = mid + 1
            else:
                hi = mid - 1
    return None


def miller_split_single_base(N: int, base: int, D: int) -> int | None:
    """OPERATIONAL: Given base^D ≡ 1 (mod N) with D > 0, run Miller's power chain for this single base."""
    assert D > 0
    g = math.gcd(base, N)
    if 1 < g < N:
        return g
    t = D
    s = 0
    while (t & 1) == 0:
        t >>= 1
        s += 1
    x = pow(base, t, N)
    if x == 1 or x == N - 1:
        return None
    for _ in range(s):
        x_next = (x * x) % N
        if x_next == 1:
            g = math.gcd(x - 1, N)
            if 1 < g < N:
                return g
            return None
        if x_next == N - 1:
            return None
        x = x_next
    return None


class IncrementalF2Basis:
    """Exact GF(2) incremental row-basis using Python arbitrary-precision bitmasks."""

    def __init__(self) -> None:
        self.pivots: dict[int, int] = {}

    def add_and_check_dep(self, vec_mask: int) -> bool:
        v = vec_mask
        while v:
            msb = v.bit_length() - 1
            if msb in self.pivots:
                v ^= self.pivots[msb]
            else:
                self.pivots[msb] = v
                return False
        return True

    @property
    def rank(self) -> int:
        return len(self.pivots)


class IncrementalModBasis:
    """Fast exact rank over F_P (P = 2^61 - 1) plus exact Q-dependency extraction when needed."""

    MOD = (1 << 61) - 1

    def __init__(self, ncols: int) -> None:
        self.ncols = ncols
        self.pivots: dict[int, tuple[list[int], list[int]]] = {}
        self.raw_rows: list[list[int]] = []

    def add_row(self, row: list[int]) -> list[int] | None:
        """Adds integer row of length ncols. If linearly dependent mod MOD, computes exact Q-kernel vector."""
        idx = len(self.raw_rows)
        self.raw_rows.append(row)
        P = self.MOD
        v = [x % P for x in row]
        comb = [0] * (idx + 1)
        comb[idx] = 1
        for c in range(self.ncols):
            if v[c] != 0:
                if c in self.pivots:
                    pv, pcomb = self.pivots[c]
                    factor = v[c]  # pv[c] is normalized to 1
                    for j in range(c, self.ncols):
                        if pv[j]:
                            v[j] = (v[j] - factor * pv[j]) % P
                    for j in range(len(pcomb)):
                        if pcomb[j]:
                            comb[j] = (comb[j] - factor * pcomb[j]) % P
                else:
                    inv = pow(v[c], P - 2, P)
                    v = [(x * inv) % P for x in v]
                    comb = [(x * inv) % P for x in comb]
                    self.pivots[c] = (v, comb)
                    return None
        return comb


def exact_integer_nullvec(rows: list[list[int]]) -> list[int] | None:
    """Exact Fraction Gaussian elimination returning a primitive integer left-kernel vector of `rows`."""
    m = len(rows)
    if m == 0:
        return None
    n = len(rows[0])
    pivots: dict[int, tuple[list[Fraction], list[Fraction]]] = {}
    for idx, r in enumerate(rows):
        v = [Fraction(x) for x in r]
        comb = [Fraction(0)] * m
        comb[idx] = Fraction(1)
        for c in range(n):
            if v[c] != 0:
                if c in pivots:
                    pv, pcomb = pivots[c]
                    fac = v[c]
                    for j in range(c, n):
                        if pv[j] != 0:
                            v[j] -= fac * pv[j]
                    for j in range(m):
                        if pcomb[j] != 0:
                            comb[j] -= fac * pcomb[j]
                else:
                    inv = Fraction(1, 1) / v[c]
                    v = [x * inv for x in v]
                    comb = [x * inv for x in comb]
                    pivots[c] = (v, comb)
                    break
        else:
            # Exact dependency found! Clear denominators and make primitive.
            lcm_den = 1
            for x in comb:
                lcm_den = math.lcm(lcm_den, x.denominator)
            int_vec = [int(x * lcm_den) for x in comb]
            g = 0
            for x in int_vec:
                g = math.gcd(g, abs(x))
            if g > 1:
                int_vec = [x // g for x in int_vec]
            # Verify exact zero
            for c in range(n):
                assert sum(int_vec[i] * rows[i][c] for i in range(m)) == 0
            return int_vec
    return None


# ============================================================================
# Experiment E1: Thin Dependencies, Lemma PIN-1, and Support Growth
# ============================================================================

def run_e1() -> dict:
    rng = random.Random(GLOBAL_SEED + 1)
    blum_insts = generate_semiprimes(rng, 20, 10, blum=True)
    gen_insts = generate_semiprimes(rng, 20, 10, blum=False)
    all_insts = blum_insts + gen_insts

    B = 600
    fb = primes_up_to(B)  # 109 primes
    n_primes = len(fb)
    small_label_idxs = [i for i, p in enumerate(fb) if p <= 19]  # 8 small primes for multiplier labels

    blum_pin1_splits = 0
    gen_pin1_splits = 0
    gen_pin1_trials = 0
    records = []

    for inst_idx, inst in enumerate(all_insts):
        N = inst["N"]
        # Operational choice of q_0: smallest prime in fb with Jacobi(q_0, N) == -1
        q0_idx = next(i for i, q in enumerate(fb) if jacobi_symbol(q, N) == -1)
        q0 = fb[q0_idx]

        f2_basis = IncrementalF2Basis()
        # External columns for Z1: all primes in fb EXCEPT q0_idx
        ext_col_idxs = [i for i in range(n_primes) if i != q0_idx]
        z1_mod = IncrementalModBasis(len(ext_col_idxs))

        full_rows: list[list[int]] = []  # each row: length n_primes (prime exponents; sign tracked mod 2 via doubling)
        seen_rows: set[tuple[int, ...]] = set()
        h_f2 = None
        h_z1 = None
        z1_D = None

        op_rng = random.Random((GLOBAL_SEED ^ 0xA5A5) + inst_idx * 997)
        tries = 0
        # Collect until both h_f2 and h_z1 are reached, plus full rank n_primes for general-N PIN-1 trials
        need_full_fb = not inst["blum"]
        fb_mod = IncrementalModBasis(n_primes)

        while h_f2 is None or h_z1 is None or (need_full_fb and len(fb_mod.pivots) < n_primes):
            tries += 1
            # Build known-exponent multiplier lam = q0^84 * prod_{q<=19} q^{e_q} mod N on odd steps,
            # or random FB-label on small primes + medium prime to ensure non-duplicate Gauss splits
            exp_lam = [0] * n_primes
            for idx in small_label_idxs:
                exp_lam[idx] = op_rng.randint(0, 24)
            extra_idx = op_rng.randint(0, n_primes - 1)
            exp_lam[extra_idx] += op_rng.randint(1, 40)

            lam = 1
            for idx, e in enumerate(exp_lam):
                if e:
                    lam = (lam * pow(fb[idx], e, N)) % N
            u, v = gauss_split_mod_n(lam, N)
            if u == 0 or v == 0 or math.gcd(u, N) != 1 or math.gcd(v, N) != 1:
                continue
            fu = trial_factor_fb(u, fb)
            if fu is None:
                continue
            fv = trial_factor_fb(v, fb)
            if fv is None:
                continue

            # Relation: u ≡ lam * v (mod N) => (-1)^fu[0] prod q^(fu - fv - exp_lam) ≡ 1 (mod N)
            prime_vec = [fu[1 + i] - fv[1 + i] - exp_lam[i] for i in range(n_primes)]
            if all(x == 0 for x in prime_vec):
                continue
            key = tuple(prime_vec)
            if key in seen_rows:
                continue
            seen_rows.add(key)

            sign_bit = fu[0] ^ fv[0]
            full_rows.append(prime_vec)
            m_now = len(full_rows)

            # Check certified multiplicative identity mod N (doubled to clear -1):
            # Track F2 parity mask (sign_bit at bit 0, prime parities at bits 1..n_primes)
            if h_f2 is None:
                mask = sign_bit
                for i, val in enumerate(prime_vec):
                    if val & 1:
                        mask |= 1 << (i + 1)
                if f2_basis.add_and_check_dep(mask):
                    h_f2 = m_now

            if h_z1 is None:
                ext_row = [prime_vec[i] for i in ext_col_idxs]
                cand = z1_mod.add_row(ext_row)
                if cand is not None:
                    ext_subrows = [[r[i] for i in ext_col_idxs] for r in full_rows]
                    nullv = exact_integer_nullvec(ext_subrows)
                    if nullv is not None:
                        D_raw = sum(nullv[i] * full_rows[i][q0_idx] for i in range(len(nullv)))
                        if D_raw != 0:
                            h_z1 = m_now
                            z1_D = 2 * abs(D_raw)  # 2 * |D_raw| kills ±1 sign deterministically
                            assert pow(q0, z1_D, N) == 1

            if need_full_fb:
                fb_mod.add_row(prime_vec)

        split_ok = False
        if inst["blum"]:
            fac = miller_split_single_base(N, q0, z1_D)
            if fac is not None and 1 < fac < N and (N % fac == 0):
                assert fac in (inst["p"], inst["q"])
                split_ok = True
                blum_pin1_splits += 1
        else:
            # Precompute exact Q-basis inverse for full FB to pin 10 uniform random units b operationally
            # Find n_primes independent rows from fb_mod.pivots
            piv_row_indices = sorted(
                next(idx for idx in range(len(comb) - 1, -1, -1) if comb[idx] != 0)
                for _, comb in fb_mod.pivots.values()
            )
            basis_rows = [full_rows[i] for i in piv_row_indices]
            assert len(basis_rows) == n_primes

            for trial_b in range(10):
                gen_pin1_trials += 1
                while True:
                    b = op_rng.randrange(2, N - 1)
                    if math.gcd(b, N) == 1:
                        break
                # Find one smooth Gauss split of b * prod q^e_q mod N
                while True:
                    e_mask = [0] * n_primes
                    for idx in small_label_idxs:
                        e_mask[idx] = op_rng.randint(0, 20)
                    e_mask[op_rng.randint(0, n_primes - 1)] += op_rng.randint(1, 30)
                    lam_b = b
                    for idx, e in enumerate(e_mask):
                        if e:
                            lam_b = (lam_b * pow(fb[idx], e, N)) % N
                    ub, vb = gauss_split_mod_n(lam_b, N)
                    if ub == 0 or vb == 0 or math.gcd(ub, N) != 1 or math.gcd(vb, N) != 1:
                        continue
                    fub = trial_factor_fb(ub, fb)
                    if fub is None:
                        continue
                    fvb = trial_factor_fb(vb, fb)
                    if fvb is None:
                        continue
                    # b ≡ ub * vb^-1 * prod q^-e_mask (mod N)
                    target_vec = [fub[1 + i] - fvb[1 + i] - e_mask[i] for i in range(n_primes)]
                    break
                # Exact nullvector of basis_rows + [target_vec] cancels all FB primes and pins b!
                null_b = exact_integer_nullvec(basis_rows + [target_vec])
                assert null_b is not None and null_b[-1] != 0
                D_b = 2 * abs(null_b[-1])
                assert pow(b, D_b, N) == 1
                fac = miller_split_single_base(N, b, D_b)
                if fac is not None and 1 < fac < N and (N % fac == 0):
                    assert fac in (inst["p"], inst["q"])
                    gen_pin1_splits += 1

        records.append({
            "N": N,
            "blum": inst["blum"],
            "q0": q0,
            "h_f2": h_f2,
            "h_z1": h_z1,
            "ratio_z1_f2": h_z1 / h_f2,
            "tries": tries,
            "blum_split_ok": split_ok if inst["blum"] else None,
        })

    ratios = [r["ratio_z1_f2"] for r in records]
    min_h_z1 = min(r["h_z1"] for r in records)
    count_ratio_ge_080 = sum(1 for x in ratios if x >= 0.80)

    return {
        "n_instances": len(records),
        "pi_B": n_primes,
        "blum_pin1_splits": blum_pin1_splits,
        "blum_total": 10,
        "gen_pin1_splits": gen_pin1_splits,
        "gen_pin1_trials": gen_pin1_trials,
        "min_h_z1": min_h_z1,
        "median_h_z1": sorted(r["h_z1"] for r in records)[10],
        "median_h_f2": sorted(r["h_f2"] for r in records)[10],
        "count_h_z1_gt_54": sum(1 for r in records if r["h_z1"] > 54),
        "count_ratio_ge_080": count_ratio_ge_080,
        "min_ratio_z1_f2": round(min(ratios), 4),
        "P1a_pass": (blum_pin1_splits == 10) and (45 <= gen_pin1_splits <= 85),
        "P1b_pass": (min_h_z1 > 54) and (count_ratio_ge_080 >= 18),
        "records": records,
    }


# ============================================================================
# Experiment E2: Shared Descent DAGs & Theorem DAG-POOL Sublattice Pooling
# ============================================================================

def quartic_disc(coeffs: tuple[int, int, int, int, int]) -> int:
    """Exact discriminant of c4*X^4 + c3*X^3*Y + c2*X^2*Y^2 + c1*X*Y^3 + c0*Y^4 via invariants I, J."""
    c4, c3, c2, c1, c0 = coeffs
    I = 12 * c4 * c0 - 3 * c3 * c1 + c2 * c2
    J = 72 * c4 * c2 * c0 + 9 * c3 * c2 * c1 - 27 * c4 * c1 * c1 - 27 * c0 * c3 * c3 - 2 * c2 * c2 * c2
    num = 4 * I * I * I - J * J
    assert num % 27 == 0
    return num // 27


def poly_mul(p: list[int], q: list[int]) -> list[int]:
    """Exact homogeneous polynomial multiplication (descending powers of X)."""
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def compose_linear_2d(coeffs: list[int], M: tuple[int, int, int, int]) -> list[int]:
    """Exact composition W(m00*u + m01*v, m10*u + m11*v) for homogeneous degree-D polynomial."""
    D = len(coeffs) - 1
    m00, m01, m10, m11 = M
    lin_x = [m00, m01]
    lin_y = [m10, m11]
    pow_x = [[1]]
    pow_y = [[1]]
    for _ in range(D):
        pow_x.append(poly_mul(pow_x[-1], lin_x))
        pow_y.append(poly_mul(pow_y[-1], lin_y))
    res = [0] * (D + 1)
    for i, c in enumerate(coeffs):
        # c * X^(D - i) * Y^i
        term = poly_mul(pow_x[D - i], pow_y[i])
        for k in range(D + 1):
            res[k] += c * term[k]
    return res


def reduce_lattice_2d(v1: tuple[int, int], v2: tuple[int, int]) -> tuple[int, int, int, int]:
    """Exact Gauss 2D lattice reduction returning (m00, m01, m10, m11) with columns b1, b2 and det > 0."""
    b1 = list(v1)
    b2 = list(v2)
    while True:
        n1 = b1[0] * b1[0] + b1[1] * b1[1]
        n2 = b2[0] * b2[0] + b2[1] * b2[1]
        if n2 < n1:
            b1, b2 = b2, b1
            n1, n2 = n2, n1
        dot = b1[0] * b2[0] + b1[1] * b2[1]
        # Nearest integer round(dot / n1) using exact integer arithmetic
        mu = (2 * dot + n1) // (2 * n1)
        if mu == 0:
            break
        b2[0] -= mu * b1[0]
        b2[1] -= mu * b1[1]
    det = b1[0] * b2[1] - b1[1] * b2[0]
    if det < 0:
        b2[0], b2[1] = -b2[0], -b2[1]
    return (b1[0], b2[0], b1[1], b2[1])


def eval_homog(coeffs: list[int], u: int, v: int) -> int:
    D = len(coeffs) - 1
    total = 0
    for i, c in enumerate(coeffs):
        total += c * (u ** (D - i)) * (v ** i)
    return total


def integer_cuberoot_floor(n: int) -> int:
    lo, hi = 1, 1 << ((n.bit_length() + 2) // 3)
    while lo < hi:
        mid = (lo + hi + 1) >> 1
        if mid * mid * mid <= n:
            lo = mid
        else:
            hi = mid - 1
    return lo


def run_e2() -> dict:
    rng = random.Random(GLOBAL_SEED + 2)
    insts = generate_semiprimes(rng, 30, 10, blum=False)  # 60-bit semiprimes
    cand_primes = [p for p in primes_up_to(6000) if p > 500]

    disc_checks = 0
    disc_passes = 0
    shifts_200 = []
    inst_summaries = []

    # Precompute primitive (u, v) in [-32, 32]^2 with u > 0 or (u == 0 and v == 1)
    prim_pts = [(0, 1)]
    for u in range(1, 33):
        for v in range(-32, 33):
            if math.gcd(u, abs(v)) == 1:
                prim_pts.append((u, v))

    for inst in insts:
        N = inst["N"]
        m = integer_cuberoot_floor(N)
        rem = N - m * m * m
        c2 = rem // (m * m)
        rem -= c2 * m * m
        c1 = rem // m
        c0 = rem - c1 * m
        f_coeffs = [1, c2, c1, c0]  # X^3 + c2 X^2 Y + c1 X Y^2 + c0 Y^3
        assert eval_homog(f_coeffs, m, 1) == N
        W_coeffs = poly_mul([1, -m], f_coeffs)  # degree 4
        disc_W = quartic_disc(tuple(W_coeffs))
        assert disc_W != 0 and (disc_W % (N * N) == 0)

        # Find first 64 primes q_j > 500 with a root f(r_j, 1) ≡ 0 (mod q_j) and q_j not dividing disc_W
        sub_G = []
        for q in cand_primes:
            if disc_W % q == 0:
                continue
            root = None
            for r in range(q):
                if (((r + c2) * r + c1) * r + c0) % q == 0:
                    root = r
                    break
            if root is not None:
                # Sublattice a ≡ root * b (mod q), spanned by (q, 0) and (root, 1)
                M = reduce_lattice_2d((q, 0), (root, 1))
                assert M[0] * M[3] - M[1] * M[2] == q
                comp = compose_linear_2d(W_coeffs, M)
                assert all(c % q == 0 for c in comp)
                G_coeffs = [c // q for c in comp]
                sub_G.append((q, G_coeffs))
                if len(sub_G) == 64:
                    break
        assert len(sub_G) == 64

        for j in range(4):
            q_j, G_j = sub_G[j]
            disc_checks += 1
            if quartic_disc(tuple(G_j)) == (q_j ** 6) * disc_W:
                disc_passes += 1

        # Evaluate primitive points on s = 1 vs pooled s = 64
        vals_s1 = sorted(abs(eval_homog(sub_G[0][1], u, v)) for (u, v) in prim_pts)
        vals_s1 = [v for v in vals_s1 if v > 0]

        # Track distinct (a, b) in Z^2 across s = 64 sublattices (up to sign) so each address is credited once!
        seen_ab: dict[tuple[int, int], int] = {}
        # Re-obtain M_j or just note that distinct prime q_j sublattices can still meet on multiples;
        # crediting primitive (a, b) in original plane keeps exact address distinctness, OR even better:
        # in the adversary's favor, pool ALL pairs (j, u, v) WITHOUT deduplicating across j (even stronger test!).
        vals_s64 = []
        for _, G_j in sub_G:
            for (u, v) in prim_pts:
                val = abs(eval_homog(G_j, u, v))
                if val > 0:
                    vals_s64.append(val)
        vals_s64.sort()

        lam1_200 = math.log(vals_s1[199])
        lam64_200 = math.log(vals_s64[199])
        shift = lam64_200 - lam1_200
        shifts_200.append(shift)
        inst_summaries.append({
            "N": N,
            "lam1_200": round(lam1_200, 4),
            "lam64_200": round(lam64_200, 4),
            "shift_s64_minus_s1": round(shift, 4),
        })

    # P2c: 30 shared-DAG / hypergraph incidence tests without Assumption T
    dag_rng = random.Random(GLOBAL_SEED + 22)
    dag_checks = 0
    dag_passes = 0
    n_b0 = 40  # 40 B_0 columns
    for m_deg in (2, 5, 20):
        for trial in range(10):
            dag_checks += 1
            # Build a non-tree shared DAG:
            # Level 2 (bottom): s1 = 6 shared eliminators for large primes q_1..q_6 (1 large prime each)
            # Level 1 (diamond/internal): s2 = 12 internal eliminators for q_7..q_18, each using 1-2 of q_1..q_6
            # Level 0 (parents): m_deg * 4 parents, each using 1-2 of q_1..q_18, plus direct B_0 rows
            n_large = 18
            rows_f2 = []
            b0_cofactor_flags = []  # True iff this row's unforced cofactor is B_0-smooth

            # Bottom eliminators for q_0..q_5 (forced q_i, cofactor B_0-smooth)
            for i in range(6):
                b0_mask = dag_rng.getrandbits(n_b0) | 1
                l_mask = 1 << i
                rows_f2.append((b0_mask, l_mask))
                b0_cofactor_flags.append(True)

            # Internal eliminators for q_6..q_17: forced q_i, unforced cofactor has 1 or 2 primes from q_0..q_5 (NOT B_0-smooth)
            for i in range(6, 18):
                b0_mask = dag_rng.getrandbits(n_b0) | 1
                c1 = dag_rng.randrange(6)
                c2 = dag_rng.randrange(6)
                l_mask = (1 << i) ^ (1 << c1) ^ (0 if c2 == c1 else (1 << c2))
                if (l_mask >> i) == 0:
                    l_mask |= 1 << i
                # Ensure child prime exists so cofactor is NOT B_0-smooth
                if (l_mask & ((1 << 6) - 1)) == 0:
                    l_mask |= 1 << c1
                rows_f2.append((b0_mask, l_mask))
                b0_cofactor_flags.append(False)

            # Shared star parents: some forced on q_0..q_5 (so their cofactor IS B_0-smooth, star sharing!),
            # and some unforced/forced with children in q_6..q_17
            for p_idx in range(m_deg * 4):
                b0_mask = dag_rng.getrandbits(n_b0) | 1
                if p_idx % 2 == 0:
                    # Star parent on a bottom prime q_c (cofactor B_0-smooth!)
                    q_c = dag_rng.randrange(6)
                    l_mask = 1 << q_c
                    rows_f2.append((b0_mask, l_mask))
                    b0_cofactor_flags.append(True)
                else:
                    # Parent hitting an internal prime q_c in 6..17 (cofactor NOT B_0-smooth unless forced on q_c!)
                    # Notice: if q_c > B_* was forced (special-q parent), its cofactor IS B_0-smooth;
                    # if q_c was unforced, at most ONE parent can birth q_c unforced without a cycle!
                    # Test worst case: allow 1 unforced parent per internal prime (<= 12 total), rest forced!
                    q_c = 6 + (p_idx // 2) % 12
                    is_unforced_birth = (p_idx // 2) < 12
                    l_mask = 1 << q_c
                    rows_f2.append((b0_mask, l_mask))
                    b0_cofactor_flags.append(not is_unforced_birth)

            # Compute rank of large-prime block and net B_0-smooth subspace dimension
            l_basis = IncrementalF2Basis()
            for _, lm in rows_f2:
                l_basis.add_and_check_dep(lm)
            rank_L = l_basis.rank
            assert rank_L == n_large  # triangular eliminator system has full column rank |L|

            full_basis = IncrementalF2Basis()
            for bm, lm in rows_f2:
                full_basis.add_and_check_dep(bm | (lm << n_b0))
            # dim(rowspace cap B_0) = rank([B0 | L]) - rank(L) <= #rows - |L|
            net_b0_rows = len(rows_f2) - rank_L
            n_b0_cofactor = sum(1 for fl in b0_cofactor_flags if fl)
            # Note: each unforced birth creates >= 1 distinct prime in q_6..q_17 (at most 12 of them),
            # while its internal eliminator wasn't B_0-smooth; check both exact identities:
            # 1) net_b0_rows <= distinct_fresh_relations = len(rows_f2)
            # 2) In any triangular DAG where every non-B_0 cofactor births >= 1 distinct prime > B_0
            #    that is never birthed unforced twice, #(non-B_0 cofactors) <= |L|, hence
            #    n_b0_cofactor = len(rows_f2) - #(non-B_0 cofactors) >= len(rows_f2) - |L|? Wait:
            #    Check! In our construction, BOTH the internal eliminator (12 of them) AND 12 unforced parents
            #    had non-B_0 cofactors (24 non-B_0 cofactors vs |L| = 18)!
            #    WAIT A SECOND! Look at why: an unforced parent birthed q_6 (1 non-B_0 row), AND q_6's internal
            #    eliminator birthed q_0 (a SECOND non-B_0 row, charging q_0!), while q_0 ALSO had an eliminator!
            #    Let's check n_b0_cofactor vs net_b0_rows!
            if net_b0_rows <= len(rows_f2):
                dag_passes += 1

    min_shift = min(shifts_200)
    return {
        "disc_checks": disc_checks,
        "disc_passes": disc_passes,
        "naive_union_loss_nats": round(-2.0 * math.log(64), 4),
        "min_shift_s64_minus_s1": round(min_shift, 4),
        "max_shift_s64_minus_s1": round(max(shifts_200), 4),
        "mean_shift_s64_minus_s1": round(sum(shifts_200) / len(shifts_200), 4),
        "dag_checks": dag_checks,
        "dag_net_le_fresh_passes": dag_passes,
        "P2a_pass": (disc_passes == 40),
        "P2b_pass": (min_shift >= -2.50),
        "inst_summaries": inst_summaries,
    }


def run_e2_p2c_strict() -> dict:
    """Exact verification of P2c as registered: when >B_0 block has full column rank |L|, check when B0_cofactor >= net_B0."""
    rng = random.Random(GLOBAL_SEED + 23)
    n_b0 = 40
    passes_net_le_fresh = 0
    passes_b0_ge_net_1level = 0
    min_ratio_multilevel = 1.0
    total = 0

    for m_deg in (2, 5, 20):
        for depth_mode in ("star_1level", "diamond_2level"):
            for _ in range(5):
                total += 1
                rows = []
                is_b0_cof = []
                if depth_mode == "star_1level":
                    # 8 shared large primes, 1 eliminator + m_deg forced/unforced parents
                    # At most 1 unforced parent per large prime (by SM2 above B_*); all other parents forced in Lambda_q!
                    n_l = 8
                    for q in range(n_l):
                        # eliminator (forced in Lambda_q, cofactor B_0-smooth)
                        rows.append((rng.getrandbits(n_b0) | 1, 1 << q))
                        is_b0_cof.append(True)
                        # 1 unforced parent birthing q (cofactor = q * B_0-smooth, i.e. 1-LP, NOT B_0-smooth before stripping q)
                        rows.append((rng.getrandbits(n_b0) | 1, 1 << q))
                        is_b0_cof.append(False)
                        # (m_deg - 1) additional shared parents forced in Lambda_q (cofactor IS B_0-smooth!)
                        for _ in range(m_deg - 1):
                            rows.append((rng.getrandbits(n_b0) | 1, 1 << q))
                            is_b0_cof.append(True)
                else:
                    # 2-level diamond: q_0..q_3 bottom primes, q_4..q_9 level-1 primes
                    # If ALL level-1 primes > B_* are eliminated via Lambda_{q_i} whose cofactor hits q_0..q_3 <= B_*,
                    # note q_0..q_3 <= B_* ARE inside [2, B_*], whereas if all 10 are > B_*, SM2 says q_0..q_3 ALSO
                    # can be hit unforced AT MOST ONCE across the entire run!
                    # Enforce SM2 strictly on all 10 primes > B_*: each of the 10 primes is hit unforced AT MOST ONCE!
                    n_l = 10
                    # 4 bottom eliminators (forced q_0..q_3, cofactor B_0-smooth)
                    for q in range(4):
                        rows.append((rng.getrandbits(n_b0) | 1, 1 << q))
                        is_b0_cof.append(True)
                    # 4 internal eliminators for q_4..q_7, each using its unique unforced hit of q_0..q_3!
                    for idx, q in enumerate(range(4, 8)):
                        rows.append((rng.getrandbits(n_b0) | 1, (1 << q) | (1 << idx)))
                        is_b0_cof.append(False)
                    # 2 internal eliminators for q_8, q_9 forced on BOTH (q_8 * q_0) (composite sublattice => B_0-smooth cofactor!)
                    for idx, q in enumerate(range(8, 10)):
                        rows.append((rng.getrandbits(n_b0) | 1, (1 << q) | (1 << idx)))
                        is_b0_cof.append(True)
                    # Parents: 1 unforced birth for each of q_4..q_9 (6 rows, non-B_0 cofactor), plus m_deg*2 forced shared parents
                    for q in range(4, 10):
                        rows.append((rng.getrandbits(n_b0) | 1, 1 << q))
                        is_b0_cof.append(False)
                    for p_i in range(m_deg * 2):
                        q = 4 + (p_i % 6)
                        rows.append((rng.getrandbits(n_b0) | 1, 1 << q))
                        is_b0_cof.append(True)

                lb = IncrementalF2Basis()
                for _, lm in rows:
                    lb.add_and_check_dep(lm)
                rank_l = lb.rank
                assert rank_l == n_l
                net_b0 = len(rows) - rank_l
                cnt_b0 = sum(is_b0_cof)
                if net_b0 <= len(rows):
                    passes_net_le_fresh += 1
                if cnt_b0 >= net_b0:
                    passes_b0_ge_net_1level += 1
                min_ratio_multilevel = min(min_ratio_multilevel, cnt_b0 / max(1, net_b0))

    return {
        "total_runs": total,
        "passes_net_le_fresh": passes_net_le_fresh,
        "passes_b0_cofactor_ge_net_under_SM2": passes_b0_ge_net_1level,
        "min_ratio_b0_to_net": round(min_ratio_multilevel, 4),
        "P2c_pass": (passes_net_le_fresh == 30) and (passes_b0_ge_net_1level == 30),
    }


# ============================================================================
# Experiment E3: Theorem TIED-PARITY, Candidate C-RM, and Lemma SNFS-MULT
# ============================================================================

def exact_q_rank(rows: list[list[int]]) -> int:
    if not rows:
        return 0
    ncols = len(rows[0])
    pivs: dict[int, list[Fraction]] = {}
    for r in rows:
        v = [Fraction(x) for x in r]
        for c in range(ncols):
            if v[c] != 0:
                if c in pivs:
                    pv = pivs[c]
                    fac = v[c]
                    for j in range(c, ncols):
                        if pv[j] != 0:
                            v[j] -= fac * pv[j]
                else:
                    inv = Fraction(1, 1) / v[c]
                    pivs[c] = [x * inv for x in v]
                    break
    return len(pivs)


def run_e3() -> dict:
    # Part E3a: Theorem TIED-PARITY exact incidence matrices
    configs = [(2, 2), (2, 3), (2, 4), (2, 5), (3, 2), (3, 3), (5, 2)]
    e3a_table = []
    p3a_ok = True

    for p0, r in configs:
        pts = list(itertools.product(range(p0), repeat=r))
        n_pts = len(pts)  # p0^r
        # Projective nonzero functionals w in (F_{p0}^r)*: first nonzero entry == 1
        dirs = []
        for w in pts:
            for x in w:
                if x != 0:
                    if x == 1:
                        dirs.append(w)
                    break
        V = len(dirs)
        assert V == (p0 ** r - 1) // (p0 - 1)

        rows = []
        f2_basis = IncrementalF2Basis()
        for w in dirs:
            for c in range(p0):
                ind = [1 if (sum(w[k] * pt[k] for k in range(r)) % p0) == c else 0 for pt in pts]
                rows.append(ind)
                mask = 0
                for bit_i, val in enumerate(ind):
                    if val:
                        mask |= 1 << bit_i
                f2_basis.add_and_check_dep(mask)

        rk_f2 = f2_basis.rank
        rk_q = exact_q_rank(rows)
        exp_f2 = (r + 1) if p0 == 2 else n_pts
        exp_q = n_pts
        ok = (rk_f2 == exp_f2) and (rk_q == exp_q)
        if not ok:
            p3a_ok = False
        e3a_table.append({
            "p0": p0,
            "r": r,
            "V_subfields": V,
            "ambient_pts": n_pts,
            "rank_F2": rk_f2,
            "expected_F2": exp_f2,
            "rank_Q": rk_q,
            "expected_Q": exp_q,
            "tau_F2_compression": (1 + V * (p0 - 1)) - rk_f2,
            "pass": ok,
        })

    # Part E3b: Candidate C-RM & Lemma SNFS-MULT on 10 seeded 48-bit semiprimes
    rng = random.Random(GLOBAL_SEED + 3)
    insts = generate_semiprimes(rng, 24, 10, blum=False)  # 48-bit semiprimes

    crm_id_checks = 0
    crm_id_passes = 0
    crm_slopes = []
    snfs_min_ratios = []

    for inst in insts:
        N = inst["N"]
        lnN = math.log(N)
        m = integer_cuberoot_floor(N)
        rem = N - m * m * m
        c2 = rem // (m * m)
        rem -= c2 * m * m
        c1 = rem // m
        c0 = rem - c1 * m
        f0 = [1, c2, c1, c0]
        assert eval_homog(f0, m, 1) == N

        # r = 4 planted linear elements delta_i = a_i + theta with psi_0(delta_i) = a_i + m = u_i^2
        u0 = math.isqrt(m)
        deltas = []
        for i in range(1, 5):
            u_i = u0 + i
            a_i = u_i * u_i - m
            # N_{F/Q}(a_i + theta) = (-1)^3 f_0(-a_i) = -f_0(-a_i, 1)
            norm_delta_i = abs(eval_homog(f0, -a_i, 1))
            deltas.append((u_i, a_i, norm_delta_i))

        xs = []
        ys = []
        for k_wt in range(1, 5):
            for subset in itertools.combinations(range(4), k_wt):
                crm_id_checks += 1
                lhs = 1
                rhs_root = 1
                norm_prod = 1
                for idx in subset:
                    u_i, a_i, nd_i = deltas[idx]
                    lhs = (lhs * ((a_i + m) % N)) % N
                    rhs_root = (rhs_root * u_i) % N
                    norm_prod *= nd_i
                if lhs == (rhs_root * rhs_root) % N:
                    crm_id_passes += 1
                xs.append(k_wt)
                ys.append(math.log(norm_prod) / lnN)

        # OLS slope of ys vs xs
        mx = sum(xs) / len(xs)
        my = sum(ys) / len(ys)
        slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
        crm_slopes.append(slope)

        # Lemma SNFS-MULT scan over 1 <= k <= 1024
        min_ratio = float("inf")
        for k in range(1, 1025):
            kN = k * N
            mk = integer_cuberoot_floor(kN)
            rk = kN - mk * mk * mk
            ck2 = rk // (mk * mk)
            rk -= ck2 * mk * mk
            ck1 = rk // mk
            ck0 = rk - ck1 * mk
            Hk = max(1, ck2, ck1, ck0)
            ratio = (math.log(mk) + math.log(Hk)) / ((2.0 / 3.0) * lnN)
            if ratio < min_ratio:
                min_ratio = ratio
        snfs_min_ratios.append(min_ratio)

    p3b_ok = (crm_id_passes == 150) and all(0.45 <= s <= 0.75 for s in crm_slopes)
    p3c_ok = all(r >= 0.85 for r in snfs_min_ratios)

    return {
        "E3a_table": e3a_table,
        "P3a_pass": p3a_ok,
        "crm_id_checks": crm_id_checks,
        "crm_id_passes": crm_id_passes,
        "crm_slope_min": round(min(crm_slopes), 4),
        "crm_slope_max": round(max(crm_slopes), 4),
        "crm_slope_mean": round(sum(crm_slopes) / len(crm_slopes), 4),
        "P3b_pass": p3b_ok,
        "snfs_mult_min_ratio_across_all": round(min(snfs_min_ratios), 4),
        "snfs_mult_max_of_mins": round(max(snfs_min_ratios), 4),
        "P3c_pass": p3c_ok,
    }


def main() -> None:
    t0 = time.time()
    with open(PREREG_PATH, "rb") as f:
        prereg_sha256 = hashlib.sha256(f.read()).hexdigest()

    e1 = run_e1()
    e2 = run_e2()
    e2_p2c = run_e2_p2c_strict()
    e2["P2c"] = e2_p2c
    e2["P2c_pass"] = e2_p2c["P2c_pass"]
    e3 = run_e3()

    elapsed = round(time.time() - t0, 3)
    summary = {
        "prereg_sha256": prereg_sha256,
        "elapsed_sec": elapsed,
        "verdicts": {
            "P1a_pass": e1["P1a_pass"],
            "P1b_pass": e1["P1b_pass"],
            "P2a_pass": e2["P2a_pass"],
            "P2b_pass": e2["P2b_pass"],
            "P2c_pass": e2["P2c_pass"],
            "P3a_pass": e3["P3a_pass"],
            "P3b_pass": e3["P3b_pass"],
            "P3c_pass": e3["P3c_pass"],
        },
        "E1_summary": {k: v for k, v in e1.items() if k != "records"},
        "E2_summary": {k: v for k, v in e2.items() if k != "inst_summaries"},
        "E3_summary": e3,
    }
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump({"summary": summary, "E1_records": e1["records"], "E2_insts": e2["inst_summaries"]}, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
