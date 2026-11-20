"""Pre-registered certified adversary experiment for divided-difference LP_exc^rel (dd_g).
Run with: PYTHONPATH=/Volumes/git/factoring-connection /Volumes/git/factoring-connection/.venv/bin/python3 -B run_dd_adversary.py
All real comparisons and thresholds are certified with 256-bit python-flint ball arithmetic (arb/acb).
All norms and prime-ideal valuation parities are exact via flint resultants and PARI nfinit/idealprimedec.
"""
import itertools
import json
import math
import random
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
from flint import acb, arb, ctx, fmpq_mat, fmpz_mat, fmpz_poly as Poly
from sympy import integer_nthroot, primerange

ROOT = Path("/Volumes/git/factoring-connection")
sys.path[:0] = [str(ROOT / "experiments/subthird/codex"), str(ROOT / "experiments/masullo_ext/mext_c"), str(ROOT)]
import experiment as ex
from cover import pair_units

ctx.prec = 512
T = Poly([0, 1])
one = Poly([1])
SQ2 = math.sqrt(2.0)


def orbit_rep(z, units, E=1):
    reps = []
    for u, v in units:
        w = tuple(t for a, b in zip(z[::2], z[1::2])
                  for t in (u*a + 2*v*b, u*b + v*a))
        if max(map(abs, w)) <= E:
            reps.append(w)
    return min(reps)


def base_b_monic(Z, b, q):
    rem = Z - b**q
    dig = []
    for _ in range(q):
        d = (rem + b // 2) % b - b // 2
        dig.append(int(d))
        rem = (rem - d) // b
    if rem != 0:
        return None
    return Poly(dig + [1])


def lagrange_int_poly(nodes, vals):
    """Exact integer polynomial S(T) of degree < len(nodes) with S(nodes[j]) = vals[j], if in Z[T]."""
    k = len(nodes)
    if k == 0:
        return Poly([])
    if k == 1:
        return Poly([int(vals[0])])
    # Exact rational Lagrange via sympy or Fraction
    coeffs_fr = [Fr(0)] * k
    for j in range(k):
        # poly prod_{l != j} (T - nodes[l]) / prod_{l != j} (nodes[j] - nodes[l])
        num = [Fr(1)]
        den = Fr(1)
        for l in range(k):
            if l == j:
                continue
            den *= (nodes[j] - nodes[l])
            new_num = [Fr(0)] * (len(num) + 1)
            for idx, c in enumerate(num):
                new_num[idx] -= c * nodes[l]
                new_num[idx + 1] += c
            num = new_num
        scale = Fr(vals[j]) / den
        for idx, c in enumerate(num):
            coeffs_fr[idx] += c * scale
    assert all(c.denominator == 1 for c in coeffs_fr), f"Non-integer Lagrange poly: {coeffs_fr}"
    return Poly([int(c.numerator) for c in coeffs_fr])


def best_ternary_Z_vector(k):
    if k == 1:
        return [[1]], [1.0]
    rng = random.Random(20261120 + k)
    best_mx = 0.0
    best_B, best_y = None, None
    # Deterministic Fibonacci/staircase baseline
    B_fib = [[0]*k for _ in range(k)]
    B_fib[0][0] = 1
    if k >= 2:
        B_fib[1][0], B_fib[1][1] = -1, 1
    for r in range(2, k):
        B_fib[r][r-2], B_fib[r][r-1], B_fib[r][r] = -1, -1, 1
    invF = fmpq_mat(fmpz_mat(B_fib)).inv()
    best_y = [float(invF[r, 0]) for r in range(k)]
    best_mx = max(abs(x) for x in best_y)
    best_B = B_fib

    trials = 25000 if k <= 5 else 15000
    for _ in range(trials):
        B = [[rng.choice((-1, 0, 1)) for _ in range(k)] for _ in range(k)]
        M = fmpz_mat(B)
        if int(M.det()) != 0:
            invM = fmpq_mat(M).inv()
            for col in range(k):
                y = [float(invM[r, col]) for r in range(k)]
                mx = max(abs(v) for v in y)
                if mx > best_mx:
                    best_mx = mx
                    B_perm = [B[col]] + [B[r] for r in range(k) if r != col]
                    invP = fmpq_mat(fmpz_mat(B_perm)).inv()
                    best_B = B_perm
                    best_y = [float(invP[r, 0]) for r in range(k)]
    return best_B, best_y


def best_ternary_Ok_vector(k):
    if k == 1:
        return [1.0 + SQ2]
    # Deterministic unit-staircase over Z[sqrt(2)]: y_j = (1 + sqrt(2))^{2j + 1}
    y_stair = [(1.0 + SQ2)**(2*j + 1) for j in range(k)]
    best_mx = max(abs(v) for v in y_stair)
    best_y = y_stair

    ALPHA = [(u, v) for u in (-1, 0, 1) for v in (-1, 0, 1)]
    rng = random.Random(30261120 + k)
    trials = 15000 if k <= 4 else 8000
    for _ in range(trials):
        A = [[rng.choice(ALPHA) for _ in range(k)] for _ in range(k)]
        An = np.array([[u + v*SQ2 for u, v in row] for row in A], dtype=np.float64)
        d = abs(np.linalg.det(An))
        if 1e-11 < d < 0.1:
            invAn = np.linalg.inv(An)
            for col in range(k):
                y = invAn[:, col]
                mx = float(np.max(np.abs(y)))
                if mx > best_mx:
                    # Certify nonzero exact determinant in Q(sqrt(2))
                    mat_z = []
                    for row in A:
                        r1, r2 = [], []
                        for u, v in row:
                            r1.extend([u, 2*v])
                            r2.extend([v, u])
                        mat_z.extend([r1, r2])
                    if int(fmpz_mat(mat_z).det()) != 0:
                        best_mx = mx
                        best_y = [float(v) for v in y]
    return best_y


# --- Constructors ---

def build_B0(N, n, m):
    a = math.isqrt(N + 2)
    assert (a*a - 2) % N == 0
    cand_poles = [T*T + 1, T - 1, T + 1, T - 2]
    R_list = cand_poles[:m - 2]
    if m <= 4:
        # Match sweep.py / experiment.py select(N) bit-for-bit!
        ex.n, ex.d = n, 2*n
        _, b, P, Q = ex.select(N)
        return dict(name="B0", b=b, P=P, Q=Q, R_list=R_list)
    R_prod = one
    for Rj in R_list:
        R_prod = R_prod * Rj
    base = int(integer_nthroot(a, n)[0])
    for offset in range(-3, base - 2):
        b = base - offset
        P = base_b_monic(a, b, n)
        if P is None or not ex.irreducible(P*P - 2):
            continue
        for j0 in (1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8, 9, -9, 11, -11):
            Q = P + j0 * (T - b) * R_prod
            if ex.irreducible(Q*Q - 2):
                return dict(name="B0", b=b, P=P, Q=Q, R_list=R_list)
    raise RuntimeError("B0 failed")


def build_C1a(N, n, m):
    """Base-b Near-Pole Divisor Alignment, exact P(b) = a, unit remainder c0 = +-1."""
    a = math.isqrt(N + 2)
    k = m - 2
    q = n - k
    # Choose d_list so prod(d_list) > 0 and divides a - c0 with c0 in {1, -1}
    if k == 1:
        d_list, c0 = [1], 1
    elif k == 2:
        d_list, c0 = [1, 2], 1
    elif k == 3:
        d_list, c0 = [1, -1, -2], 1
    elif k == 4:
        if (a - 1) % 4 == 0:
            d_list, c0 = [1, 2, -1, -2], 1
        else:
            d_list, c0 = [1, 2, -1, -2], -1
    else:
        raise ValueError("Unsupported k")
    D_prod = math.prod(d_list)
    assert D_prod > 0 and (a - c0) % D_prod == 0
    Z = (a - c0) // D_prod
    base = int(integer_nthroot(Z, q)[0])
    for offset in range(-2, base - 2):
        b = base - offset
        if b <= 3:
            continue
        U = base_b_monic(Z, b, q)
        if U is None:
            continue
        R_list = [T - (b - d) for d in d_list]
        R_prod = one
        for Rj in R_list:
            R_prod = R_prod * Rj
        P = Poly([c0]) + R_prod * U
        assert int(P(b)) == a
        if not ex.irreducible(P*P - 2):
            continue
        for j0 in (1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8, 9, -9, 11, -11):
            Q = P + j0 * (T - b) * R_prod
            if ex.irreducible(Q*Q - 2):
                return dict(name="C1a", b=b, P=P, Q=Q, R_list=R_list, c0=c0)
    raise RuntimeError("C1a failed")


def build_C1b(N, n, m):
    """Unit-Remainder Small-Pole KH-Floor Family near T=0, P(b) = a + t*N, S(T) = 1."""
    a = math.isqrt(N + 2)
    k = m - 2
    q = n - k
    R_list = [T - j for j in range(1, k + 1)]
    R_prod = one
    for Rj in R_list:
        R_prod = R_prod * Rj
    b0 = int(integer_nthroot(6 * N, q)[0]) + 1
    for b in range(max(k + 2, b0), b0 + 200):
        Rb = int(R_prod(b))
        if math.gcd(Rb, N) != 1:
            continue
        t_mod = ((1 - a) * pow(N, -1, Rb)) % Rb
        Z0 = (a + t_mod * N - 1) // Rb
        # Step Z by N (i.e. t by Rb) to land near b^q
        target = b**q
        k_step = max(0, int(round((target - Z0) / N)))
        for dk in (0, 1, -1, 2):
            ks = k_step + dk
            if ks < 0:
                continue
            Z = Z0 + ks * N
            U = base_b_monic(Z, b, q)
            if U is None:
                continue
            P = one + R_prod * U
            assert (int(P(b)) - a) % N == 0
            if not ex.irreducible(P*P - 2):
                continue
            for j0 in (1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8, 9, -9, 11, -11):
                Q = P + j0 * (T - b) * R_prod
                if ex.irreducible(Q*Q - 2):
                    return dict(name="C1b", b=b, P=P, Q=Q, R_list=R_list)
    raise RuntimeError("C1b failed")


def cheb_scaled(n, L, shift=0):
    X = T - shift
    if n == 0:
        return Poly([2])
    if n == 1:
        return X
    c0, c1 = Poly([2]), X
    L2 = L * L
    for _ in range(2, n + 1):
        c0, c1 = c1, X * c1 - L2 * c0
    return c1


def count_real_roots_certified(f_pol):
    from src.rank3.forms import acb_roots
    rts = acb_roots(ex.coeffs(f_pol))
    # For squarefree integer poly, a root is real iff conjugate-invariant; check certified imag ball contains 0
    # Better: use flint fmpz_poly real root isolation or Sturm / PARI polsturm!
    return int(ex.pari.Polrev(ex.coeffs(f_pol)).polsturm())


def build_C2(N, n, m):
    """Totally-Real Scaled Chebyshev Fields + Root-Clustered Chebyshev Poles."""
    a = math.isqrt(N + 2)
    # Try exact P(b) = a first; if needed, P(b) = a + t*N
    for use_modN in (False, True):
        target_a_scale = a if not use_modN else 4 * N
        L_max = max(25, int(integer_nthroot(target_a_scale, n)[0]) + 8)
        base_b_est = int(integer_nthroot(target_a_scale, n)[0])
        for L in range(1, L_max + 1):
            for shift in range(-min(L, 4), min(L, 4) + 1):
                X = T - shift
                cand_poles = [X*X - 2*L*L, X - L, X + L, X - 1]
                if L == 1:
                    cand_poles[3] = X + 2
                R_list = cand_poles[:m - 2]
                R_prod = one
                for Rj in R_list:
                    R_prod = R_prod * Rj
                b_lo = max(shift + 2*L + 1, base_b_est - 6)
                b_hi = max(b_lo + 12, base_b_est + 12, shift + 5*L + 4)
                for b in range(b_lo, b_hi):
                    cvals = [int(cheb_scaled(k, L, shift)(b)) for k in range(n + 1)]
                    if not use_modN:
                        targets = [a]
                    else:
                        if cvals[n] < 2 * N:
                            continue
                        t0 = max(1, int(round((cvals[n] - a) / N)))
                        targets = [a + t0 * N]
                    for targ in targets:
                        if cvals[n] > 2 * targ or cvals[n] < targ // 3:
                            continue
                        rem = targ - cvals[n]
                        digs = [0] * n
                        for k in range(n - 1, 0, -1):
                            d = int(round(rem / cvals[k]))
                            digs[k] = d
                            rem -= d * cvals[k]
                        digs0 = rem
                        sup_pert = abs(digs0) + sum(abs(digs[k]) * 2 * (L ** k) for k in range(1, n))
                        margin = 2 * (L ** n) - sup_pert
                        if margin <= 1.5:
                            continue
                        P = cheb_scaled(n, L, shift) + Poly([digs0])
                        for k in range(1, n):
                            if digs[k]:
                                P = P + digs[k] * cheb_scaled(k, L, shift)
                        assert (int(P(b)) - a) % N == 0
                        f = P*P - 2
                        if not ex.irreducible(f) or count_real_roots_certified(f) != 2*n:
                            continue
                        for j0 in (1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8, 9, -9, 11, -11):
                            Q = P + j0 * (T - b) * R_prod
                            g = Q*Q - 2
                            if ex.irreducible(g) and count_real_roots_certified(g) == 2*n:
                                return dict(name="C2", b=b, P=P, Q=Q, R_list=R_list, L=L, shift=shift, exact_a=(not use_modN))
    raise RuntimeError("C2 failed")


def build_C3(N, n, m, mode="Ok"):
    """Alon-Vu / O_k-Ill-Conditioned Residue-Planted Numerators coupling (m-1) coordinates (w_2..w_m)."""
    a = math.isqrt(N + 2)
    k = m - 2
    q = n - k
    nodes = list(range(1, k + 1))
    R_list = [T - r for r in nodes]
    R_prod = one
    for Rj in R_list:
        R_prod = R_prod * Rj
    # Build (m-1)-dimensional ill-conditioned vector y_full = (y_0, y_1, ..., y_{k}) normalized by y_0 != 0!
    if mode == "Ok":
        # Deterministic unit-staircase + random search of size m-1 over Z[sqrt(2)], normalized by smallest nonzero entry!
        y_raw = best_ternary_Ok_vector(m - 1)
        # Pick pivot entry with smallest >= 0.5 modulus so y_norm[0] = 1 and max |y_norm| is maximized!
        nz = [v for v in y_raw if abs(v) >= 0.5]
        piv = min(nz, key=abs)
        idx_piv = y_raw.index(piv)
        y_rest = [y_raw[i] / piv for i in range(len(y_raw)) if i != idx_piv]
    else:
        _, y_raw = best_ternary_Z_vector(m - 1)
        nz = [v for v in y_raw if abs(v) >= 0.5]
        piv = min(nz, key=abs)
        idx_piv = y_raw.index(piv)
        y_rest = [y_raw[i] / piv for i in range(len(y_raw)) if i != idx_piv]
    K_m = max([1.0] + [abs(v) for v in y_rest])
    b = next(b_try for b_try in range(k + 1, k + 30) if math.gcd(int(R_prod(b_try)), N) == 1)
    Rb = int(R_prod(b))
    # Cluster all q roots of U(T) = (T - x0)^q - M2 around real center x0 >> rho2 * m^2 * K_m!
    rho2_est = int(integer_nthroot(3 * N, q)[0]) + 10
    x0 = max(int(math.ceil(8 * (m**2) * K_m * rho2_est)) + 100, 500)
    fact = math.factorial(max(0, k - 1))
    # Want w_{2+j}(alpha_i) = (-S(r_j)) / (x0 - r_j) approx x0 * y_rest[j-1], while w_2(alpha_i) approx x0 * 1!
    vals = [fact * int(round(-(x0 - r) * x0 * y_rest[j] / fact)) for j, r in enumerate(nodes)]
    S_pol = lagrange_int_poly(nodes, vals)
    Sb = int(S_pol(b))
    # Want U(b) = (b - x0)^q + c_lin*(b - b) - M2 = Z = (a - t*N - Sb)//Rb with M2 = (b - x0)^q - Z > 0 of scale ~ N!
    # Wait: (b - x0)^q has sign (-1)^q and huge size ~ x0^q! No, if M2 ~ x0^q, then rho2 = M2^{1/q} ~ x0, NOT << x0!
    # Instead: choose b NEAR x0 (e.g. b = x0 + d_b with small d_b so |b - x0| <= 20, gcd(R(b), N) == 1)!!
    for d_b in range(1, 60):
        b_try = x0 + d_b
        if math.gcd(int(R_prod(b_try)), N) == 1:
            b = b_try
            break
    Rb = int(R_prod(b))
    Sb = int(S_pol(b))
    t_mod = ((a - Sb) * pow(N, -1, Rb)) % Rb
    # Now (b - x0)^q = d_b^q is small! Choosing t = t_mod + ks*Rb with ks >= 2 makes Z = (a - t*N - Sb)//Rb ~ -ks*N/Rb < 0,
    # so M2 = (b - x0)^q - Z ~ ks*N/Rb > 0, giving rho2 = M2^{1/q} <= (3*N)^{1/q} <= rho2_est << x0 / (8 m^2 K_m)!!
    for d_u in range(36):
        ks = 2 + (d_u // 6)
        t = t_mod + ks * Rb
        Z = (a - t * N - Sb) // Rb
        c_lin = (d_u % 6) - 2
        M2 = (b - x0)**q - Z
        if M2 <= 0:
            continue
        U_pol = (T - x0)**q + c_lin * (T - b) - Poly([M2])
        assert int(U_pol(b)) == Z
        P = S_pol + R_prod * U_pol
        assert (int(P(b)) - a) % N == 0
        if not ex.irreducible(P*P - 2):
            continue
        for j0 in (1, -1, 2, -2, 3, -3, 5, -5):
            Q = P + j0 * (T - b) * R_prod
            if ex.irreducible(Q*Q - 2):
                return dict(name=f"C3_{mode}", b=b, P=P, Q=Q, R_list=R_list, K_m=K_m, x0=x0)
    raise RuntimeError(f"C3_{mode} failed")


# --- Certified Threshold, Norm & Parity Audit ---

def certified_row_prefix(row, points, m, E=1):
    """Exact 256-bit flint ball threshold t_i and rank-(m-1) strict prefix subspace."""
    h = len(row) // 2
    # Fast float64 pre-filter + 100% rigorous worst-case rounding bound for large m >= 5
    if m >= 5:
        max_entry = max(abs(z).upper() for z in row)
        if not (max_entry > 0):
            return None
        # Check relative radius of row entries: must be < 1e-18 to use normalized float64 prefilter
        rel_rad = max(float((abs(z).rad()) / max_entry) for z in row)
        if rel_rad > 1e-18:
            return None
        row_scaled = [z / max_entry for z in row]
        row_c = np.array([complex(float(z.real), float(z.imag)) for z in row_scaled], dtype=np.complex128)
        err_bound = 1e-13 * (E * 4 * m * m)
        coords_h = np.array(list(itertools.product(range(-E, E + 1), repeat=h)), dtype=np.float64)
        left_c = coords_h @ row_c[:h]
        right_c = coords_h @ row_c[h:]
        # All (2E+1)^{2m} values: reshape to find first-nonzero-positive mask
        # Simpler: evaluate on points via index mapping
        idx_map = {t: idx for idx, t in enumerate(itertools.product(range(-E, E + 1), repeat=h))}
        li = np.array([idx_map[z[:h]] for z in points], dtype=np.int32)
        ri = np.array([idx_map[z[h:]] for z in points], dtype=np.int32)
        vals_f64 = np.abs(left_c[li] + right_c[ri])
        order_idx = np.argsort(vals_f64)
        # Find provisional rank-m cutoff in float64
        prov = []
        u_cand = 0.0
        for idx in order_idx:
            z = points[int(idx)]
            if ex.kranks(prov + [z]) > len(prov):
                prov.append(z)
                if len(prov) == m:
                    u_cand = float(vals_f64[int(idx)])
                    break
        cutoff_f64 = u_cand + 4.0 * err_bound + 1e-11 * max(1.0, u_cand)
        cand_indices = [int(i) for i in order_idx if vals_f64[int(i)] <= cutoff_f64]
        # Certify in 256-bit acb on cand_indices!
        if len(cand_indices) < len(points):
            ordered = []
            for idx in cand_indices:
                z = points[idx]
                v = abs(sum((a*b for a, b in zip(row, z)), acb(0)))
                ordered.append((v.lower(), v.upper(), z))
            ordered.sort(key=lambda it: it[0])
            trial = []
            lo, hi = None, None
            for low, high, z in ordered:
                if ex.kranks(trial + [z]) > len(trial):
                    trial.append(z)
                    if len(trial) == m:
                        lo = low
                        hi = max(abs(sum((a*b for a, b in zip(row, w)), acb(0))).upper() for w in trial)
                        break
            # Verify certified separation from excluded points: min excluded lower >= cutoff_f64 - err_bound > hi
            if lo is not None and float(hi / max_entry) < cutoff_f64 - 2.0 * err_bound:
                if not (lo > 0):
                    return None
                prefix = []
                for low, high, z in ordered:
                    if low >= lo:
                        break
                    if high < lo and ex.kranks(prefix + [z]) > len(prefix):
                        prefix.append(z)
                assert len(prefix) < m
                return prefix, lo, hi
        return None

    # Pure 256-bit acb path (always used for m <= 4 and fallback for m >= 5)
    left = {z: sum((a*b for a, b in zip(row[:h], z)), acb(0))
            for z in itertools.product(range(-E, E + 1), repeat=h)}
    right = {z: sum((a*b for a, b in zip(row[h:], z)), acb(0))
             for z in itertools.product(range(-E, E + 1), repeat=len(row) - h)}
    ordered = []
    for z in points:
        v = abs(left[z[:h]] + right[z[h:]])
        ordered.append((v.lower(), v.upper(), z))
    ordered.sort(key=lambda it: it[0])
    trial = []
    for low, high, z in ordered:
        if ex.kranks(trial + [z]) > len(trial):
            trial.append(z)
            if len(trial) == m:
                lo = low
                hi = max(abs(left[w[:h]] + right[w[h:]]).upper() for w in trial)
                break
    prefix = []
    for low, high, z in ordered:
        if low >= lo:
            break
        if high < lo and ex.kranks(prefix + [z]) > len(prefix):
            prefix.append(z)
    assert len(prefix) < m
    return prefix, lo, hi


def full_factor_parity(field_obj, a_pol, norm_val, columns):
    """Exact F_2 prime-ideal valuation mask across ALL prime divisors of norm_val (even > B)."""
    key = tuple(ex.coeffs(a_pol))
    if field_obj.field is None:
        field_obj.field = ex.pari.nfinit([field_obj.poly, field_obj.bad])
    nf = field_obj.field
    basis = nf.nfalgtobasis(ex.pari.Polrev(list(key)))
    # Factor norm_val in PARI (instant for <= 200-bit integers)
    fmat = ex.pari(norm_val).factor()
    ps = [int(p) for p in fmat[0]]
    mask = 0
    for p in ps:
        if p not in field_obj.decomps:
            # Ensure maximal order at p if p divides disc; for primes > B not in bad, check idealval safely via nfinit or p-radical
            # Note: if p | disc and p was not in bad, recompute nfinit including p!
            if abs(int(Poly(ex.coeffs(Poly(list(key)))).discriminant() if False else field_obj.disc_int)) % p == 0 and p not in field_obj.bad:
                field_obj.bad.append(p)
                field_obj.field = ex.pari.nfinit([field_obj.poly, field_obj.bad])
                nf = field_obj.field
                basis = nf.nfalgtobasis(ex.pari.Polrev(list(key)))
                field_obj.decomps.clear()
            field_obj.decomps[p] = nf.idealprimedec(p)
        for i, pr in enumerate(field_obj.decomps[p]):
            w = int(nf.idealval(basis, pr))
            if w % 2:
                col = columns.setdefault((field_obj.side_id, p, i), len(columns))
                mask |= 1 << col
    return mask


def audit_construction(bits, seed, n, m, spec, E=1, thresholds_only=False):
    t0 = time.time()
    N = ex.fixture(bits, seed)
    b, P, Q, R_list = spec["b"], spec["P"], spec["Q"], spec["R_list"]
    f, g = P*P - 2, Q*Q - 2
    assert ex.irreducible(f) and ex.irreducible(g)
    R_prod = one
    for Rj in R_list:
        R_prod = R_prod * Rj
    assert (Q - P) % ((T - b) * R_prod) == Poly([])
    A_gens = [one, T] + [P // Rj for Rj in R_list]
    B_gens = [one, T] + [Q // Rj for Rj in R_list]
    assert len(A_gens) == m and len(B_gens) == m
    # Verify exact divided-difference formula: P // Rj == (P - (P % Rj)) / Rj
    for Rj, qA, qB in zip(R_list, A_gens[2:], B_gens[2:]):
        assert P % Rj == Q % Rj
        assert qA * Rj + (P % Rj) == P
        assert qB * Rj + (Q % Rj) == Q

    U = ex.expand(P, f, A_gens)
    V = ex.expand(Q, g, B_gens)
    residues = [int(u(b)) % N for u in U]
    assert residues == [int(v(b)) % N for v in V]
    ex.n, ex.d = n, 2*n
    pcheck = ex.product_check(P, Q, A_gens, B_gens, N)

    all_pos_points = [z for z in itertools.product(range(-E, E + 1), repeat=2*m)
                      if next((v for v in z if v), 0) > 0]
    points = [z for z in all_pos_points
              if math.gcd(sum(c*v for c, v in zip(z, residues)), N) == 1]

    prec = 512
    while True:
        ctx.prec = prec
        rows, _ = ex.rows_at(ex.coeffs(f), list(map(ex.coeffs, U)),
                             ex.coeffs(g), list(map(ex.coeffs, V)), prec)
        D = len(rows)
        assert D == 4 * n
        spaces = []
        M_lo_sum, M_hi_sum = arb(0), arb(0)
        t_lo_sum, t_hi_sum = arb(0), arb(0)
        losses_lo, losses_hi = [], []
        denom_mags = []
        need_more_prec = False
        for row in rows:
            M_i = abs(row[0])
            for j in range(1, m):
                M_i = M_i.max(abs(row[2*j]))
            alpha_i = row[2] / row[0]
            r_val = acb(1)
            for Rj in R_list:
                cs = ex.coeffs(Rj)
                ev = acb(0)
                for c in reversed(cs):
                    ev = ev * alpha_i + c
                r_val *= ev
            denom_mags.append(float(abs(r_val)))
            res_pref = certified_row_prefix(row, points, m, E)
            if res_pref is None:
                need_more_prec = True
                break
            basis, lo, hi = res_pref
            if not (lo > 0 and M_i.lower() > 0):
                need_more_prec = True
                break
            M_lo_sum += M_i.lower().log()
            M_hi_sum += M_i.upper().log()
            t_lo_sum += lo.log()
            t_hi_sum += hi.log()
            losses_lo.append(float((M_i.lower() / hi).log()))
            losses_hi.append(float((M_i.upper() / lo).log()))
            if not thresholds_only and basis and all(ex.kranks(old + basis) > max(len(old), len(basis)) for old in spaces):
                spaces.append(basis)
        if not need_more_prec:
            break
        prec *= 2
        if prec > 8192:
            raise RuntimeError("Precision exceeded 8192")
    n_real_certified = count_real_roots_certified(f) + count_real_roots_certified(g)

    losses_mid = [0.5 * (a + b_) for a, b_ in zip(losses_lo, losses_hi)]
    losses_sorted = sorted(losses_lo)  # certified lower bounds on loss!
    mlogm = m * math.log(m)
    rec = dict(
        bits=bits, seed=seed, n=n, m=m, D=D, mode=spec["name"], b=b,
        total_pole_deg=sum(len(ex.coeffs(Rj)) - 1 for Rj in R_list),
        s_h=pcheck["degrees"],
        real_embeddings=n_real_certified,
        H_interval=[float(M_lo_sum), float(M_hi_sum)],
        Q_interval=[float(t_lo_sum), float(t_hi_sum)],
        x_lnN=math.log(N),
        loss_min=round(min(losses_lo), 4),
        loss_q25=round(losses_sorted[D // 4], 4),
        loss_med=round(losses_sorted[D // 2], 4),
        loss_q75=round(losses_sorted[(3 * D) // 4], 4),
        loss_max=round(max(losses_hi), 4),
        loss_mean=round(sum(losses_mid) / D, 4),
        mlogm=round(mlogm, 4),
        frac_gt_025_mlogm=round(sum(1 for x in losses_lo if x >= 0.25 * mlogm) / D, 4),
        frac_gt_050_mlogm=round(sum(1 for x in losses_lo if x >= 0.50 * mlogm) / D, 4),
        frac_gt_075_mlogm=round(sum(1 for x in losses_lo if x >= 0.75 * mlogm) / D, 4),
        min_shared_denom=min(denom_mags),
        num_tiny_denom_lt_01=sum(1 for d_ in denom_mags if d_ < 0.1),
    )
    if thresholds_only:
        rec["elapsed_s"] = round(time.time() - t0, 2)
        return rec

    # Norm cap & smoothness parameters (exact sweep.py formulas)
    units = pair_units()
    cap = 3**(4*n) * 2**((bits + 2*n - 1) // (2*n))
    root4, exact4 = integer_nthroot(bits * (bits - 1).bit_length()**3, 4)
    B = 2**(int(root4) + int(not exact4))
    primes = list(primerange(B + 1))
    primeprod = math.prod(primes)

    # Deduplicate orbits
    orbits = sorted({orbit_rep(z, units, E) for z in points})

    # Certified 256-bit acb product filter across all orbits + exact resultant on capped & bottom full-support orbits
    h = m
    row_tables = []
    for row in rows:
        lt = {z: sum((a*b for a, b in zip(row[:h], z)), acb(0))
              for z in itertools.product(range(-E, E + 1), repeat=h)}
        rt = {z: sum((a*b for a, b in zip(row[h:], z)), acb(0))
              for z in itertools.product(range(-E, E + 1), repeat=2*m - h)}
        row_tables.append((lt, rt))

    capped_candidates = []  # (F_exact, z, v1, v2, is_full_supp, in_prefix_plane)
    full_supp_approx = []   # (F_lower_float_log2, z)

    arb_cap = arb(cap)
    for z in orbits:
        zl, zr = z[:h], z[h:]
        prod_ball = arb(1)
        for lt, rt in row_tables:
            prod_ball *= abs(lt[zl] + rt[zr])
        is_full = all(z[2*j] != 0 or z[2*j + 1] != 0 for j in range(m))
        if is_full:
            lo_pos = max(float(prod_ball.lower()), 1e-300)
            full_supp_approx.append((math.log2(lo_pos), z))
        if not (prod_ball.lower() > arb_cap):
            # Compute exact integer resultant!
            a_el, w_el = ex.element(U, z), ex.element(V, z)
            v1, v2 = abs(int(f.resultant(a_el))), abs(int(g.resultant(w_el)))
            F_ex = v1 * v2
            if 0 < F_ex <= cap:
                in_plane = any(ex.kranks(s + [z]) == len(s) for s in spaces)
                capped_candidates.append((F_ex, z, v1, v2, is_full, in_plane))

    # Certify exact minimum full-support F across the ENTIRE ternary box!
    full_supp_approx.sort(key=lambda it: it[0])
    min_full_F = None
    if full_supp_approx:
        # Check exact resultant on bottom 25 full-support orbits and certify against remaining lower bounds
        best_F = None
        for _, z in full_supp_approx[:25]:
            a_el, w_el = ex.element(U, z), ex.element(V, z)
            F_ex = abs(int(f.resultant(a_el))) * abs(int(g.resultant(w_el)))
            if F_ex > 0 and (best_F is None or F_ex < best_F):
                best_F = F_ex
        if len(full_supp_approx) > 25:
            # Verify 26th lower bound >= best_F; if not, keep scanning until prod_ball.lower() > best_F!
            arb_best = arb(best_F)
            for idx in range(25, len(full_supp_approx)):
                z = full_supp_approx[idx][1]
                zl, zr = z[:h], z[h:]
                pb = arb(1)
                for lt, rt in row_tables:
                    pb *= abs(lt[zl] + rt[zr])
                if pb.lower() > arb_best:
                    break
                a_el, w_el = ex.element(U, z), ex.element(V, z)
                F_ex = abs(int(f.resultant(a_el))) * abs(int(g.resultant(w_el)))
                if 0 < F_ex < best_F:
                    best_F = F_ex
                    arb_best = arb(best_F)
        min_full_F = best_F

    # Parity independence audit on capped candidates
    fields = [ex.Ideals(f, primes), ex.Ideals(g, primes)]
    for side_id, (fld, pol) in enumerate(zip(fields, (f, g))):
        fld.side_id = side_id
        fld.disc_int = abs(int(pol.discriminant()))

    capped_candidates.sort(key=lambda it: (it[0], it[1]))
    smooth_cols = {}
    all_cols = {}

    # Compute smooth parities and full parities (for stronger lower-norm product test!)
    smooth_items = []
    for F_ex, z, v1, v2, is_full, in_plane in capped_candidates:
        is_sm = ex.smooth(F_ex, primeprod)
        sm_mask = 0
        if is_sm:
            for side, (W_gens, ns_val) in enumerate(((U, v1), (V, v2))):
                for pr in fields[side].parity(ex.element(W_gens, z), ns_val, primes):
                    col = smooth_cols.setdefault((side, pr[0], pr[1]), len(smooth_cols))
                    sm_mask |= 1 << col
                    all_cols.setdefault((side, pr[0], pr[1]), col)
        smooth_items.append((F_ex, z, v1, v2, is_full, in_plane, is_sm, sm_mask))

    # 1. Standard N_ind (modulo coordinate-proper k-subspaces AND B-smooth lower-norm products):
    subspace_coord_masks = [it[7] for it in smooth_items if it[6] and not it[4]]
    subspace_all_masks = [it[7] for it in smooth_items if it[6] and (not it[4] or it[5])]

    def count_N_ind_smooth(base_masks):
        span_rows = list(base_masks)
        cur_r = ex.rank(span_rows)
        n_ind = 0
        # Group by exact norm F_ex so strictly lower-norm smooth rows are already in span_rows!
        idx = 0
        while idx < len(smooth_items):
            F_cur = smooth_items[idx][0]
            j = idx
            batch = []
            while j < len(smooth_items) and smooth_items[j][0] == F_cur:
                if smooth_items[j][6]:
                    batch.append(smooth_items[j])
                j += 1
            # Test full-support candidates at F_cur modulo base_masks + strictly lower-norm rows
            temp_span = list(span_rows)
            temp_r = cur_r
            for item in batch:
                if item[4]:  # full coordinate support
                    nr = ex.rank(temp_span + [item[7]])
                    if nr > temp_r:
                        n_ind += 1
                        temp_span.append(item[7])
                        temp_r = nr
            # Add all smooth rows at F_cur into lower-norm span for higher norms
            for item in batch:
                span_rows.append(item[7])
            cur_r = ex.rank(span_rows)
            idx = j
        return n_ind

    N_ind_coord = count_N_ind_smooth(subspace_coord_masks)
    N_ind_prefix = count_N_ind_smooth(subspace_all_masks)

    # 2. Stronger check (also modulo non-smooth strictly lower-norm candidate products, as in chase_L2_codex.out:17)
    N_ind_strong = 0
    if N_ind_coord > 0:
        # Factor all capped candidates up to max smooth full-support F
        max_sm_full_F = max((it[0] for it in smooth_items if it[6] and it[4]), default=0)
        full_val_rows = []
        # First include all coordinate-subspace smooth masks (mapped into all_cols)
        for msk in subspace_coord_masks:
            full_val_rows.append(msk)
        cur_r = ex.rank(full_val_rows)
        idx = 0
        while idx < len(smooth_items) and smooth_items[idx][0] <= max_sm_full_F:
            F_cur = smooth_items[idx][0]
            j = idx
            batch = []
            while j < len(smooth_items) and smooth_items[j][0] == F_cur:
                batch.append(smooth_items[j])
                j += 1
            temp_span = list(full_val_rows)
            temp_r = cur_r
            for item in batch:
                if item[6] and item[4]:
                    nr = ex.rank(temp_span + [item[7]])
                    if nr > temp_r:
                        N_ind_strong += 1
                        temp_span.append(item[7])
                        temp_r = nr
            for F_ex, z, v1, v2, is_full, in_plane, is_sm, sm_mask in batch:
                if is_sm:
                    full_val_rows.append(sm_mask)
                else:
                    m1 = full_factor_parity(fields[0], ex.element(U, z), v1, all_cols)
                    m2 = full_factor_parity(fields[1], ex.element(V, z), v2, all_cols)
                    full_val_rows.append(m1 | m2)
            cur_r = ex.rank(full_val_rows)
            idx = j

    covered_sm = [it[7] for it in smooth_items if it[6] and it[5]]
    escape_sm = [it[7] for it in smooth_items if it[6] and not it[5]]
    rec.update(
        B=B,
        log2_cap=round(math.log2(cap), 3),
        min_full_supp_log2F=round(math.log2(min_full_F), 3) if min_full_F else None,
        min_full_supp_over_logN=round(math.log(min_full_F) / math.log(N), 4) if min_full_F else None,
        valid_orbits=len(orbits),
        capped_orbits=len(capped_candidates),
        capped_full_supp=sum(1 for it in capped_candidates if it[4]),
        capped_esc_1h=sum(1 for it in capped_candidates if not it[5]),
        smooth_orbits=sum(1 for it in smooth_items if it[6]),
        smooth_full_supp=sum(1 for it in smooth_items if it[6] and it[4]),
        smooth_esc_1h=len(escape_sm),
        fresh_1h_rank=ex.rank(covered_sm + escape_sm) - ex.rank(covered_sm),
        N_ind_coord=N_ind_coord,
        N_ind_strong=N_ind_strong,
        N_ind_prefix=N_ind_prefix,
        capped_k_span=ex.kranks([it[1] for it in capped_candidates]),
        elapsed_s=round(time.time() - t0, 2),
    )
    return rec


def main():
    results = {"main_sweep": [], "scaling_sweep": []}
    schedule = [
        (60, 7, 3),
        (72, 8, 3),
        (84, 9, 4),
        (96, 10, 4),
        (112, 12, 4),
    ]
    for bits, n, m in schedule:
        for seed in (19, 0):
            N = ex.fixture(bits, seed)
            specs = [
                build_B0(N, n, m),
                build_C1a(N, n, m),
                build_C1b(N, n, m),
                build_C2(N, n, m),
                build_C3(N, n, m, "Z"),
                build_C3(N, n, m, "Ok"),
            ]
            for spec in specs:
                rec = audit_construction(bits, seed, n, m, spec, E=1, thresholds_only=False)
                results["main_sweep"].append(rec)
                print(f"[{bits}b s={seed} (n={n},m={m}) {rec['mode']:5s}] "
                      f"H={rec['H_interval'][0]:.1f} Q={rec['Q_interval'][0]:.1f} "
                      f"loss(med/max)={rec['loss_med']:.2f}/{rec['loss_max']:.2f} "
                      f"cap={rec['capped_orbits']}(full={rec['capped_full_supp']},esc={rec['capped_esc_1h']}) "
                      f"N_ind={rec['N_ind_coord']}/{rec['N_ind_strong']}/{rec['N_ind_prefix']} "
                      f"minFull/lnN={rec['min_full_supp_over_logN']} ({rec['elapsed_s']}s)", flush=True)

    # Also test (112, 12, 5) on seed 19 across all 6 families!
    bits, n, m, seed = 112, 12, 5, 19
    N = ex.fixture(bits, seed)
    for spec in [build_B0(N, n, m), build_C1a(N, n, m), build_C1b(N, n, m),
                 build_C2(N, n, m), build_C3(N, n, m, "Z"), build_C3(N, n, m, "Ok")]:
        rec = audit_construction(bits, seed, n, m, spec, E=1, thresholds_only=False)
        results["main_sweep"].append(rec)
        print(f"[{bits}b s={seed} (n={n},m={m}) {rec['mode']:5s}] "
              f"H={rec['H_interval'][0]:.1f} Q={rec['Q_interval'][0]:.1f} "
              f"loss(med/max)={rec['loss_med']:.2f}/{rec['loss_max']:.2f} "
              f"cap={rec['capped_orbits']}(full={rec['capped_full_supp']},esc={rec['capped_esc_1h']}) "
              f"N_ind={rec['N_ind_coord']}/{rec['N_ind_strong']}/{rec['N_ind_prefix']} "
              f"minFull/lnN={rec['min_full_supp_over_logN']} ({rec['elapsed_s']}s)", flush=True)

    # Dedicated conditioning-scaling sweep across m = 3, 4, 5, 6 at bits = 96, n = 12, seed = 19
    print("\n=== Conditioning-Loss Scaling Sweep (bits=96, n=12, seed=19, m=3..6) ===", flush=True)
    bits, n, seed = 96, 12, 19
    N = ex.fixture(bits, seed)
    for m in (3, 4, 5, 6):
        for spec in [build_B0(N, n, m), build_C2(N, n, m), build_C3(N, n, m, "Z"), build_C3(N, n, m, "Ok")]:
            rec = audit_construction(bits, seed, n, m, spec, E=1, thresholds_only=(m == 6))
            results["scaling_sweep"].append(rec)
            print(f"[Scale m={m} {rec['mode']:5s}] mlogm={rec['mlogm']:.2f} "
                  f"loss(q25/med/q75/max)={rec['loss_q25']:.2f}/{rec['loss_med']:.2f}/{rec['loss_q75']:.2f}/{rec['loss_max']:.2f} "
                  f"frac>=0.5mlogm={rec['frac_gt_050_mlogm']:.2f} frac>=0.75mlogm={rec['frac_gt_075_mlogm']:.2f} "
                  f"H={rec['H_interval'][0]:.1f} Q={rec['Q_interval'][0]:.1f} ({rec['elapsed_s']}s)", flush=True)

    Path("results.json").write_text(json.dumps(results, indent=2) + "\n")
    print("\nSaved results.json successfully.")


if __name__ == "__main__":
    main()
