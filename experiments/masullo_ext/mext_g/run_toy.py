"""Pre-registered diagnostic toy checks for masullo_ext (mext_g).
Strictly read-only with respect to the main repo; writes results.json in cwd."""
from __future__ import annotations

import itertools
import json
import math
import random
from functools import reduce
from pathlib import Path

from flint import acb, arb, ctx, fmpq_mat, fmpz_poly as Poly
from sympy import integer_nthroot, isprime, nextprime, sqrt_mod

HERE = Path(__file__).resolve().parent


def predicted_P(m: int, D_rel: int) -> int:
    return m * (m - 1) * math.factorial(D_rel - 2) // math.factorial(D_rel - m)


def poly_eval_acb(coeffs: list[int], z: acb) -> acb:
    s = acb(0)
    for c in reversed(coeffs):
        s = s * z + int(c)
    return s


def det_acb(M: list[list[acb]]) -> acb:
    n = len(M)
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    if n == 3:
        return (
            M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])
        )
    s = acb(0)
    for j in range(n):
        minor = [row[:j] + row[j + 1 :] for row in M[1:]]
        term = M[0][j] * det_acb(minor)
        s = s + term if j % 2 == 0 else s - term
    return s


def irreducible(f: Poly) -> bool:
    facs = f.factor()[1]
    return len(facs) == 1 and facs[0][1] == 1


def make_qsqrt2_semiprime(bits: int, seed: int) -> tuple[int, int, int, int]:
    rng = random.Random(91001 * bits + seed)
    lo = 1 << (bits // 2 - 1)
    hi = (1 << (bits // 2)) - 1
    while True:
        p = int(nextprime(rng.randrange(lo, hi)))
        if p % 8 not in (1, 7):
            continue
        q = int(nextprime(rng.randrange(lo, hi)))
        if q == p or q % 8 not in (1, 7):
            continue
        N = p * q
        if N.bit_length() != bits:
            continue
        rp = int(sqrt_mod(2, p))
        rq = int(sqrt_mod(2, q))
        # CRT lift for sqrt(2) mod N
        inv_p = pow(p, -1, q)
        A = (rp + p * ((rq - rp) * inv_p % q)) % N
        if A < 10:
            continue
        assert (A * A - 2) % N == 0
        return N, p, q, A


def make_qi_semiprime(bits: int, seed: int) -> tuple[int, int, int, int]:
    rng = random.Random(77003 * bits + seed)
    lo = 1 << (bits // 2 - 1)
    hi = (1 << (bits // 2)) - 1
    while True:
        p = int(nextprime(rng.randrange(lo, hi)))
        if p % 4 != 1:
            continue
        q = int(nextprime(rng.randrange(lo, hi)))
        if q == p or q % 4 != 1:
            continue
        N = p * q
        if N.bit_length() != bits:
            continue
        rp = int(sqrt_mod(-1, p))
        rq = int(sqrt_mod(-1, q))
        inv_p = pow(p, -1, q)
        A = (rp + p * ((rq - rp) * inv_p % q)) % N
        assert (A * A + 1) % N == 0
        return N, p, q, A


def make_qcbrt3_semiprime(bits: int, seed: int) -> tuple[int, int, int, int]:
    rng = random.Random(55007 * bits + seed)
    lo = 1 << (bits // 2 - 1)
    hi = (1 << (bits // 2)) - 1
    while True:
        p = int(nextprime(rng.randrange(lo, hi)))
        # p = 2 mod 3 gives a unique cube root of 3 mod p!
        if p % 3 != 2:
            continue
        q = int(nextprime(rng.randrange(lo, hi)))
        if q == p or q % 3 != 2:
            continue
        N = p * q
        if N % 4 != 3 or N.bit_length() != bits:
            continue
        rp = pow(3, (2 * p - 1) // 3, p)
        rq = pow(3, (2 * q - 1) // 3, q)
        assert (rp**3 - 3) % p == 0 and (rq**3 - 3) % q == 0
        inv_p = pow(p, -1, q)
        B = (rp + p * ((rq - rp) * inv_p % q)) % N
        assert (B**3 - 3) % N == 0
        return N, p, q, B


def build_divdiff_pair(N: int, A: int, n: int, d_val: int) -> tuple[int, Poly, Poly, list[Poly], list[Poly], int]:
    """Build monic degree-n P, Q with P(b)==Q(b)==A, P^2 - d_val and Q^2 - d_val irreducible,
    and divided-difference non-power generators."""
    T = Poly([0, 1])
    one = Poly([1])
    base = max(3, int(integer_nthroot(A, n)[0]))
    for b in range(base, base + 80):
        rem = A - b**n
        dig = []
        ok = True
        for _ in range(n - 1):
            digit = (rem + b // 2) % b - b // 2
            dig.append(digit)
            rem = (rem - digit) // b
        P = Poly(dig + [rem, 1])
        assert int(P(b)) == A
        f = P * P - d_val
        if not irreducible(f):
            continue
        for j in (1, -1, 2, -2, 3, -3, 5, -5):
            if n == 3:
                Q = P + j * (T - b) * (T - 1)
                gen_P = [one, T, P // (T - 1)]
                gen_Q = [one, T, Q // (T - 1)]
            else:
                Q = P + j * (T - b) * (T - 1) * (T + 1)
                gen_P = [one, T, P // (T - 1), P // (T + 1)]
                gen_Q = [one, T, Q // (T - 1), Q // (T + 1)]
            g = Q * Q - d_val
            if Q.degree() == n and int(Q[n]) == 1 and irreducible(g):
                assert int(Q(b)) == A
                for uP, uQ in zip(gen_P, gen_Q):
                    assert int(uP(b)) == int(uQ(b))
                return b, P, Q, gen_P, gen_Q, j
    raise RuntimeError(f"Could not build pair for N={N}, n={n}, d_val={d_val}")


def certify_quad_base_deltas(
    P: Poly, Q: Poly, gen_P: list[Poly], gen_Q: list[Poly], d_val: int
) -> dict:
    """Certify absolute Delta_F (sampling check + exact block-rank proof + full check when D<=12)
    and exact relative Q2 = N_{k/Q}(delta_rel) via certified acb balls."""
    f = P * P - d_val
    g = Q * Q - d_val
    m = len(gen_P)
    r = 2 * m
    P_c = [int(c) for c in P]
    Q_c = [int(c) for c in Q]
    gP_c = [[int(c) for c in u] for u in gen_P]
    gQ_c = [[int(c) for c in u] for u in gen_Q]

    for prec in (512, 1024, 2048, 4096, 8192, 16384, 32768, 65536):
        ctx.prec = prec
        w_pos = acb(d_val).sqrt()
        w_neg = -w_pos
        f_roots = [z for z, mult in f.complex_roots()]
        g_roots = [z for z, mult in g.complex_roots()]
        n = len(P_c) - 1

        # Partition roots into sigma_+ (P(z) == +w) and sigma_- (P(z) == -w)
        f_pos, f_neg = [], []
        ambig = False
        for z in f_roots:
            vz = poly_eval_acb(P_c, z)
            in_p = (vz - w_pos).contains(0)
            in_n = (vz - w_neg).contains(0)
            if in_p == in_n:
                ambig = True
                break
            (f_pos if in_p else f_neg).append(z)
        g_pos, g_neg = [], []
        for z in g_roots:
            vz = poly_eval_acb(Q_c, z)
            in_p = (vz - w_pos).contains(0)
            in_n = (vz - w_neg).contains(0)
            if in_p == in_n:
                ambig = True
                break
            (g_pos if in_p else g_neg).append(z)
        if ambig or len(f_pos) != n or len(g_pos) != n:
            continue

        rows_pos = [[poly_eval_acb(u, z) for u in gP_c] for z in f_pos] + [
            [poly_eval_acb(u, z) for u in gQ_c] for z in g_pos
        ]
        rows_neg = [[poly_eval_acb(u, z) for u in gP_c] for z in f_neg] + [
            [poly_eval_acb(u, z) for u in gQ_c] for z in g_neg
        ]

        prod_pos = acb(1)
        for S in itertools.combinations(range(2 * n), m):
            prod_pos *= det_acb([rows_pos[i] for i in S])
        prod_neg = acb(1)
        for S in itertools.combinations(range(2 * n), m):
            prod_neg *= det_acb([rows_neg[i] for i in S])

        q2 = (prod_pos * prod_neg) ** 2
        z_int = q2.real.unique_fmpz()
        im_ok = q2.imag.contains(0)

        # Also certify absolute Delta_F == 0 by constructing an explicit (m+1, m-1) absolute r-minor
        # In Z-basis (u_1, w*u_1, ..., u_m, w*u_m), row at z_pos has w = +w_pos, row at z_neg has w = -w_pos
        abs_rows_pos = []
        for row in rows_pos:
            ar = []
            for val in row:
                ar.extend([val, w_pos * val])
            abs_rows_pos.append(ar)
        abs_rows_neg = []
        for row in rows_neg:
            ar = []
            for val in row:
                ar.extend([val, w_neg * val])
            abs_rows_neg.append(ar)
        # Take m+1 rows from pos and m-1 rows from neg (valid since 2n > m => 2n >= m+1)
        test_abs_minor = det_acb(abs_rows_pos[: m + 1] + abs_rows_neg[: m - 1])
        abs_zero_ok = test_abs_minor.contains(0)
        # If D = 4n <= 12, also check the full product over all (4n choose r) minors!
        if 4 * n <= 12:
            all_abs = abs_rows_pos + abs_rows_neg
            full_abs = acb(1)
            for S in itertools.combinations(range(4 * n), r):
                full_abs *= det_acb([all_abs[i] for i in S])
            abs_zero_ok = abs_zero_ok and full_abs.contains(0)

        # Generator envelope H = sum_{all 4n rows} log max_j |abs_row_j|
        H_val = 0.0
        for ar in abs_rows_pos + abs_rows_neg:
            mx = max(float(abs(x).mid()) for x in ar)
            H_val += math.log(max(mx, 1e-300))

        if z_int is not None and im_ok and abs_zero_ok:
            return {
                "Q2": int(z_int),
                "abs_Delta_F_zero": True,
                "prec": prec,
                "H": H_val,
            }
    raise RuntimeError("Precision exhausted in certify_quad_base_deltas")


def certify_cubic_c_example(N: int) -> dict:
    """Consultant C's Q(cbrt(3)) norms-1-and-8 pair f=T^2-(N+5)T/2+1, g=T^2-(N+3)T+2
    extended to m=3 < 2n=4 via w_1=(1,1), w_2=(s,s'), w_3=(s-2,0)."""
    assert N % 2 == 1
    f = Poly([1, -(N + 5) // 2, 1])
    g = Poly([2, -(N + 3), 1])
    m, n, e = 3, 2, 3
    for prec in (512, 1024, 2048):
        ctx.prec = prec
        s_roots = [z for z, _ in f.complex_roots()]
        sp_roots = [z for z, _ in g.complex_roots()]
        # Relative rows over k (identical across all 3 embeddings of k since f, g have Z coeffs!):
        # Side 1: (1, s_i, s_i - 2); Side 2: (1, sp_j, 0)
        rel_rows = [[acb(1), z, z - 2] for z in s_roots] + [
            [acb(1), z, acb(0)] for z in sp_roots
        ]
        prod_one_emb = acb(1)
        for S in itertools.combinations(range(4), 3):
            prod_one_emb *= det_acb([rel_rows[i] for i in S])
        # Across all e=3 embeddings of k=Q(cbrt(3)), prod is prod_one_emb^3, so Q2 = (prod_one_emb^3)^2 = prod_one_emb^6
        q2 = prod_one_emb**6
        z_int = q2.real.unique_fmpz()
        im_ok = q2.imag.contains(0)

        # Absolute rows in Z-basis {1, theta, theta^2} x {w_1, w_2, w_3} (r = 9, D = 12):
        # Take 4 rows from embedding sigma_1 of k and 5 rows from sigma_2, sigma_3: all 4 rows at sigma_1
        # lie in a 3D subspace of C^9, so any 9x9 minor picking all 4 rows at sigma_1 vanishes!
        th_roots = [z for z, _ in Poly([-3, 0, 0, 1]).complex_roots()]
        abs_blocks = []
        for th in th_roots:
            blk = []
            for rrow in rel_rows:
                ar = []
                for val in rrow:
                    ar.extend([val, th * val, th * th * val])
                blk.append(ar)
            abs_blocks.append(blk)
        # 4 rows from block 0, 3 from block 1, 2 from block 2 => rank <= 3 + 3 + 2 = 8 < 9 => det == 0
        sample_abs = det_acb(abs_blocks[0][:4] + abs_blocks[1][:3] + abs_blocks[2][:2])
        if z_int is not None and im_ok and sample_abs.contains(0):
            return {"Q2": int(z_int), "abs_Delta_F_zero": True, "prec": prec}
    raise RuntimeError("Precision exhausted in certify_cubic_c_example")


def vp(val: int, p: int) -> int:
    val = abs(val)
    if val == 0:
        return 10**9
    c = 0
    while val % p == 0:
        val //= p
        c += 1
    return c


def qsqrt2_orbit_key(coords: tuple[int, ...]) -> tuple[int, ...]:
    """Canonical O_k^*-orbit representative for k = Q(sqrt(2)), O_k^* = {+-1} * <1 + sqrt(2)>.
    coords = (a_1, b_1, ..., a_m, b_m) representing alpha_j = a_j + b_j*sqrt(2)."""
    m = len(coords) // 2
    pairs = [(coords[2 * j], coords[2 * j + 1]) for j in range(m)]
    # Find first nonzero component alpha_p = a + b*sqrt(2)
    idx = next(j for j, (a, b) in enumerate(pairs) if a != 0 or b != 0)
    # Use log|alpha_idx|_{sigma_+} - 0.5*log|N(alpha_idx)| in [0, log(1+sqrt(2)))!
    # Multiplying by (1+sqrt(2)) shifts log|alpha|_{sigma_+} - 0.5*log|N(alpha)| by +log(1+sqrt(2)).
    # Do exact integer unit steps until 1 <= |a + b*sqrt(2)| / sqrt(|a^2 - 2*b^2|) < 1 + sqrt(2),
    # i.e., |a^2 - 2*b^2| <= (a + b*sqrt(2))^2 < (3 + 2*sqrt(2))*|a^2 - 2*b^2|!
    ctx.prec = 256
    s2 = arb(2).sqrt()
    u_log = (arb(1) + s2).log()
    a0, b0 = pairs[idx]
    norm0 = abs(a0 * a0 - 2 * b0 * b0)
    cur = abs(arb(a0) + arb(b0) * s2).log() - arb(norm0).log() / 2
    t = -math.floor(float((cur / u_log).mid()))

    def mul_unit_pow(pairs_in: list[tuple[int, int]], k_step: int) -> list[tuple[int, int]]:
        out = list(pairs_in)
        if k_step > 0:
            for _ in range(k_step):
                out = [(a + 2 * b, a + b) for (a, b) in out]
        elif k_step < 0:
            for _ in range(-k_step):
                out = [(-a + 2 * b, a - b) for (a, b) in out]
        return out

    red = mul_unit_pow(pairs, t)
    # Center half-open fundamental domain (-u_log/2, +u_log/2] so unit shifts k*u_log never hit endpoints
    u_f = float(u_log.mid())
    while True:
        ra, rb = red[idx]
        r_log = float((abs(arb(ra) + arb(rb) * s2).log() - arb(norm0).log() / 2).mid())
        if r_log <= -0.5 * u_f:
            red = mul_unit_pow(red, 1)
        elif r_log > 0.5 * u_f:
            red = mul_unit_pow(red, -1)
        else:
            break
    ra, rb = red[idx]
    if arb(ra) + arb(rb) * s2 < 0:
        red = [(-a, -b) for (a, b) in red]
    return tuple(x for ab in red for x in ab)


def qsqrt2_k_rank_is_full(coords: tuple[int, ...], span_check_vecs: list[tuple[int, ...]]) -> int:
    """Compute dim_k span_k(span_check_vecs) for vectors in O_k^m (k=Q(sqrt(2)))."""
    if not span_check_vecs:
        return 0
    m = len(coords) // 2
    rows = []
    for v in span_check_vecs:
        # v defines two Z-independent rows in Q^(2m): v and sqrt(2)*v
        r1 = list(v)
        r2 = []
        for j in range(m):
            a, b = v[2 * j], v[2 * j + 1]
            r2.extend([2 * b, a])
        rows.append(r1)
        rows.append(r2)
    _, rnk = fmpq_mat(rows).rref()
    assert rnk % 2 == 0
    return rnk // 2


def analyze_qsqrt2_box(
    N: int, P: Poly, Q: Poly, gen_P: list[Poly], gen_Q: list[Poly], E: int, H_val: float
) -> dict:
    m = len(gen_P)
    n = P.degree()
    d = 2 * n
    D = 4 * n
    r = 2 * m
    f = P * P - 2
    g = Q * Q - 2

    # Basis over Z: [u_1, P*u_1 mod f, ..., u_m, P*u_m mod f]
    bas_K = [v for u in gen_P for v in (u, (u * P) % f)]
    bas_Kp = [v for u in gen_Q for v in (u, (u * Q) % g)]

    # Precompute roots at sigma_+ and sigma_- for exact Masullo A.2 cone assignment
    ctx.prec = 256
    s2 = acb(2).sqrt()
    P_c = [int(c) for c in P]
    Q_c = [int(c) for c in Q]
    f_pos = [z for z, _ in f.complex_roots() if (poly_eval_acb(P_c, z) - s2).contains(0)]
    g_pos = [z for z, _ in g.complex_roots() if (poly_eval_acb(Q_c, z) - s2).contains(0)]
    bK_c = [[int(c) for c in u] for u in bas_K]
    bKp_c = [[int(c) for c in u] for u in bas_Kp]

    Y0 = ((2 * E + 1) ** D) * math.ceil(N ** (1.0 / d))
    u_log = math.log(1.0 + math.sqrt(2.0))
    R_unit_bound = H_val / (2 * n) + 2 * math.log(r * E) + math.log(Y0) / (2 * n) + u_log

    raw_nonzero = 0
    prim_z_count = 0
    orbits = {}
    cones_seen = set()
    max_t_ulog = 0.0

    rng_coords = range(-E, E + 1)
    for z in itertools.product(rng_coords, repeat=r):
        if not any(z):
            continue
        # Canonical sign: first nonzero > 0 to avoid trivial +-1 double count in raw scan
        first_nz = next(x for x in z if x != 0)
        if first_nz < 0:
            continue

        # Construct algebraic elements A_K(T), A_Kp(T)
        AK = Poly([0])
        AKp = Poly([0])
        for c_z, uK, uKp in zip(z, bas_K, bas_Kp):
            if c_z:
                AK = AK + c_z * uK
                AKp = AKp + c_z * uKp
        if AK == 0 or AKp == 0:
            continue
        W1 = int(f.resultant(AK))
        if W1 == 0 or abs(W1) > Y0:
            continue
        W2 = int(g.resultant(AKp))
        if W2 == 0:
            continue
        F_val = abs(W1 * W2)
        if F_val > Y0:
            continue

        raw_nonzero += 1
        g_all = reduce(math.gcd, (abs(x) for x in z))
        if g_all == 1:
            prim_z_count += 1

        # Compute |F_rel(z)|_{sigma_+} via roots in f_pos, g_pos
        prod_pos = acb(1)
        for rt in f_pos:
            s_val = acb(0)
            for c_z, u_c in zip(z, bK_c):
                if c_z:
                    s_val += c_z * poly_eval_acb(u_c, rt)
            prod_pos *= s_val
        for rt in g_pos:
            s_val = acb(0)
            for c_z, u_c in zip(z, bKp_c):
                if c_z:
                    s_val += c_z * poly_eval_acb(u_c, rt)
            prod_pos *= s_val

        log_F_pos = float(abs(prod_pos).log().mid())
        v_shift = log_F_pos / (2 * n) - math.log(F_val) / (4 * n)
        t_cone = -math.floor(v_shift / u_log)
        cones_seen.add(t_cone)
        max_t_ulog = max(max_t_ulog, abs(t_cone) * u_log)

        # Check Z-primitivity vs O_k-primitivity & O_k^*-orbit key
        okey = qsqrt2_orbit_key(z)
        if okey not in orbits or F_val < orbits[okey][0]:
            orbits[okey] = (F_val, z)

    # Check how many distinct O_k^*-orbits lie outside proper coordinate sub-modules
    # and compute the k-dimension spanned by all surviving small-norm orbits!
    rep_vecs = [z for (_, z) in orbits.values()]
    spanned_k_dim = qsqrt2_k_rank_is_full((0,) * r, rep_vecs)
    full_support_orbits = sum(
        1
        for (_, z) in orbits.values()
        if all(z[2 * j] != 0 or z[2 * j + 1] != 0 for j in range(m))
    )
    # Even stricter: orbits NOT contained in any proper coordinate subspace {alpha_j = 0}
    # vs span_k dimension of ALL orbits <= Y0: if spanned_k_dim < m, ALL orbits lie in a single proper k-subspace!
    # Also check orbits outside the span of lower-coordinate subspaces (alpha_m == 0)
    sub_vecs = [z for (_, z) in orbits.values() if z[2 * (m - 1)] == 0 and z[2 * m - 1] == 0]
    outside_last_coord_subspace = 0
    for _, z in orbits.values():
        if z[2 * (m - 1)] != 0 or z[2 * m - 1] != 0:
            outside_last_coord_subspace += 1

    return {
        "Y0": int(Y0),
        "box_half_pop": ((2 * E + 1) ** r - 1) // 2,
        "raw_up_to_pm1": raw_nonzero,
        "z_primitive_up_to_pm1": prim_z_count,
        "distinct_Ok_unit_orbits": len(orbits),
        "occupied_A2_cones": len(cones_seen),
        "max_t_ulog": round(max_t_ulog, 4),
        "analytic_R_unit_bound": round(R_unit_bound, 4),
        "R_unit_bound_holds": max_t_ulog <= R_unit_bound + 1e-9,
        "spanned_k_dim_of_all_survivors": spanned_k_dim,
        "orbits_with_last_generator_nonzero": outside_last_coord_subspace,
        "full_support_orbits": full_support_orbits,
        "B3_sqrt_Y0_floor": round(math.sqrt(Y0), 1),
    }


def analyze_cubic_c_box(N: int, E: int = 2, Y_cap: int = 64) -> dict:
    """Scan [-E, E]^9 on Consultant C's Q(cbrt(3)) m=3, n=2 module w_1=(1,1), w_2=(s,s'), w_3=(s-2,0)."""
    f = Poly([1, -(N + 5) // 2, 1])
    g = Poly([2, -(N + 3), 1])
    # Represent elements of Q(cbrt(3))[T]/(f) via Z[theta]-combinations:
    # Actually N_{K/Q}(alpha_1*1 + alpha_2*s + alpha_3*(s-2)) = N_{k/Q}(N_{K/k}((alpha_1 - 2*alpha_3) + (alpha_2 + alpha_3)*s))!
    # Since f(T) = T^2 - c1*T + 1 (c1 = (N+5)//2), N_{K/k}(U + V*s) = U^2 + c1*U*V + V^2 in O_k = Z[theta] (theta^3 = 3)!
    # And g(T) = T^2 - c2*T + 2 (c2 = N+3), N_{K'/k}(U' + V'*s') = U'^2 + c2*U'*V' + 2*V'^2 in Z[theta]!
    c1 = (N + 5) // 2
    c2 = N + 3

    def mul_ok(x: tuple[int, int, int], y: tuple[int, int, int]) -> tuple[int, int, int]:
        a0, a1, a2 = x
        b0, b1, b2 = y
        return (
            a0 * b0 + 3 * (a1 * b2 + a2 * b1),
            a0 * b1 + a1 * b0 + 3 * a2 * b2,
            a0 * b2 + a1 * b1 + a2 * b0,
        )

    def add_ok(x: tuple[int, int, int], y: tuple[int, int, int], scale: int = 1) -> tuple[int, int, int]:
        return (x[0] + scale * y[0], x[1] + scale * y[1], x[2] + scale * y[2])

    def norm_ok(x: tuple[int, int, int]) -> int:
        a, b, c = x
        # N_{Q(cbrt(3))/Q}(a + b*theta + c*theta^2) = a^3 + 3*b^3 + 9*c^3 - 9*a*b*c
        return a**3 + 3 * b**3 + 9 * c**3 - 9 * a * b * c

    # Canonical O_k^*-orbit reduction using fundamental unit eps = theta^2 - 2 = (-2, 0, 1), eps^-1 = (4, 3, 2)
    eps = (-2, 0, 1)
    eps_inv = (4, 3, 2)
    assert mul_ok(eps, eps_inv) == (1, 0, 0)
    ctx.prec = 256
    th_arb = arb(3) ** (arb(1) / 3)
    eps_re = abs(arb(-2) + th_arb**2)  # ~ 0.0800838 < 1
    u_log = float((-eps_re.log()).mid())
    th_real = float(th_arb.mid())

    def reduce_cubic_orbit(triples: list[tuple[int, int, int]]) -> tuple[int, ...]:
        idx = next(i for i, t in enumerate(triples) if any(t))
        a, b, c = triples[idx]
        nrm = abs(norm_ok((a, b, c)))
        v_re = abs(a + b * th_real + c * th_real**2)
        shift = math.log(v_re) - math.log(nrm) / 3.0
        k_pow = int(math.floor(shift / u_log))
        cur = list(triples)
        if k_pow > 0:
            for _ in range(k_pow):
                cur = [mul_ok(t, eps) for t in cur]
        elif k_pow < 0:
            for _ in range(-k_pow):
                cur = [mul_ok(t, eps_inv) for t in cur]
        # Put representative in half-open log-domain (-u_log/2, +u_log/2], where (1,0,0) has log-shift 0 (interior!)
        while True:
            ra, rb, rc = cur[idx]
            ratio_log = float((abs(arb(ra) + arb(rb) * th_arb + arb(rc) * th_arb**2).log() - arb(nrm).log() / 3).mid())
            if ratio_log <= -0.5 * u_log:
                cur = [mul_ok(t, eps_inv) for t in cur]
            elif ratio_log > 0.5 * u_log:
                cur = [mul_ok(t, eps) for t in cur]
            else:
                break
        a, b, c = cur[idx]
        if arb(a) + arb(b) * th_arb + arb(c) * th_arb**2 < 0:
            cur = [(-x, -y, -z) for (x, y, z) in cur]
        return tuple(v for t in cur for v in t)

    # Scan [-E, E]^9! Wait: 5^9 = 1,953,125 -> fast if pruned, or note that for |F| <= 64 and N >= 10^6,
    # U^2 + c1*U*V + V^2 with c1 ~ N/2 can ONLY have tiny norm if U == 0 or V == 0 (or U, V tiny)!
    # Let's do an exact pruned scan over alpha_1, alpha_3 (which fix U = alpha_1 - 2*alpha_3 and U' = alpha_1)
    # Wait: E=1 has 3^9 = 19,683 points; E=2 on w_1, w_2, w_3: let's check E=2 with exact algebraic evaluation!
    # Even faster: precompute all 5^3 = 125 triples in [-2, 2]^3!
    rng3 = list(itertools.product(range(-E, E + 1), repeat=3))
    hits_up_to_pm1 = 0
    orbits = {}
    for al1 in rng3:
        for al3 in rng3:
            U = add_ok(al1, al3, -2)
            Up = al1
            U2 = mul_ok(U, U)
            Up2 = mul_ok(Up, Up)
            for al2 in rng3:
                # Enforce first nonzero in (al1, al2, al3) > 0
                flat = al1 + al2 + al3
                nz = next((x for x in flat if x != 0), 0)
                if nz <= 0:
                    continue
                V = add_ok(al2, al3, 1)
                Vp = al2
                # N_{K/k} = U^2 + c1*U*V + V^2
                UV = mul_ok(U, V)
                V2 = mul_ok(V, V)
                nK_rel = (
                    U2[0] + c1 * UV[0] + V2[0],
                    U2[1] + c1 * UV[1] + V2[1],
                    U2[2] + c1 * UV[2] + V2[2],
                )
                W1 = norm_ok(nK_rel)
                if W1 == 0 or abs(W1) > Y_cap:
                    continue
                UpVp = mul_ok(Up, Vp)
                Vp2 = mul_ok(Vp, Vp)
                nKp_rel = (
                    Up2[0] + c2 * UpVp[0] + 2 * Vp2[0],
                    Up2[1] + c2 * UpVp[1] + 2 * Vp2[1],
                    Up2[2] + c2 * UpVp[2] + 2 * Vp2[2],
                )
                W2 = norm_ok(nKp_rel)
                if W2 == 0:
                    continue
                F_val = abs(W1 * W2)
                if F_val <= Y_cap:
                    hits_up_to_pm1 += 1
                    okey = reduce_cubic_orbit([al1, al2, al3])
                    orbits.setdefault(okey, []).append((F_val, flat))

    # Compute k-rank of all surviving orbits via Q-span under multiplication by 1, theta, theta^2
    rows = []
    for okey in orbits:
        t1, t2, t3 = okey[:3], okey[3:6], okey[6:9]
        for shift_u in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            r_vec = list(mul_ok(t1, shift_u) + mul_ok(t2, shift_u) + mul_ok(t3, shift_u))
            rows.append(r_vec)
    _, rnk_q = fmpq_mat(rows).rref() if rows else (None, 0)
    assert rnk_q % 3 == 0
    full_3coord_orbits = sum(
        1 for okey in orbits if any(okey[:3]) and any(okey[3:6]) and any(okey[6:9])
    )
    return {
        "Y_cap": Y_cap,
        "E": E,
        "hits_up_to_pm1": hits_up_to_pm1,
        "distinct_Ok_unit_orbits": len(orbits),
        "full_3coord_orbits": full_3coord_orbits,
        "orbit_sizes_in_box": {str(k): len(v) for k, v in orbits.items()},
        "orbit_norms": {str(k): v[0][0] for k, v in orbits.items()},
        "spanned_k_rank": rnk_q // 3,
    }


def main() -> None:
    results = {"P1_instances": [], "P2_boxes": [], "P3_60bit": [], "P3_cubic": []}

    # Cell A: k = Q(sqrt(2)), (n, m) in {(3,3), (4,3), (4,4)}
    for n, m, bits, E_box in ((3, 3, 24, 2), (4, 3, 26, 1), (4, 4, 28, 1)):
        P_req = predicted_P(m, 2 * n)
        for seed in range(1, 4):
            N, p, q, A = make_qsqrt2_semiprime(bits, seed)
            b, P, Q, gP, gQ, j_shift = build_divdiff_pair(N, A, n, 2)
            cert = certify_quad_base_deltas(P, Q, gP[:m], gQ[:m], 2)
            Q2 = cert["Q2"]
            exp_mult = math.factorial(m) // 2
            vp_val = exp_mult * vp(Q2, p)
            vq_val = exp_mult * vp(Q2, q)
            p1_ok = cert["abs_Delta_F_zero"] and Q2 != 0 and vp_val >= P_req and vq_val >= P_req
            rec = {
                "cell": "A_Qsqrt2",
                "n": n,
                "m": m,
                "N": N,
                "p": p,
                "q": q,
                "abs_Delta_F_zero": cert["abs_Delta_F_zero"],
                "Q2_nonzero": Q2 != 0,
                "vp_N_Delta_rel": vp_val,
                "vq_N_Delta_rel": vq_val,
                "P_m_2n": P_req,
                "P1_pass": p1_ok,
            }
            results["P1_instances"].append(rec)

            box_res = analyze_qsqrt2_box(N, P, Q, gP[:m], gQ[:m], E_box, cert["H"])
            box_res.update({"n": n, "m": m, "E": E_box, "N": N})
            results["P2_boxes"].append(box_res)

    # Cell B: k = Q(i), (n, m) = (3, 3)
    P_req_B = predicted_P(3, 6)
    for seed in range(1, 4):
        N, p, q, A = make_qi_semiprime(24, seed)
        b, P, Q, gP, gQ, j_shift = build_divdiff_pair(N, A, 3, -1)
        cert = certify_quad_base_deltas(P, Q, gP[:3], gQ[:3], -1)
        Q2 = cert["Q2"]
        exp_mult = math.factorial(3) // 2
        vp_val = exp_mult * vp(Q2, p)
        vq_val = exp_mult * vp(Q2, q)
        p1_ok = cert["abs_Delta_F_zero"] and Q2 != 0 and vp_val >= P_req_B and vq_val >= P_req_B
        results["P1_instances"].append(
            {
                "cell": "B_Qi",
                "n": 3,
                "m": 3,
                "N": N,
                "p": p,
                "q": q,
                "abs_Delta_F_zero": cert["abs_Delta_F_zero"],
                "Q2_nonzero": Q2 != 0,
                "vp_N_Delta_rel": vp_val,
                "vq_N_Delta_rel": vq_val,
                "P_m_2n": P_req_B,
                "P1_pass": p1_ok,
            }
        )

    # Cell C: k = Q(cbrt(3)), (n, m) = (2, 3) on Consultant C's norms-1-and-8 pair
    P_req_C = predicted_P(3, 4)
    for seed in range(1, 4):
        N, p, q, B_root = make_qcbrt3_semiprime(24, seed)
        cert = certify_cubic_c_example(N)
        Q2 = cert["Q2"]
        exp_mult = math.factorial(3) // 2
        vp_val = exp_mult * vp(Q2, p)
        vq_val = exp_mult * vp(Q2, q)
        p1_ok = cert["abs_Delta_F_zero"] and Q2 != 0 and vp_val >= P_req_C and vq_val >= P_req_C
        results["P1_instances"].append(
            {
                "cell": "C_Qcbrt3",
                "n": 2,
                "m": 3,
                "N": N,
                "p": p,
                "q": q,
                "abs_Delta_F_zero": cert["abs_Delta_F_zero"],
                "Q2_nonzero": Q2 != 0,
                "vp_N_Delta_rel": vp_val,
                "vq_N_Delta_rel": vq_val,
                "P_m_2n": P_req_C,
                "P1_pass": p1_ok,
            }
        )
        c_box = analyze_cubic_c_box(N, E=2, Y_cap=64)
        c_box["N"] = N
        results["P3_cubic"].append(c_box)

    # 60-bit Q(sqrt(2)) (n, m) = (4, 3), E = 1 for P2 + P3
    for seed in range(1, 4):
        N, p, q, A = make_qsqrt2_semiprime(60, seed)
        b, P, Q, gP, gQ, j_shift = build_divdiff_pair(N, A, 4, 2)
        cert = certify_quad_base_deltas(P, Q, gP[:3], gQ[:3], 2)
        box60 = analyze_qsqrt2_box(N, P, Q, gP[:3], gQ[:3], 1, cert["H"])
        box60.update({"n": 4, "m": 3, "E": 1, "N": N, "bits": 60})
        results["P2_boxes"].append(box60)
        results["P3_60bit"].append(box60)

    # Summary scores
    p1_passes = sum(1 for r in results["P1_instances"] if r["P1_pass"])
    p2_passes = sum(1 for r in results["P2_boxes"] if r["R_unit_bound_holds"])
    e2_dedup = [
        (r["z_primitive_up_to_pm1"], r["distinct_Ok_unit_orbits"])
        for r in results["P2_boxes"]
        if r["E"] == 2
    ]
    results["summary"] = {
        "P1_pass_count": f"{p1_passes}/{len(results['P1_instances'])}",
        "P2_R_unit_bound_pass_count": f"{p2_passes}/{len(results['P2_boxes'])}",
        "P2_E2_zprim_vs_Ok_orbits": e2_dedup,
        "P3_60bit_spanned_k_dims": [r["spanned_k_dim_of_all_survivors"] for r in results["P3_60bit"]],
        "P3_60bit_last_gen_orbits": [r["orbits_with_last_generator_nonzero"] for r in results["P3_60bit"]],
        "P3_60bit_sqrt_Y0_vs_3r": [(r["B3_sqrt_Y0_floor"], 3**6) for r in results["P3_60bit"]],
        "P3_cubic_hits_vs_orbits_vs_krank": [
            (r["hits_up_to_pm1"], r["distinct_Ok_unit_orbits"], r["spanned_k_rank"])
            for r in results["P3_cubic"]
        ],
    }

    out_path = HERE / "results.json"
    out_path.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results["summary"], indent=2))


if __name__ == "__main__":
    main()
