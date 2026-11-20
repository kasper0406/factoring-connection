"""Pre-registered adversary experiments for Chase Lemma 2 (ch_L2g).
Strictly read-only with respect to /Volumes/git/factoring-connection; writes results.json in ch_L2g."""
from __future__ import annotations

import itertools
import json
import math
import random
import sys
from functools import reduce
from pathlib import Path

ROOT = Path("/Volumes/git/factoring-connection")
HERE = Path("/Volumes/git/factoring-connection/codex/scratch/ch_L2g")
sys.path[:0] = [
    str(ROOT / "experiments/masullo_ext/mext_c"),
    str(ROOT / "experiments/subthird/codex"),
    str(ROOT),
]

import cover as cov
import experiment as ex
from flint import acb, arb, ctx, fmpq_mat, fmpz_mat, fmpz_poly as Poly
from sympy import cyclotomic_poly, factorint, integer_nthroot, isprime, nextprime, primerange, sqrt_mod, symbols


def eval_row_acb(row: list[acb], z: tuple[int, ...]) -> acb:
    s = acb(0)
    for r_val, c_val in zip(row, z):
        if c_val:
            s += r_val * int(c_val)
    return s


def cross_norm_lower(r1: list[acb], r2: list[acb]) -> arb:
    """Lower bound on ||L_1 x L_2||_2 where L_a = (r_a[0], r_a[2], r_a[4]) in C^3."""
    a0, a1, a2 = r1[0], r1[2], r1[4]
    b0, b1, b2 = r2[0], r2[2], r2[4]
    c0 = a1 * b2 - a2 * b1
    c1 = a2 * b0 - a0 * b2
    c2 = a0 * b1 - a1 * b0
    sq = abs(c0) ** 2 + abs(c1) ** 2 + abs(c2) ** 2
    return sq.sqrt()


def prefix_pair(row1: list[acb], row2: list[acb], points: list[tuple[int, ...]]) -> list[tuple[int, ...]]:
    """Parameter-free codimension-2 analogue of cover.py:prefix(), sorting by |L_i(z) * L_j(z)|."""
    m1 = [[v * j for j in range(-2, 3)] for v in row1]
    m2 = [[v * j for j in range(-2, 3)] for v in row2]
    ordered = []
    for z in points:
        v1 = sum((tab[j + 2] for tab, j in zip(m1, z)), acb(0))
        v2 = sum((tab[j + 2] for tab, j in zip(m2, z)), acb(0))
        prod_val = abs(v1 * v2)
        ordered.append((prod_val.lower(), prod_val.upper(), z))
    ordered.sort(key=lambda item: item[0])
    provisional: list[tuple[int, ...]] = []
    for lower, upper, z in ordered:
        cov.add_independent(provisional, z)
        if len(provisional) == 3:
            threshold = lower
            break
    else:
        raise ValueError("Coefficient box has deficient rank")
    basis: list[tuple[int, ...]] = []
    for lower, upper, z in ordered:
        if lower >= threshold:
            break
        if upper < threshold:
            cov.add_independent(basis, z)
    assert len(basis) < 3
    return basis


def prefix_pair_balanced(row1: list[acb], row2: list[acb], points: list[tuple[int, ...]]) -> list[tuple[int, ...]]:
    """Balanced L-infty tube max(|L_i(z)|/||L_i||_2, |L_j(z)|/||L_j||_2) around L_i = L_j = 0."""
    n1 = (abs(row1[0]) ** 2 + abs(row1[2]) ** 2 + abs(row1[4]) ** 2).sqrt()
    n2 = (abs(row2[0]) ** 2 + abs(row2[2]) ** 2 + abs(row2[4]) ** 2).sqrt()
    m1 = [[v * j for j in range(-2, 3)] for v in row1]
    m2 = [[v * j for j in range(-2, 3)] for v in row2]
    ordered = []
    for z in points:
        v1 = abs(sum((tab[j + 2] for tab, j in zip(m1, z)), acb(0))) / n1
        v2 = abs(sum((tab[j + 2] for tab, j in zip(m2, z)), acb(0))) / n2
        mx = v1.max(v2)
        ordered.append((mx.lower(), mx.upper(), z))
    ordered.sort(key=lambda item: item[0])
    provisional: list[tuple[int, ...]] = []
    for lower, upper, z in ordered:
        cov.add_independent(provisional, z)
        if len(provisional) == 3:
            threshold = lower
            break
    else:
        raise ValueError("Coefficient box has deficient rank")
    basis: list[tuple[int, ...]] = []
    for lower, upper, z in ordered:
        if lower >= threshold:
            break
        if upper < threshold:
            cov.add_independent(basis, z)
    assert len(basis) < 3
    return basis


# ==============================================================================
# T1: 72-bit Escaping Witness Anatomy & Rank-3 Two-Hyperplane Liouville Lemma
# ==============================================================================
def run_T1() -> dict:
    ctx.prec = 256
    units = cov.pair_units()
    N72 = ex.fixture(72, 19)
    _, b72, P72, Q72 = ex.select(N72)
    f72, g72 = P72 * P72 - 2, Q72 * Q72 - 2
    U72, V72 = [ex.expand(X, X * X - 2, ex.generators(X, "divided", 3)) for X in (P72, Q72)]
    residues = [int(u(b72)) % N72 for u in U72]
    points = [
        z
        for z in itertools.product(range(-2, 3), repeat=6)
        if next((j for j in z if j), 0) > 0
        and math.gcd(sum(t * c for t, c in zip(z, residues)), N72) == 1
    ]
    rows, _ = ex.rows_at(ex.coeffs(f72), list(map(ex.coeffs, U72)), ex.coeffs(g72), list(map(ex.coeffs, V72)), 256)

    # Verify relative block structure: row[2*j+1] == row[1]*row[2*j] with row[1] == +-sqrt(2)
    pos_idx, neg_idx = [], []
    for idx, r in enumerate(rows):
        for j in range(3):
            assert (r[2 * j + 1] - r[1] * r[2 * j]).contains(0)
        if r[1].real > 0:
            pos_idx.append(idx)
        else:
            neg_idx.append(idx)
    assert len(pos_idx) == 16 and len(neg_idx) == 16

    # 1. Inspect single-hyperplane spaces W_i and thresholds t_i
    spaces_1h = []
    t_1h_lowers = []
    plane_a3_zero_count = 0
    for r in rows:
        m_tab = [[v * j for j in range(-2, 3)] for v in r]
        ordered = []
        for z in points:
            val = abs(sum((tab[j + 2] for tab, j in zip(m_tab, z)), acb(0)))
            ordered.append((val.lower(), val.upper(), z))
        ordered.sort(key=lambda x: x[0])
        prov = []
        for lo, up, z in ordered:
            cov.add_independent(prov, z)
            if len(prov) == 3:
                t_1h_lowers.append(float(lo))
                thr = lo
                break
        bas = []
        for lo, up, z in ordered:
            if lo >= thr:
                break
            if up < thr:
                cov.add_independent(bas, z)
        spaces_1h.append(bas)
        if len(bas) == 2 and all(v[4] == 0 and v[5] == 0 for v in bas):
            plane_a3_zero_count += 1

    log2_prod_t1h = sum(math.log2(x) for x in t_1h_lowers)

    # Witness z* = (-2, 1, -1, 1, -1, 0)
    z_star = (-2, 1, -1, 1, -1, 0)
    vals_side1 = sorted(float(abs(eval_row_acb(rows[i], z_star)).mid()) for i in range(16))
    vals_side2 = sorted(float(abs(eval_row_acb(rows[i], z_star)).mid()) for i in range(16, 32))

    # 2. Analytic Rank-3 Two-Hyperplane Relative Liouville Lemma check across all 240 same-sigma pairs & dyadic scales
    # Constant for e=2, m=3, E=2: denom = 27*sqrt(3) * ((1 + sqrt(2)) * 2)^4
    ctx.prec = 256
    denom_const = arb(27) * arb(3).sqrt() * ((arb(1) + arb(2).sqrt()) * arb(2)) ** 4
    all_box_pts = [z for z in itertools.product(range(-2, 3), repeat=6) if any(z) and next(j for j in z if j) > 0]
    row_evals = [[abs(eval_row_acb(r, z)).upper() for z in all_box_pts] for r in rows]

    dyadic_scales = [arb(2) ** k for k in range(-12, 13)]
    analytic_checks = 0
    analytic_violations = 0
    same_sigma_pairs = list(itertools.combinations(pos_idx, 2)) + list(itertools.combinations(neg_idx, 2))
    assert len(same_sigma_pairs) == 240

    for i, j in same_sigma_pairs:
        cn = cross_norm_lower(rows[i], rows[j])
        two_tau_star = (arb(2) * cn / denom_const).lower()
        ev_i = row_evals[i]
        ev_j = row_evals[j]
        for s_scale in dyadic_scales:
            s_lo = s_scale.lower()
            u_lo = (arb(two_tau_star) / s_scale).lower()
            # Filter points strictly inside |L_i(z)| < s_lo and |L_j(z)| < u_lo (using upper() <= bound)
            # Wait: if upper() < bound, |L(z)| < bound is 100% certified; to test the theorem AGAINST violations,
            # any true violation |L_i| < s and |L_j| < 2*tau^*/s at high precision has mid/upper < bound too!
            # Use certified mid/lower check: if lower() < bound could include points, check exact lower() < bound!
            cand_vecs = []
            for p_idx, z in enumerate(all_box_pts):
                # Use lower bound of |L_i(z)| so we OVER-approximate the set (stronger test: even borderline points included!)
                # Wait: compute exact ball if close; here 256-bit precision has radius < 1e-70, so lower vs upper never splits
                if ev_i[p_idx] < s_lo and ev_j[p_idx] < u_lo:
                    cov.add_independent(cand_vecs, z)
                    if len(cand_vecs) == 3:
                        analytic_violations += 1
                        break
            analytic_checks += 1

    # 3. Parameter-free Two-Hyperplane Prefix Cover on 72-bit and 96-bit divided & power fixtures
    spaces_2h_prod = [prefix_pair(rows[i], rows[j], points) for i, j in same_sigma_pairs]
    spaces_2h_bal = [prefix_pair_balanced(rows[i], rows[j], points) for i, j in same_sigma_pairs]

    distinct_2h_planes = set()
    for bas in spaces_2h_prod:
        if len(bas) == 2:
            # Canonical normal in P^2(Q(sqrt(2)))? Just count how many equal {alpha_3 = 0}
            pass
    non_a3_zero_2h = sum(1 for bas in spaces_2h_prod if not all(v[4] == 0 and v[5] == 0 for v in bas))

    z_star_cov_1h = any(ex.kranks(bas + [z_star]) == len(bas) for bas in spaces_1h)
    z_star_cov_2h_prod = any(ex.kranks(bas + [z_star]) == len(bas) for bas in spaces_2h_prod)
    z_star_cov_2h_bal = any(ex.kranks(bas + [z_star]) == len(bas) for bas in spaces_2h_bal)
    z_star_cov_union = z_star_cov_1h or z_star_cov_2h_prod or z_star_cov_2h_bal

    # Full smooth parity rank comparison on 72-bit divided
    cap72 = 3**32 * 2 ** ((72 + 15) // 16)
    B72 = 8192
    primes72 = list(primerange(B72 + 1))
    primorial72 = math.prod(primes72)
    fields72 = [ex.Ideals(X, primes72) for X in (f72, g72)]

    keys_seen = set()
    uncov_1h, uncov_2h_prod, uncov_union = 0, 0, 0
    rows_1h, rows_2h_prod, rows_union, rows_all = [], [], [], []
    columns: dict = {}
    for z in points:
        a_el, w_el = ex.element(U72, z), ex.element(V72, z)
        v1, v2 = abs(int(f72.resultant(a_el))), abs(int(g72.resultant(w_el)))
        val = v1 * v2
        if not 0 < val <= cap72:
            continue
        key = cov.unit_key(z, units)
        if key in keys_seen:
            continue
        keys_seen.add(key)
        c1 = any(ex.kranks(bas + [z]) == len(bas) for bas in spaces_1h)
        c2 = any(ex.kranks(bas + [z]) == len(bas) for bas in spaces_2h_prod)
        c2b = any(ex.kranks(bas + [z]) == len(bas) for bas in spaces_2h_bal)
        cu = c1 or c2 or c2b
        if not c1:
            uncov_1h += 1
        if not c2:
            uncov_2h_prod += 1
        if not cu:
            uncov_union += 1
        if ex.smooth(val, primorial72):
            mask = 0
            for side, (alpha, nrm) in enumerate(((a_el, v1), (w_el, v2))):
                for pr in fields72[side].parity(alpha, nrm, primes72):
                    col = columns.setdefault((side, pr), len(columns))
                    mask |= 1 << col
            rows_all.append(mask)
            if c1:
                rows_1h.append(mask)
            if c2:
                rows_2h_prod.append(mask)
            if cu:
                rows_union.append(mask)

    # Also check scaling of min full-support |F| vs cap across bits = 72, 96, 120 on divided-difference n=8, m=3
    scaling_div = []
    for bits_test in (72, 96, 120):
        Nt = ex.fixture(bits_test, 19)
        _, bt, Pt, Qt = ex.select(Nt)
        ft, gt = Pt * Pt - 2, Qt * Qt - 2
        Ut, Vt = [ex.expand(X, X * X - 2, ex.generators(X, "divided", 3)) for X in (Pt, Qt)]
        cap_t = 3**32 * 2 ** ((bits_test + 15) // 16)
        min_full_F = None
        full_under_cap = 0
        for z in itertools.product(range(-2, 3), repeat=6):
            if next((j for j in z if j), 0) <= 0:
                continue
            if not all(z[2 * j] != 0 or z[2 * j + 1] != 0 for j in range(3)):
                continue
            a_el, w_el = ex.element(Ut, z), ex.element(Vt, z)
            v1, v2 = abs(int(ft.resultant(a_el))), abs(int(gt.resultant(w_el)))
            prod_v = v1 * v2
            if prod_v > 0:
                if min_full_F is None or prod_v < min_full_F:
                    min_full_F = prod_v
                if prod_v <= cap_t:
                    full_under_cap += 1
        scaling_div.append(
            {
                "bits": bits_test,
                "b": bt,
                "log2_cap": round(math.log2(cap_t), 2),
                "log2_min_full_F": round(math.log2(min_full_F), 2),
                "full_support_vectors_under_cap": full_under_cap,
            }
        )

    return {
        "N72": str(N72),
        "b72": b72,
        "plane_a3_zero_in_32_1h_spaces": plane_a3_zero_count,
        "log2_prod_t1h": round(log2_prod_t1h, 2),
        "log2_witness_F": round(math.log2(108241 * 3175924591), 2),
        "log2_cap72": round(math.log2(cap72), 2),
        "witness_side1_sorted_abs": [round(x, 4) for x in vals_side1],
        "witness_side2_sorted_abs": [round(x, 4) for x in vals_side2],
        "analytic_2h_checks": analytic_checks,
        "analytic_2h_violations": analytic_violations,
        "non_a3_zero_in_240_2h_prod_spaces": non_a3_zero_2h,
        "witness_covered_1h": z_star_cov_1h,
        "witness_covered_2h_prod": z_star_cov_2h_prod,
        "witness_covered_2h_balanced": z_star_cov_2h_bal,
        "capped_unit_classes": len(keys_seen),
        "uncovered_1h": uncov_1h,
        "uncovered_2h_prod": uncov_2h_prod,
        "uncovered_union": uncov_union,
        "smooth_rank_all": ex.rank(rows_all),
        "smooth_rank_1h": ex.rank(rows_1h),
        "smooth_rank_2h_prod": ex.rank(rows_2h_prod),
        "smooth_rank_union": ex.rank(rows_union),
        "incremental_rank_after_2h_prod": ex.rank(rows_all) - ex.rank(rows_2h_prod),
        "incremental_rank_after_union": ex.rank(rows_all) - ex.rank(rows_union),
        "scaling_divided_72_96_120": scaling_div,
    }


def make_qsqrt2_semiprime(bits: int, seed: int) -> tuple[int, int, int, int]:
    rng = random.Random(91001 * bits + seed)
    lo = 1 << (bits // 2 - 1)
    hi = (1 << (bits // 2)) - 1
    while True:
        p = int(nextprime(rng.randrange(lo, hi)))
        if p % 8 not in (1, 7) or p % 15 == 0:
            continue
        q = int(nextprime(rng.randrange(lo, hi)))
        if q == p or q % 8 not in (1, 7) or q % 15 == 0:
            continue
        N = p * q
        if N.bit_length() != bits:
            continue
        rp = int(sqrt_mod(2, p))
        rq = int(sqrt_mod(2, q))
        inv_p = pow(p, -1, q)
        A = (rp + p * ((rq - rp) * inv_p % q)) % N
        if A < 10:
            continue
        assert (A * A - 2) % N == 0
        return N, p, q, A


def qsqrt2_orbit_rep(z: tuple[int, ...]) -> tuple[int, ...]:
    """Canonical O_k^*-orbit representative for k = Q(sqrt(2))."""
    m = len(z) // 2
    pairs = [(z[2 * j], z[2 * j + 1]) for j in range(m)]
    idx = next(j for j, (a, b) in enumerate(pairs) if a != 0 or b != 0)
    ctx.prec = 256
    s2 = arb(2).sqrt()
    u_log = (arb(1) + s2).log()
    u_f = float(u_log.mid())
    a0, b0 = pairs[idx]
    norm0 = abs(a0 * a0 - 2 * b0 * b0)
    cur = abs(arb(a0) + arb(b0) * s2).log() - arb(norm0).log() / 2
    t = -math.floor(float((cur / u_log).mid()))

    def step_u(ps: list[tuple[int, int]], k_step: int) -> list[tuple[int, int]]:
        out = list(ps)
        if k_step > 0:
            for _ in range(k_step):
                out = [(a + 2 * b, a + b) for (a, b) in out]
        elif k_step < 0:
            for _ in range(-k_step):
                out = [(-a + 2 * b, a - b) for (a, b) in out]
        return out

    red = step_u(pairs, t)
    while True:
        ra, rb = red[idx]
        r_log = float((abs(arb(ra) + arb(rb) * s2).log() - arb(norm0).log() / 2).mid())
        if r_log <= -0.5 * u_f:
            red = step_u(red, 1)
        elif r_log > 0.5 * u_f:
            red = step_u(red, -1)
        else:
            break
    ra, rb = red[idx]
    if arb(ra) + arb(rb) * s2 < 0:
        red = [(-a, -b) for (a, b) in red]
    return tuple(x for ab in red for x in ab)


# ==============================================================================
# T2: Multi-Cone Lagrange Relative-Unit / Small-Norm Seeds (2^m - 2 Orbits)
# ==============================================================================
def run_T2() -> list[dict]:
    T_var = Poly([0, 1])
    R0 = Poly([3, 2])  # 2T + 3: at r in (0, 1, -1, 2), R0(r)^2 - 2 in {7, 23, -1, 47}
    out = []
    for n, m, bits in ((6, 3, 48), (8, 4, 48)):
        N, p, q, A0 = make_qsqrt2_semiprime(bits, 7)
        nodes = (0, 1, -1, 2)[:m]
        Qm = reduce(lambda a, b: a * b, (T_var - Poly([r_val]) for r_val in nodes))
        deg_A = n - m
        b_base = max(6, int(integer_nthroot(N, deg_A)[0]))
        found = None
        for b in range(b_base, b_base + 200):
            qb = int(Qm(b))
            if qb == 0 or math.gcd(abs(qb), N) != 1:
                continue
            # Solve A0 - R0(b) + t*N == 0 (mod qb)
            rem0 = (A0 - int(R0(b))) % abs(qb)
            t_shift = ((-rem0) * pow(N % abs(qb), -1, abs(qb))) % abs(qb)
            target = (A0 - int(R0(b)) + t_shift * N) // qb
            # Expand target as monic degree deg_A in base b
            rem = target - b**deg_A
            digs = []
            for _ in range(deg_A - 1):
                d_i = (rem + b // 2) % b - b // 2
                digs.append(d_i)
                rem = (rem - d_i) // b
            A_poly = Poly(digs + [rem, 1])
            assert int(R0(b) + Qm(b) * A_poly(b)) % N == A0 % N
            P = R0 + Qm * A_poly
            f = P * P - 2
            if not ex.irreducible(f):
                continue
            for j_sh in (1, -1, 2, -2, 3, -3, 5, -5):
                B_poly = A_poly + j_sh * (T_var - Poly([b]))
                Q = R0 + Qm * B_poly
                g = Q * Q - 2
                if Q.degree() == n and int(Q[n]) == 1 and ex.irreducible(g):
                    found = (b, P, Q, f, g)
                    break
            if found:
                break
        assert found is not None
        b, P, Q, f, g = found
        # Lagrange basis w_j = Q_m(T) // (T - r_j) on both sides
        lag_gens = [Qm // (T_var - Poly([r_val])) for r_val in nodes]
        U = ex.expand(P, f, lag_gens)
        V = ex.expand(Q, g, lag_gens)
        residues = [int(u(b)) % N for u in U]
        assert residues == [int(v(b)) % N for v in V]

        # Scan [-1, 1]^{2m}
        small_cap = 60**m
        small_orbits: dict[tuple[int, ...], tuple[int, tuple[int, ...], Poly, Poly, int, int]] = {}
        min_full_F = None
        all_primes = set()
        for z in itertools.product((-1, 0, 1), repeat=2 * m):
            if next((x for x in z if x), 0) <= 0:
                continue
            if math.gcd(sum(t * c for t, c in zip(z, residues)), N) != 1:
                continue
            a_el, w_el = ex.element(U, z), ex.element(V, z)
            v1, v2 = abs(int(f.resultant(a_el))), abs(int(g.resultant(w_el)))
            F_val = v1 * v2
            if F_val == 0:
                continue
            is_full = all(z[2 * j] != 0 or z[2 * j + 1] != 0 for j in range(m))
            if is_full:
                if min_full_F is None or F_val < min_full_F:
                    min_full_F = F_val
            if F_val <= small_cap:
                okey = qsqrt2_orbit_rep(z)
                if okey not in small_orbits or F_val < small_orbits[okey][0]:
                    small_orbits[okey] = (F_val, z, a_el, w_el, v1, v2)

        for _, _, _, _, v1, v2 in small_orbits.values():
            all_primes.update(factorint(v1).keys())
            all_primes.update(factorint(v2).keys())
        sorted_prs = sorted(all_primes)
        ideals = [ex.Ideals(X, sorted_prs) for X in (f, g)]
        cols: dict = {}
        masks = []
        orbit_vecs = []
        for F_val, z, a_el, w_el, v1, v2 in small_orbits.values():
            orbit_vecs.append(z)
            msk = 0
            for side, (al, nrm) in enumerate(((a_el, v1), (w_el, v2))):
                for pr in ideals[side].parity(al, nrm, sorted_prs):
                    msk |= 1 << cols.setdefault((side, pr), len(cols))
            masks.append(msk)

        out.append(
            {
                "n": n,
                "m": m,
                "bits": bits,
                "b": b,
                "small_cap": small_cap,
                "distinct_small_Ok_orbits": len(small_orbits),
                "predicted_2m_minus_2": 2**m - 2,
                "spanned_k_rank_of_small_orbits": ex.kranks(orbit_vecs),
                "ideal_parity_rank_of_all_small_orbits": ex.rank(masks),
                "full_support_small_orbits": sum(
                    1 for _, z, *_ in small_orbits.values() if all(z[2 * j] or z[2 * j + 1] for j in range(m))
                ),
                "log_min_full_F_over_log_N": round(math.log(min_full_F) / math.log(N), 4),
                "predicted_4_m_minus_1_log_b_over_log_N": round(4 * (m - 1) * math.log(b) / math.log(N), 4),
            }
        )
    return out


# ==============================================================================
# T3: Common-Denominator Zero-Remainder Simultaneous Multi-Embedding Cancellation
# ==============================================================================
def run_T3() -> list[dict]:
    ctx.prec = 256
    T_var = Poly([0, 1])
    D_poly = Poly([1, 0, 1])  # T^2 + 1
    out = []
    for n, m, bits in ((6, 3, 48), (8, 4, 60)):
        N, p, q, A0 = make_qsqrt2_semiprime(bits, 11)
        deg_A = n - 2
        b_base = max(6, int(integer_nthroot(N, deg_A)[0]))
        found = None
        for b in range(b_base, b_base + 200):
            db = b * b + 1
            if math.gcd(db, N) != 1:
                continue
            t_shift = ((-A0) * pow(N % db, -1, db)) % db
            target = (A0 + t_shift * N) // db
            rem = target - b**deg_A
            digs = []
            for _ in range(deg_A - 1):
                d_i = (rem + b // 2) % b - b // 2
                digs.append(d_i)
                rem = (rem - d_i) // b
            A_poly = Poly(digs + [rem, 1])
            P = D_poly * A_poly
            f = P * P - 2
            if not ex.irreducible(f):
                continue
            for j_sh in (1, -1, 2, -2, 3, -3):
                B_poly = A_poly + j_sh * (T_var - Poly([b]))
                Q = D_poly * B_poly
                g = Q * Q - 2
                if Q.degree() == n and int(Q[n]) == 1 and ex.irreducible(g):
                    found = (b, P, Q, A_poly, B_poly, f, g)
                    break
            if found:
                break
        assert found is not None
        b, P, Q, A_poly, B_poly, f, g = found
        gens_P = [Poly([1])] + [(T_var ** (j - 1)) * A_poly for j in range(1, m)]
        gens_Q = [Poly([1])] + [(T_var ** (j - 1)) * B_poly for j in range(1, m)]
        U = ex.expand(P, f, gens_P)
        V = ex.expand(Q, g, gens_Q)
        residues = [int(u(b)) % N for u in U]
        assert residues == [int(v(b)) % N for v in V]

        rows, _ = ex.rows_at(ex.coeffs(f), list(map(ex.coeffs, U)), ex.coeffs(g), list(map(ex.coeffs, V)), 256)
        D_abs = 4 * n
        X0 = (3**D_abs) * (2 ** ((bits + 2 * n - 1) // (2 * n)))

        max_small_embs = 0
        min_full_F = None
        full_under_X0 = 0
        for z in itertools.product((-1, 0, 1), repeat=2 * m):
            if next((x for x in z if x), 0) <= 0:
                continue
            if not all(z[2 * j] != 0 or z[2 * j + 1] != 0 for j in range(m)):
                continue
            if math.gcd(sum(t * c for t, c in zip(z, residues)), N) != 1:
                continue
            a_el, w_el = ex.element(U, z), ex.element(V, z)
            v1, v2 = abs(int(f.resultant(a_el))), abs(int(g.resultant(w_el)))
            F_val = v1 * v2
            if F_val == 0:
                continue
            if min_full_F is None or F_val < min_full_F:
                min_full_F = F_val
                small_cnt = sum(1 for r in rows if float(abs(eval_row_acb(r, z)).mid()) <= 10.0)
                max_small_embs = max(max_small_embs, small_cnt)
            if F_val <= X0:
                full_under_X0 += 1

        out.append(
            {
                "n": n,
                "m": m,
                "D": D_abs,
                "bits": bits,
                "b": b,
                "small_embeddings_at_best_full_point": max_small_embs,
                "predicted_diffuse_small_embeddings_4_n_minus_2": 4 * (n - 2),
                "log2_X0": round(math.log2(X0), 2),
                "log2_min_full_F": round(math.log2(min_full_F), 2),
                "log_min_full_F_over_log_N": round(math.log(min_full_F) / math.log(N), 4),
                "predicted_4_m_minus_2_log_b_over_log_N": round(4 * (m - 2) * math.log(b) / math.log(N), 4),
                "full_support_under_X0": full_under_X0,
            }
        )
    return out


# ==============================================================================
# T4: Factoring-Free Multiquadratic Relative-Unit Tower with Exact N-Cancellation
# ==============================================================================
def run_T4() -> list[dict]:
    ctx.prec = 512
    out = []
    for t_dim in (2, 3):
        N, p, q, _ = make_qsqrt2_semiprime(48, 20 + t_dim)
        # Ensure gcd(N, 30) == 1
        while math.gcd(N, 30) != 1:
            N, p, q, _ = make_qsqrt2_semiprime(48, random.randint(50, 500))
        b_list = (2, 3, 5)[:t_dim]
        u_list = (2, 7, 23)[:t_dim]
        # Build f_i(T) = T^2 - c_i T + 1, g_i(T) = T^2 - c'_i T + u_i with c_i, c'_i in Z via modular inverse ONLY!
        c_list, cp_list, r_fracs = [], [], []
        s_roots, sp_roots = [], []
        for b_i, u_i in zip(b_list, u_list):
            # Solve b_i * c_i == b_i^2 + 1 (mod N) with c_i = (a_i * N + b_i^2 + 1) // b_i
            a_i = (-(b_i * b_i + 1) * pow(N % b_i, -1, b_i)) % b_i
            if a_i == 0:
                a_i = b_i
            c_i = (a_i * N + b_i * b_i + 1) // b_i
            assert (b_i * b_i - c_i * b_i + 1) % N == 0
            ap_i = (-(b_i * b_i + u_i) * pow(N % b_i, -1, b_i)) % b_i
            if ap_i == 0:
                ap_i = b_i
            cp_i = (ap_i * N + b_i * b_i + u_i) // b_i
            assert (b_i * b_i - cp_i * b_i + u_i) % N == 0
            c_list.append(c_i)
            cp_list.append(cp_i)
            r_fracs.append((a_i, b_i))

            # Exact arb roots of T^2 - c_i T + 1 = 0: (c_i +- sqrt(c_i^2 - 4))/2
            ci_arb = arb(c_i)
            disc_s = (ci_arb * ci_arb - 4).sqrt()
            s_big = (ci_arb + disc_s) / 2
            s_small = arb(1) / s_big
            s_roots.append((s_small, s_big))

            cpi_arb = arb(cp_i)
            disc_sp = (cpi_arb * cpi_arb - 4 * u_i).sqrt()
            sp_big = (cpi_arb + disc_sp) / 2
            sp_small = arb(u_i) / sp_big
            sp_roots.append((sp_small, sp_big))

        # Construct explicit integer N-cancelling full-support vector (alpha_0, alpha_1, ..., alpha_t) in Z^{t+1} subset O_k^{t+1}:
        # r_1 = a_1/b_1, r_2 = a_2/b_2: for t=2, alpha_1 = a_2 * b_1, alpha_2 = -a_1 * b_2 => alpha_1*r_1 + alpha_2*r_2 == 0!
        # For t=3: alpha_1 = b_1 * a_3, alpha_2 = b_2 * a_3, alpha_3 = -2 * a_1 * ... wait:
        # alpha_1 * (a_1/b_1) + alpha_2 * (a_2/b_2) + alpha_3 * (a_3/b_3) == 0:
        # take alpha_1 = b_1 * a_2 * a_3, alpha_2 = b_2 * a_1 * a_3, alpha_3 = -2 * b_3 * a_1 * a_2!
        if t_dim == 2:
            (a1, b1), (a2, b2) = r_fracs
            g_cd = math.gcd(a2 * b1, a1 * b2)
            al_vec = [1, (a2 * b1) // g_cd, -(a1 * b2) // g_cd]
        else:
            (a1, b1), (a2, b2), (a3, b3) = r_fracs
            v1, v2, v3 = b1 * a2 * a3, b2 * a1 * a3, -2 * b3 * a1 * a2
            g_cd = reduce(math.gcd, (abs(v1), abs(v2), abs(v3)))
            al_vec = [1, v1 // g_cd, v2 // g_cd, v3 // g_cd]

        # Evaluate Side-1 and Side-2 norms over k (then square for e=2 since alpha_j in Z subset O_k!)
        n_rel = 1 << t_dim
        emb_vals_s1 = []
        prod_s1 = arb(1)
        prod_s2 = arb(1)
        small_embs_s1 = 0
        large_embs_s1 = 0
        for eps_bits in itertools.product((0, 1), repeat=t_dim):
            val1 = arb(al_vec[0])
            val2 = arb(al_vec[0])
            for idx_i, bit in enumerate(eps_bits):
                val1 += arb(al_vec[idx_i + 1]) * s_roots[idx_i][bit]
                val2 += arb(al_vec[idx_i + 1]) * sp_roots[idx_i][bit]
            av1 = abs(val1)
            av2 = abs(val2)
            prod_s1 *= av1
            prod_s2 *= av2
            f_v1 = float(av1.mid())
            emb_vals_s1.append(round(math.log2(max(f_v1, 1e-30)), 2))
            if f_v1 < 1000.0:
                small_embs_s1 += 1
            else:
                large_embs_s1 += 1

        # Exact integer certification of N_{K/k} and N_{K'/k}!
        norm_s1_rel = int(prod_s1.unique_fmpz())
        norm_s2_rel = int(prod_s2.unique_fmpz())
        # Over k = Q(sqrt(2)) (e = 2), absolute norm is rel_norm^2
        log_F_over_log_N = (2 * math.log(abs(norm_s1_rel)) + 2 * math.log(abs(norm_s2_rel))) / math.log(N)

        # Also check the 2^t monomials s^S: on side 1 norm is 1; on side 2 norm is prod_{i in S} u_i^(2^t) over Q(sqrt(2))?
        # Wait: in K' of degree 2^t over k (degree 2^(t+1) over Q), N_{K'/Q}(s'_i) = (u_i)^(2 * 2^(t-1)) = u_i^(2^t).
        # More importantly, every monomial prod_{i in S} (s_i, s'_i) is a product of the t generators (s_i, s'_i),
        # so over ANY order/field its ideal valuation vector is the sum over i in S of v(s_i, s'_i) -> rank <= t!
        out.append(
            {
                "t": t_dim,
                "n_rel": n_rel,
                "m": t_dim + 1,
                "N": N,
                "r_fracs": r_fracs,
                "cancelling_full_support_alpha": al_vec,
                "side1_log2_embedding_moduli": emb_vals_s1,
                "side1_O1_embeddings": small_embs_s1,
                "side1_ON_embeddings": large_embs_s1,
                "sperner_min_ON_embeddings": n_rel - math.comb(t_dim, t_dim // 2),
                "certified_exact_rel_norm_side1": str(abs(norm_s1_rel)),
                "log_F_over_log_N": round(log_F_over_log_N, 4),
                "monomial_small_orbits": n_rel,
                "max_multiplicative_parity_rank_of_monomials": t_dim,
            }
        )
    return out


# ==============================================================================
# T5: CM / Cyclotomic Zero-Height Side 1 (K = Q(sqrt(2), zeta_15)) + Perturbed Side 2
# ==============================================================================
def run_T5() -> list[dict]:
    ctx.prec = 256
    x_sym = symbols("x")
    phi15_coeffs = [int(c) for c in reversed(cyclotomic_poly(15, x_sym).as_poly(x_sym).all_coeffs())]
    Phi15 = Poly(phi15_coeffs)  # degree 8, irreducible over Q(sqrt(2))
    assert Phi15.degree() == 8
    T_var = Poly([0, 1])

    # Build Side 1 absolute poly: Phi_15(T) (Wait: K = Q(sqrt(2))[T]/(Phi_15(T)) has Z-basis {1, sqrt(2)} x {1, T, ..., T^7})
    # Better: use exact relative resultant over k = Q(sqrt(2)) via certified acb roots of Phi_15(T) = 0 (at both +-sqrt(2))
    # AND Side 2: g_rel(T) = Phi_15(T) + (T - b)*(1 + sqrt(2)), whose Q-norm is:
    # G_abs(T) = (Phi_15(T) + (T - b))^2 - 2*(T - b)^2!
    # Check: if A0^2 == 2 (mod N) and Phi_15(b) == 0 (mod N), then at T = b, Phi_15(b) == 0 (mod N) AND G_abs(b) = Phi_15(b)^2 == 0 (mod N)!
    # Even simpler: pick b such that M_val = int(Phi_15(b)) has primes splitting in Q(sqrt(2)), or take N = int(Phi_15(b)) directly
    # (squarefree SNFS/CM modulus N = Phi_15(b) where 2 is a QR mod every prime p | N, OR test N = Phi_15(b) directly!).
    chosen_b_N = []
    for b in range(20, 400):
        val = abs(int(Phi15(b)))
        facs = factorint(val)
        # Require squarefree and every prime factor p % 8 in (1, 7) so sqrt(2) exists mod N!
        if all(exp == 1 and pr > 5 and pr % 8 in (1, 7) for pr, exp in facs.items()) and len(facs) >= 2:
            chosen_b_N.append((b, val))
            if len(chosen_b_N) == 2:
                break
    assert len(chosen_b_N) >= 1
    b, N = chosen_b_N[-1]  # Larger b for cleaner asymptotic separation

    G_abs = (Phi15 + (T_var - Poly([b]))) ** 2 - 2 * (T_var - Poly([b])) ** 2
    assert ex.irreducible(G_abs)

    # Precompute roots of Phi_15 (8 unit-circle roots) and G_abs (8 roots at +sqrt(2), 8 roots at -sqrt(2))
    s2 = acb(2).sqrt()
    phi_roots = [rt for rt, _ in Phi15.complex_roots()]
    g_pos_roots = [rt for rt, _ in G_abs.complex_roots() if (cov_eval := sum(acb(int(Phi15[k])) * (rt**k) for k in range(9)) + (rt - b) * (1 + s2)).contains(0)]
    g_neg_roots = [rt for rt, _ in G_abs.complex_roots() if (sum(acb(int(Phi15[k])) * (rt**k) for k in range(9)) + (rt - b) * (1 - s2)).contains(0)]
    assert len(phi_roots) == 8 and len(g_pos_roots) == 8 and len(g_neg_roots) == 8

    out = []
    for m in (3, 4):
        D_abs = 32
        bits = N.bit_length()
        X0 = (3**D_abs) * (2 ** ((bits + 15) // 16))
        max_W1 = 0
        min_full_W2 = None
        min_full_F = None
        full_under_X0 = 0
        for z in itertools.product((-1, 0, 1), repeat=2 * m):
            if next((x for x in z if x), 0) <= 0:
                continue
            if not all(z[2 * j] != 0 or z[2 * j + 1] != 0 for j in range(m)):
                continue
            # Evaluate W1 = N_{K/Q}(sum_{j=0}^{m-1} (a_j + b_j*sqrt(2)) zeta_15^j)
            p_w1 = acb(1)
            for sgn in (s2, -s2):
                for rt in phi_roots:
                    val = acb(0)
                    for j in range(m):
                        val += (acb(z[2 * j]) + acb(z[2 * j + 1]) * sgn) * (rt**j)
                    p_w1 *= val
            w1_int = abs(int(p_w1.real.unique_fmpz()))
            if w1_int == 0:
                continue
            p_w2 = acb(1)
            for rt in g_pos_roots:
                val = acb(0)
                for j in range(m):
                    val += (acb(z[2 * j]) + acb(z[2 * j + 1]) * s2) * (rt**j)
                p_w2 *= val
            for rt in g_neg_roots:
                val = acb(0)
                for j in range(m):
                    val += (acb(z[2 * j]) - acb(z[2 * j + 1]) * s2) * (rt**j)
                p_w2 *= val
            w2_int = abs(int(p_w2.real.unique_fmpz()))
            if w2_int == 0:
                continue
            F_val = w1_int * w2_int
            max_W1 = max(max_W1, w1_int)
            if min_full_W2 is None or w2_int < min_full_W2:
                min_full_W2 = w2_int
            if min_full_F is None or F_val < min_full_F:
                min_full_F = F_val
            if F_val <= X0:
                full_under_X0 += 1

        out.append(
            {
                "m": m,
                "b": b,
                "N": N,
                "bits": bits,
                "max_Side1_norm_W1": max_W1,
                "log2_min_full_Side2_W2": round(math.log2(min_full_W2), 2),
                "predicted_log2_b_pow_2_m_minus_1": round(2 * (m - 1) * math.log2(b), 2),
                "log_min_full_F_over_log_N": round(math.log(min_full_F) / math.log(N), 4),
                "predicted_2_m_minus_1_over_8": round(2 * (m - 1) / 8.0, 4),
                "log2_X0": round(math.log2(X0), 2),
                "log2_min_full_F": round(math.log2(min_full_F), 2),
                "full_support_under_X0": full_under_X0,
            }
        )
    return out


def main() -> None:
    res = {
        "T1_two_hyperplane_and_witness": run_T1(),
        "T2_multicone_lagrange_seeds": run_T2(),
        "T3_diffuse_common_denom_divdiff": run_T3(),
        "T4_multiquadratic_unit_cancellation": run_T4(),
        "T5_cm_cyclotomic": run_T5(),
    }
    out_file = HERE / "results.json"
    out_file.write_text(json.dumps(res, indent=2) + "\n")
    print("DONE. Summary:")
    print("T1:", json.dumps({k: v for k, v in res["T1_two_hyperplane_and_witness"].items() if "sorted_abs" not in k}, indent=2))
    print("T2:", json.dumps(res["T2_multicone_lagrange_seeds"], indent=2))
    print("T3:", json.dumps(res["T3_diffuse_common_denom_divdiff"], indent=2))
    print("T4:", json.dumps(res["T4_multiquadratic_unit_cancellation"], indent=2))
    print("T5:", json.dumps(res["T5_cm_cyclotomic"], indent=2))


if __name__ == "__main__":
    main()
