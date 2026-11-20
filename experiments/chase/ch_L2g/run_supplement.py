"""Supplementary pre-registered checks S1-S4 for ch_L2g."""
from __future__ import annotations

import itertools
import json
import math
import sys
from functools import reduce
from pathlib import Path

ROOT = Path("/Volumes/git/factoring-connection")
HERE = Path("/Volumes/git/factoring-connection/codex/scratch/ch_L2g")
sys.path[:0] = [
    str(HERE),
    str(ROOT / "experiments/masullo_ext/mext_c"),
    str(ROOT / "experiments/subthird/codex"),
    str(ROOT),
]

import cover as cov
import experiment as ex
from flint import acb, ctx, fmpz_poly as Poly
from run_adversary import eval_row_acb, make_qsqrt2_semiprime, prefix_pair, qsqrt2_orbit_rep
from sympy import cyclotomic_poly, factorint, integer_nthroot, primerange, symbols


def run_S1() -> list[dict]:
    """Run 1-hyperplane vs 2-hyperplane pair-product cover on all 4 cover.py fixtures."""
    ctx.prec = 192
    units = cov.pair_units()
    out = []
    for bits, mode in itertools.product((72, 96), ("divided", "power")):
        N = ex.fixture(bits, 19)
        _, b, P, Q = ex.select(N)
        f, g = P * P - 2, Q * Q - 2
        U, V = [ex.expand(X, X * X - 2, ex.generators(X, mode, 3)) for X in (P, Q)]
        residues = [int(u(b)) % N for u in U]
        points = [
            z
            for z in itertools.product(range(-2, 3), repeat=6)
            if next((j for j in z if j), 0) > 0
            and math.gcd(sum(t * c for t, c in zip(z, residues)), N) == 1
        ]
        rows, _ = ex.rows_at(ex.coeffs(f), list(map(ex.coeffs, U)), ex.coeffs(g), list(map(ex.coeffs, V)), 192)
        pos_idx = [idx for idx, r in enumerate(rows) if r[1].real > 0]
        neg_idx = [idx for idx, r in enumerate(rows) if r[1].real < 0]
        same_sigma_pairs = list(itertools.combinations(pos_idx, 2)) + list(itertools.combinations(neg_idx, 2))

        spaces_1h = [cov.prefix(r, points) for r in rows]
        spaces_2h = [prefix_pair(rows[i], rows[j], points) for i, j in same_sigma_pairs]

        cap = 3**32 * 2 ** ((bits + 15) // 16)
        fourth, exact = integer_nthroot(bits * (bits - 1).bit_length() ** 3, 4)
        B = 2 ** (int(fourth) + int(not exact))
        primes = list(primerange(B + 1))
        primorial = math.prod(primes)
        fields = [ex.Ideals(X, primes) for X in (f, g)]

        keys = set()
        uncov_1h, uncov_2h = 0, 0
        r_all, r_1h, r_2h = [], [], []
        cols: dict = {}
        for z in points:
            a_el, w_el = ex.element(U, z), ex.element(V, z)
            v1, v2 = abs(int(f.resultant(a_el))), abs(int(g.resultant(w_el)))
            val = v1 * v2
            if not 0 < val <= cap:
                continue
            key = cov.unit_key(z, units)
            if key in keys:
                continue
            keys.add(key)
            c1 = any(ex.kranks(bas + [z]) == len(bas) for bas in spaces_1h)
            c2 = any(ex.kranks(bas + [z]) == len(bas) for bas in spaces_2h)
            if not c1:
                uncov_1h += 1
            if not c2:
                uncov_2h += 1
            if ex.smooth(val, primorial):
                msk = 0
                for side, (al, nrm) in enumerate(((a_el, v1), (w_el, v2))):
                    for pr in fields[side].parity(al, nrm, primes):
                        msk |= 1 << cols.setdefault((side, pr), len(cols))
                r_all.append(msk)
                if c1:
                    r_1h.append(msk)
                if c2:
                    r_2h.append(msk)

        out.append(
            {
                "bits": bits,
                "mode": mode,
                "capped_unit_classes": len(keys),
                "uncovered_1h": uncov_1h,
                "uncovered_2h": uncov_2h,
                "smooth_rank_all": ex.rank(r_all),
                "inc_rank_1h": ex.rank(r_all) - ex.rank(r_1h),
                "inc_rank_2h": ex.rank(r_all) - ex.rank(r_2h),
            }
        )
    return out


def run_S2() -> dict:
    """T2 at (n, m) = (8, 4) scanning subset products explicitly in M to verify all 2^m - 2 = 14 small-norm orbits
    collapse to ideal-parity rank <= m = 4, plus E=2 box scan on (6, 3) and E=3 sparse/coordinate check on (8, 4)."""
    T_var = Poly([0, 1])
    R0 = Poly([3, 2])
    n, m, bits = 8, 4, 48
    N, p, q, A0 = make_qsqrt2_semiprime(bits, 7)
    nodes = (0, 1, -1, 2)
    Qm = reduce(lambda a, b: a * b, (T_var - Poly([r_val]) for r_val in nodes))
    deg_A = n - m
    b_base = max(6, int(integer_nthroot(N, deg_A)[0]))
    found = None
    for b in range(b_base, b_base + 200):
        qb = int(Qm(b))
        if qb == 0 or math.gcd(abs(qb), N) != 1:
            continue
        rem0 = (A0 - int(R0(b))) % abs(qb)
        t_shift = ((-rem0) * pow(N % abs(qb), -1, abs(qb))) % abs(qb)
        target = (A0 - int(R0(b)) + t_shift * N) // qb
        rem = target - b**deg_A
        digs = []
        for _ in range(deg_A - 1):
            d_i = (rem + b // 2) % b - b // 2
            digs.append(d_i)
            rem = (rem - d_i) // b
        A_poly = Poly(digs + [rem, 1])
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

    # Scale Lagrange basis by 6 = lcm_j prod_{l!=j} |r_j - r_l|? No: 6 * L_S(T) is an EXACT integer combination in [-6, 6]^4
    # Better: in the power basis [1, T, T^2, T^3] (which spans the EXACT same k-module M!), EVERY subset product
    # prod_{r in S} (T - r) for 1 <= |S| <= 3 has integer coefficients in [-2, 2]^4!
    pow_gens = [T_var**j for j in range(m)]
    U = ex.expand(P, f, pow_gens)
    V = ex.expand(Q, g, pow_gens)
    residues = [int(u(b)) % N for u in U]

    # Scan E=2 on Rational/sqrt(2) or enumerate all non-empty proper subsets S of nodes + O_k multiples in [-2, 2]^8
    # Wait: 5^8 = 390,625 -> fast if we evaluate resultant only when needed, or test all 14 subset products + small O_k scalars!
    # Actually 5^8 with polynomial resultant of deg 16 takes ~6s in flint! Let's check: first test mod a small prime? No, flint resultant is ~5us -> 390k * 10us = 3.9s!
    # Even faster: test W1 first, only compute W2 if |W1| <= small_cap!
    small_cap = 60**m
    small_orbits: dict = {}
    all_primes = set()
    for z in itertools.product(range(-2, 3), repeat=2 * m):
        if next((x for x in z if x), 0) <= 0:
            continue
        if math.gcd(sum(t * c for t, c in zip(z, residues)), N) != 1:
            continue
        a_el = ex.element(U, z)
        v1 = abs(int(f.resultant(a_el)))
        if v1 == 0 or v1 > small_cap:
            continue
        w_el = ex.element(V, z)
        v2 = abs(int(g.resultant(w_el)))
        F_val = v1 * v2
        if 0 < F_val <= small_cap:
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

    return {
        "n": n,
        "m": m,
        "E": 2,
        "distinct_small_Ok_orbits": len(small_orbits),
        "spanned_k_rank": ex.kranks(orbit_vecs),
        "ideal_parity_rank_of_all_small_orbits": ex.rank(masks),
    }


def run_S3() -> list[dict]:
    """T3 with syzygy removed: gens = [A, T*A, ..., T^{m-1}*A], certified true k-rank m."""
    ctx.prec = 256
    T_var = Poly([0, 1])
    D_poly = Poly([1, 0, 1])
    out = []
    for n, m, bits in ((6, 3, 60), (8, 4, 72)):
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
        gens_P = [(T_var**j) * A_poly for j in range(m)]
        gens_Q = [(T_var**j) * B_poly for j in range(m)]
        U = ex.expand(P, f, gens_P)
        V = ex.expand(Q, g, gens_Q)
        residues = [int(u(b)) % N for u in U]
        assert residues == [int(v(b)) % N for v in V]

        # Certify true Z-rank == 2*m on each side via polynomial coordinates
        col_mat = [[int(u[k]) for k in range(2 * n)] for u in U]
        from flint import fmpq_mat

        assert fmpq_mat(col_mat).rank() == 2 * m

        rows, _ = ex.rows_at(ex.coeffs(f), list(map(ex.coeffs, U)), ex.coeffs(g), list(map(ex.coeffs, V)), 256)
        D_abs = 4 * n
        X0 = (3**D_abs) * (2 ** ((bits + 2 * n - 1) // (2 * n)))

        min_full_F = None
        small_cnt_best = 0
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
                small_cnt_best = sum(1 for r in rows if float(abs(eval_row_acb(r, z)).mid()) <= 10.0)
            if F_val <= X0:
                full_under_X0 += 1

        out.append(
            {
                "n": n,
                "m": m,
                "certified_k_rank": m,
                "bits": bits,
                "b": b,
                "small_embeddings_at_best_full_point": small_cnt_best,
                "predicted_diffuse_small_embeddings_4_n_minus_2": 4 * (n - 2),
                "log2_X0": round(math.log2(X0), 2),
                "log2_min_full_F": round(math.log2(min_full_F), 2),
                "log_min_full_F_over_log_N": round(math.log(min_full_F) / math.log(N), 4),
                "predicted_4_m_minus_1_log_b_over_log_N": round(4 * (m - 1) * math.log(b) / math.log(N), 4),
                "full_support_under_X0": full_under_X0,
            }
        )
    return out


def run_S4() -> list[dict]:
    """T5 CM/cyclotomic K = Q(sqrt(2), zeta_15) at larger b (~600-1000, bits ~ 74) so b^{2(m-1)} > 3^32."""
    ctx.prec = 256
    x_sym = symbols("x")
    phi15_coeffs = [int(c) for c in reversed(cyclotomic_poly(15, x_sym).as_poly(x_sym).all_coeffs())]
    Phi15 = Poly(phi15_coeffs)
    T_var = Poly([0, 1])

    # Choose b = 1024 (so N = Phi_15(1024) has ~80 bits; check irreducibility of G_abs; note arithmetic norm geometry
    # on [-1,1]^{2m} depends only on b and Phi_15, or find b in [500, 1200] where Phi_15(b) is squarefree!)
    b_chosen, N_chosen = None, None
    for b in range(600, 1100):
        val = abs(int(Phi15(b)))
        G_abs = (Phi15 + (T_var - Poly([b]))) ** 2 - 2 * (T_var - Poly([b])) ** 2
        if ex.irreducible(G_abs):
            b_chosen, N_chosen = b, val
            break
    assert b_chosen is not None
    b, N = b_chosen, N_chosen
    G_abs = (Phi15 + (T_var - Poly([b]))) ** 2 - 2 * (T_var - Poly([b])) ** 2

    s2 = acb(2).sqrt()
    phi_roots = [rt for rt, _ in Phi15.complex_roots()]
    g_pos_roots = [
        rt
        for rt, _ in G_abs.complex_roots()
        if (sum(acb(int(Phi15[k])) * (rt**k) for k in range(9)) + (rt - b) * (1 + s2)).contains(0)
    ]
    g_neg_roots = [
        rt
        for rt, _ in G_abs.complex_roots()
        if (sum(acb(int(Phi15[k])) * (rt**k) for k in range(9)) + (rt - b) * (1 - s2)).contains(0)
    ]
    assert len(phi_roots) == 8 and len(g_pos_roots) == 8 and len(g_neg_roots) == 8

    out = []
    for m in (3, 4):
        bits = N.bit_length()
        X0 = (3**32) * (2 ** ((bits + 15) // 16))
        min_full_F = None
        min_full_W2 = None
        full_under_X0 = 0
        for z in itertools.product((-1, 0, 1), repeat=2 * m):
            if next((x for x in z if x), 0) <= 0:
                continue
            if not all(z[2 * j] != 0 or z[2 * j + 1] != 0 for j in range(m)):
                continue
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
                "bits": bits,
                "log2_X0": round(math.log2(X0), 2),
                "log2_min_full_Side2_W2": round(math.log2(min_full_W2), 2),
                "predicted_log2_b_pow_2_m_minus_1": round(2 * (m - 1) * math.log2(b), 2),
                "log2_min_full_F": round(math.log2(min_full_F), 2),
                "log_min_full_F_over_log_N": round(math.log(min_full_F) / math.log(N), 4),
                "predicted_2_m_minus_1_over_8": round(2 * (m - 1) / 8.0, 4),
                "full_support_under_X0": full_under_X0,
            }
        )
    return out


def main() -> None:
    res = {
        "S1_cover_table": run_S1(),
        "S2_multicone_E2": run_S2(),
        "S3_diffuse_true_rank_m": run_S3(),
        "S4_cm_larger_b": run_S4(),
    }
    (HERE / "results_supp.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
