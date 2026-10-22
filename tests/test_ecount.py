"""Tests for src/ecount (agent ecount, phase 23)."""
import ast
import random
from math import gcd

import gmpy2
import pytest

from src.ecount import harness as H
from src.ecount.zn_curve import (Split, coppersmith_linear, du_affine_twist, du_malleability,
                                 du_projective_lemma_as_written, du_projective_twist_repaired, ec_add, ec_mul,
                                 find_divisor_near, sr2_supersingular, tw_factor, tw_trial)


def naive_count(a, b, p):
    c = 1
    for x in range(p):
        v = (x ** 3 + a * x + b) % p
        c += 1 if v == 0 else (2 if pow(v, (p - 1) // 2, p) == 1 else 0)
    return c


def naive_points(a, b, p):
    sq = {}
    for y in range(p):
        sq.setdefault(y * y % p, []).append(y)
    pts = [None]
    for x in range(p):
        for y in sq.get((x ** 3 + a * x + b) % p, []):
            pts.append((x, y))
    return pts


def test_operational_imports_clean():
    src = open("src/ecount/zn_curve.py").read()
    tree = ast.parse(src)
    mods = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods |= {a.name for a in node.names}
        elif isinstance(node, ast.ImportFrom):
            mods.add(node.module)
    assert not any(m and ("harness" in m or "sympy" in m or "cypari" in m) for m in mods)
    assert "factor(" not in src.replace("find_divisor", "").replace("_factor", "").replace("factor (", "")


@pytest.mark.parametrize("curve", list(H.CURVES))
def test_pari_counts_vs_naive(curve):
    a, b = H.CURVES[curve]
    for p in (101, 103, 1009, 1013, 4099):
        if H.disc_part(a, b) % p:
            assert H.card(a, b, p) == naive_count(a, b, p)
            d = H.local_data(a, b, p)
            assert d["A"] + d["At"] == 2 * p + 2
            at, bt = H.twist_coeffs(a, b, p)
            assert naive_count(at, bt, p) == d["At"]


def test_group_law_crt_consistency():
    """ec_add over Z/N agrees with the group law computed separately mod p and mod q (when it does not split)."""
    a, b = 2, 3
    p, q = 97, 103
    N = p * q
    Pp, Pq = naive_points(a, b, p), naive_points(a, b, q)
    rng = random.Random(0)

    def lift(u, v):
        if u is None and v is None:
            return None
        return (int(gmpy2.f_mod(u[0] * q * gmpy2.invert(q, p) + v[0] * p * gmpy2.invert(p, q), N)),
                int(gmpy2.f_mod(u[1] * q * gmpy2.invert(q, p) + v[1] * p * gmpy2.invert(p, q), N)))

    def add_mod(P, Q, r):
        if P is None:
            return Q
        if Q is None:
            return P
        try:
            R = ec_add(P, Q, a, r)
        except Split:
            raise AssertionError
        return R

    splits = 0
    for _ in range(3000):
        u1, v1 = rng.choice(Pp[1:]), rng.choice(Pq[1:])
        u2, v2 = rng.choice(Pp[1:]), rng.choice(Pq[1:])
        if rng.random() < 0.1:
            u2, v2 = u1, v1
        if rng.random() < 0.05:
            u2, v2 = (u1[0], (-u1[1]) % p), (v1[0], (-v1[1]) % q)
        P, Q = lift(u1, v1), lift(u2, v2)
        sp, sq_ = add_mod(u1, u2, p), add_mod(v1, v2, q)
        try:
            R = ec_add(P, Q, a, N)
        except Split as s:
            splits += 1
            assert s.factor in (p, q)
            assert (sp is None) != (sq_ is None) or (u1[0] == u2[0]) != (v1[0] == v2[0])
            continue
        if R is None:
            assert sp is None and sq_ is None
        else:
            assert (R[0] % p, R[1] % p) == sp and (R[0] % q, R[1] % q) == sq_
    assert splits > 0


@pytest.mark.parametrize("curve", list(H.NONCM))
@pytest.mark.parametrize("n", [32, 64, 128])
def test_tw_factors(curve, n):
    a, b = H.CURVES[curve]
    for seed in range(5):
        N, p, q = H.instance(curve, n, seed)
        M = H.count(a, b, p, q)
        res = tw_factor(N, a, b, M, random.Random(seed), 64)
        assert res["factor"] in (p, q)


def test_tw_random_M_does_not_split():
    a, b = H.CURVES["C2_37a1"]
    N, p, q = H.instance("C2_37a1", 64, 0)
    rng = random.Random(1)
    hits = 0
    for _ in range(200):
        fac, out = tw_trial(N, a, b, rng.randrange(N // 2, N), rng)
        hits += fac is not None and out == "split"
    assert hits == 0


def test_tw_rate_matches_rho_small():
    """Empirical per-trial structural split rate vs exact rho at small primes (exact enumeration over x0)."""
    a, b = H.CURVES["C1_11a1"]
    p, q = 1009, 1013
    N = p * q
    M = H.count(a, b, p, q)
    rho, bonus, exc = H.rho_exact(H.local_data(a, b, p), H.local_data(a, b, q))
    assert not exc
    rng = random.Random(3)
    T = 4000
    s = 0
    for _ in range(T):
        fac, out = tw_trial(N, a, b, M, rng)
        s += fac is not None
    emp = s / T
    assert abs(emp - (rho + bonus)) < 0.04 or emp > rho + bonus


def test_cm1_doubly_supersingular():
    a, b = H.CURVES["CM1_x3mx"]
    for seed in range(5):
        N, p, q = H.instance_cm1_33(48, seed)
        M = H.count(a, b, p, q)
        assert M == (p + 1) * (q + 1)
        assert sr2_supersingular(N, M) in (p, q)
        assert H.affine_count(a, b, p, q) == N
        assert H.twist_count(a, b, 2, p, q, True) == N
        rho, bonus, exc = H.rho_exact(H.local_data(a, b, p), H.local_data(a, b, q))
        assert exc and rho == 0  # twist trick proper cannot fire; any TW split here is the 2-adic (Miller) descent
        r = tw_factor(N, a, b, M, random.Random(seed), 16)
        assert r["factor"] in (None, p, q)


@pytest.mark.parametrize("curve", list(H.NONCM))
def test_du_affine_both_cases(curve):
    a, b = H.CURVES[curve]
    for n in (32, 64):
        for seed in range(4):
            N, p, q = H.instance(curve, n, seed)
            Ea = H.affine_count(a, b, p, q)
            d = 2
            while gmpy2.jacobi(d, N) != -1:
                d += 1
            assert du_affine_twist(N, Ea, H.twist_count(a, b, d, p, q, True), d) in (p, q)
            d1 = 2
            while not (gmpy2.legendre(d1, p) == -1 and gmpy2.legendre(d1, q) == -1):
                d1 += 1
            assert du_affine_twist(N, Ea, H.twist_count(a, b, d1, p, q, True), d1) in (p, q)


def test_du_projective_lemma_needs_PQ():
    a, b = H.CURVES["C2_37a1"]
    N, p, q = H.instance("C2_37a1", 48, 1)
    d1 = 2
    while not (gmpy2.legendre(d1, p) == -1 and gmpy2.legendre(d1, q) == -1):
        d1 += 1
    E, Ed = H.count(a, b, p, q), H.twist_count(a, b, d1, p, q, False)
    r1, r2 = du_projective_lemma_as_written(N, E, Ed, (p + 1) * (q + 1))
    others = {H.twist_count(a, b, d, p, q, False) for d in range(2, 60)} - {E, Ed}
    assert {r1, r2} <= others | {r1, r2} and r1 + r2 + E + Ed == 4 * (p + 1) * (q + 1)
    # with the wrong PQ (e.g. N + 1) the "roots" are not the twist counts
    try:
        s1, s2 = du_projective_lemma_as_written(N, E, Ed, N + 1)
        assert {s1, s2} != {r1, r2}
    except ValueError:  # negative discriminant: no real roots at all
        pass


@pytest.mark.parametrize("curve", list(H.NONCM))
def test_du_projective_repaired(curve):
    a, b = H.CURVES[curve]
    for n in (32, 64):
        for seed in range(3):
            N, p, q = H.instance(curve, n, seed)
            d = 2
            while gmpy2.jacobi(d, N) != -1:
                d += 1
            g = du_projective_twist_repaired(N, H.count(a, b, p, q), H.twist_count(a, b, d, p, q, False), d)
            assert g in (p, q)


def test_coppersmith_known_msb():
    rng = random.Random(7)
    from src.common.instances import balanced_semiprime
    for seed in range(3):
        inst = balanced_semiprime(96, seed)
        p, q = inst.factors
        N = inst.N
        err = rng.randrange(-2 ** 20, 2 ** 20)
        g = find_divisor_near(N, q + err, 2 ** 20, (96 - 2) / 2.0)
        assert g in (p, q)


def test_du_malleability_small():
    a, b = H.CURVES["C2_37a1"]
    for seed in range(3):
        N, p, q = H.instance("C2_37a1", 48, seed)
        M = H.count(a, b, p, q)
        r = du_malleability(N, H.factorint(M))
        assert r["factor"] in (p, q)


def test_chebotarev_l2_value():
    assert abs(H.chebotarev_adv_2() - 2 / 9) < 1e-12


@pytest.mark.parametrize("n", [64, 96])
def test_du_projective_repaired_lattice_path(n):
    a, b = H.CURVES["C2_37a1"]
    for seed in range(3):
        N, p, q = H.instance("C2_37a1", n, seed)
        d = 2
        while gmpy2.jacobi(d, N) != -1:
            d += 1
        g = du_projective_twist_repaired(N, H.count(a, b, p, q), H.twist_count(a, b, d, p, q, False), d, small_k=0)
        assert g in (p, q)
