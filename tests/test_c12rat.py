"""Tests for c12rat (notes/c12rat_results.md): exact checks of the finite lemmas LS, MULT, TF, BT-S, CP, TYPE and
of the E3 decomposition identities."""
from math import gcd, isqrt

import numpy as np
import pytest

from src.c12rat import bts, e3, geom
from src.pc3 import toy as pc3toy


def test_shape_cdf_bound_exact_small_q():
    for q in (101, 1009):
        sig = np.array([geom.shape(q, r)[0] for r in range(q)])
        for s in (1 / 64, 1 / 16, 1 / 4, 1.0):
            assert (sig <= s).mean() <= geom.shape_cdf_bound(s, q) + 1e-12


def test_box_det_mult_thin_bounds():
    rng = np.random.default_rng(5)
    for q in (503, 1009, 10007):
        for _ in range(6):
            r = int(rng.integers(0, q))
            for S in (5, 8):
                ij, ab = geom.box_arrays(q, r, S)
                iu, ju, D, det, kap = geom.pair_stats(ij, ab)
                assert np.all(det == q * D)
                assert np.all(D > 0)  # distinct classes are non-proportional
                assert geom.mult_max(D) <= 5 * S * S + S
                sigma, t2 = geom.shape(q, r)
                n = len(ij)
                for K in (1.5, 4, 16, 64):
                    assert (kap > K).mean() <= geom.thin_bound(S, n, t2, K) + 1e-12
                # kappa bound used in Lemma TF
                assert np.all(kap <= 2 * np.sqrt(3) * S * S / (t2 * D) + 1e-9)


def test_reduced_norm_comparison():
    rng = np.random.default_rng(6)
    for _ in range(200):
        q = 10007
        r = int(rng.integers(0, q))
        e1, e2 = geom.label_basis(q, r)
        l1, l2 = np.hypot(*e1), np.hypot(*e2)
        for _ in range(20):
            i, j = (int(v) for v in rng.integers(-30, 31, size=2))
            P = (i * e1[0] + j * e2[0], i * e1[1] + j * e2[1])
            n2 = P[0] ** 2 + P[1] ** 2
            base = i * i * l1 * l1 + j * j * l2 * l2
            assert 0.5 * base - 1e-6 <= n2 <= 1.5 * base + 1e-6
        assert l1 * l2 <= 2 / np.sqrt(3) * q + 1e-6


def test_bts_inequality_bruteforce():
    rng = np.random.default_rng(7)
    X, y = 2_000_000, 30
    smooth, _ = pc3toy.build_tables(X, y, y + 1)
    for _ in range(300):
        D = int(rng.integers(1, 400))
        a = int(rng.integers(0, D)) if D > 1 else 0
        if gcd(a, D) != 1:
            continue
        Z = int(rng.integers(D, 200 * D))
        X0 = int(rng.integers(1, X - Z - 1))
        W = int(rng.integers(1, max(2, min(X0, Z // D + 2))))
        if W > X0:
            continue
        L = bts.lhs(smooth, X0, Z, a, D)
        brute = sum(1 for n in range(X0 + 1, X0 + Z + 1) if n % D == a % D and smooth[n])
        assert L == brute
        assert L <= bts.rhs(smooth, X0, Z, D, y, W) <= bts.rhs_analytic(smooth, Z, D, y, W) + 1e-9


def _units(M):
    return [r for r in range(M) if gcd(r, M) == 1]


def test_coupled_pairs_lemma_CP():
    """Lemma CP: C diagonal-invariant in (Z/M)^x squared; error <= (|C|/phi^2) sqrt(sum_{chi!=1}|a^|^2 sum|b^|^2),
    where sum_{chi != chi0}|a^(chi)|^2 = phi * sum_r |a(r) - A/phi|^2 (Parseval on the group)."""
    rng = np.random.default_rng(8)
    for M in (12, 35, 60, 97):
        U = _units(M)
        phi = len(U)
        idx = {u: k for k, u in enumerate(U)}
        for _ in range(20):
            # C = union of diagonal orbits: pairs (u, u*c) for c in a random subset of units
            cs = [c for c in U if rng.random() < 0.3] or [1]
            C = {(u, (u * c) % M) for u in U for c in cs}
            n = 400
            m1 = rng.integers(1, 10_000, size=n)
            m2 = rng.integers(1, 10_000, size=n)
            m1 = m1[np.gcd(m1, M) == 1]
            m2 = m2[np.gcd(m2, M) == 1]
            a = rng.random(len(m1))
            b = rng.random(len(m2))
            ar = np.zeros(phi)
            br = np.zeros(phi)
            np.add.at(ar, [idx[int(v) % M] for v in m1], a)
            np.add.at(br, [idx[int(v) % M] for v in m2], b)
            lhs = sum(ar[idx[r1]] * br[idx[r2]] for (r1, r2) in C)
            main = len(C) * ar.sum() * br.sum() / phi ** 2
            Ea = phi * ((ar - ar.sum() / phi) ** 2).sum()
            Eb = phi * ((br - br.sum() / phi) ** 2).sum()
            assert abs(lhs - main) <= len(C) / phi ** 2 * np.sqrt(Ea * Eb) + 1e-9


def test_parseval_cyclic():
    """sum over non-principal characters of |a^|^2 = phi * sum_r |a(r) - A/phi|^2, explicit characters mod a prime."""
    p = 31
    g = next(x for x in range(2, p) if all(pow(x, (p - 1) // f, p) != 1 for f in (2, 3, 5)))
    log = {pow(g, k, p): k for k in range(p - 1)}
    rng = np.random.default_rng(9)
    a = rng.random(p - 1)  # a(r) for r = g^k
    tot = 0.0
    for j in range(1, p - 1):
        chi = np.exp(2j * np.pi * j * np.arange(p - 1) / (p - 1))
        tot += abs((a * chi).sum()) ** 2
    assert abs(tot - (p - 1) * ((a - a.mean()) ** 2).sum()) < 1e-8


def test_type_coupling_sets_TYPE():
    """Lemma TYPE: for n = (G'(g), G''(g)) with g P0-primitive, the reduced pair (m1, m2) = (n1/s1, n2/s2) mod M_tau
    ranges exactly over the diagonal-invariant set C_tau; the fibre count over m1 depends only on the type."""
    P0 = 3
    small = (2, 3)
    for (Pp, Pq) in [((1, 2), (5, 2)), ((3, 4), (-1, 8)), ((2, 3), (8, 3)), ((1, 6), (5, 6))]:
        a1, b1 = Pp
        a2, b2 = Pq
        D = abs(a1 * b2 - a2 * b1)
        T = 4 * 9 * D * 6
        seen = {}
        for g0 in range(1, T + 1):
            for g1 in range(1, T + 1):
                if g0 % 2 == 0 and g1 % 2 == 0:
                    continue
                if g0 % 3 == 0 and g1 % 3 == 0:
                    continue
                n1 = g1 * a1 - g0 * b1  # signed: the lattice relation is between signed values
                n2 = g1 * a2 - g0 * b2
                if n1 == 0 or n2 == 0:
                    continue
                s1 = s2 = 1
                e1 = []
                for p in small:
                    v = 0
                    while n1 % p ** (v + 1) == 0:
                        v += 1
                    w = 0
                    while n2 % p ** (w + 1) == 0:
                        w += 1
                    if v > 4 or w > 4:
                        break
                    e1.append((v, w))
                    s1 *= p ** v
                    s2 *= p ** w
                else:
                    # mid/large primes of D: gcd with D for p > P0
                    d1 = gcd(n1 // s1, D)
                    d2 = gcd(n2 // s2, D)
                    tau = (tuple(e1), d1, d2)
                    delta = gcd(n1, D)
                    M = D // delta
                    for p in small:
                        if delta % p == 0:
                            M *= p
                    m1, m2 = (n1 // s1) % M if M > 1 else 0, (n2 // s2) % M if M > 1 else 0
                    seen.setdefault((tau, M), {}).setdefault(m1, set()).add(m2)
        for (tau, M), fib in seen.items():
            sizes = {len(v) for v in fib.values()}
            assert len(sizes) == 1, (Pp, Pq, tau, M, sizes)  # fibre size depends only on the type
            # diagonal invariance: realised set closed under (m1, m2) -> (u m1, u m2)
            pairs = {(k, v) for k, vs in fib.items() for v in vs}
            for u in _units(M)[:6]:
                for (x1, x2) in list(pairs)[:50]:
                    assert ((u * x1) % M, (u * x2) % M) in pairs


def test_e3_identities_and_cs():
    smooth, _ = pc3toy.build_tables(600_000, 60, 61)
    g = e3.all_primitive_g(120)
    out = e3.label_run(211, 17, 4, g, smooth, 60, 0.5, 16.0, prim_model=True)
    assert abs(out["C_nu"]) <= out["B1"] + 1e-12
    assert abs(out["phi_ex"] - out["C_nu"] - out["R_lab"]) < 1e-9


def test_local_ratio_formula():
    """ratio_p = (p+1)(2AB + (p-1)B^2)/(A + pB)^2 equals pair/single^2 for (n1, n2) uniform on F_p^2 minus 0 with
    weights A on p | n and B on p not dividing n."""
    for p in (11, 13, 101):
        for A in (0.5, 1.7, 3.0):
            B = (1 - A / p) / (1 - 1 / p)
            pts = [(x, z) for x in range(p) for z in range(p) if (x, z) != (0, 0)]
            w = lambda t: A if t == 0 else B
            pair = sum(w(x) * w(z) for x, z in pts) / len(pts)
            single = sum(w(x) for x, z in pts) / len(pts)
            ratio = (p + 1) * (2 * A * B + (p - 1) * B * B) / (A + p * B) ** 2
            assert abs(pair / single ** 2 - ratio) < 1e-12
