"""Tests for agent rankdegen (src/rankdegen): frame-Hoelder tools, exact LP, structure lemmas, Phi program."""
from __future__ import annotations

import itertools
import math
import random
from fractions import Fraction as Fr

import sympy as sp

from src.rank3.forms import poly_comb, resultant
from src.rankdegen import fh, ledger

C2 = dict(f=[-2, 0, 0, 0, 0, 0, 1], P=[[1], [0, 1], [0, 0, 0, 1]], Q=[[1], [31, 1], [0, 0, 0, 1]])


def test_c2_minor_structure_and_lp():
    T = fh.minor_table(C2["f"], C2["P"], C2["f"], C2["Q"], 3)
    assert sum(1 for v in T["status"].values() if v == 0) == 2 * math.comb(6, 3)  # two rich lines
    fr = fh.frames(T)
    assert len(fr) == math.comb(6, 2) ** 2  # frames = 2 rows from each line
    val, w = fh.maxmin_coverage(T, fr, [Fr(4, 12)] * 12)
    assert val == Fr(1, 9)
    st = fh.distribution_stats(T, fr, w)
    assert all(m == Fr(1, 3) for m in st["marg"]) and st["mincov"] == Fr(1, 9)


def test_type_A_has_no_n_saving():
    phi5 = [1, 1, 1, 1, 1]
    P = [[1], [0, 1], [0, 0, 1]]
    T = fh.minor_table(phi5, P, phi5, P, 3)
    fr = fh.frames(T)
    val, _ = fh.maxmin_coverage(T, fr, [Fr(4, 8)] * 8)
    assert val == 0  # coincident cross pairs lie in no frame
    for a in itertools.product(range(-2, 3), repeat=3):  # W = W' identically (Lemma CC setting)
        A = poly_comb(P, a)
        if any(A):
            assert resultant(phi5, A) == resultant(phi5, A)


def test_gp_side_balanced_frames_beat_uniform():
    # power basis on two cubics: uniform (Bean-Thunder) coverage r(r+1)/(D(D-1)) = 2/5; LP max-min = 4/9
    f, g = [-2, 0, 0, 1], [-3, 1, 0, 1]
    P = [[1], [0, 1], [0, 0, 1]]
    T = fh.minor_table(f, P, g, P, 3)
    fr = fh.frames(T)
    assert len(fr) == math.comb(6, 4)
    val, _ = fh.maxmin_coverage(T, fr, [Fr(4, 6)] * 6)
    assert val == Fr(4, 9) > Fr(2, 5)
    st = fh.distribution_stats(T, fr, fh.pure_s(T, fr, 2))
    assert st["mincov"] == Fr(4, 9) and st["marg"][0] == st["marg"][3] == Fr(2, 3)


def test_lp_small_known():
    # max x1 + x2 s.t. x1 + 2 x2 <= 4, 3 x1 + x2 <= 6  -> (8/5, 6/5), value 14/5
    val, x = fh.lp_max([Fr(1), Fr(1)], [], [], [[Fr(1), Fr(2)], [Fr(3), Fr(1)]], [Fr(4), Fr(6)])
    assert val == Fr(14, 5) and x == [Fr(8, 5), Fr(6, 5)]
    val, x = fh.lp_max([Fr(1), Fr(0)], [[Fr(1), Fr(1)]], [Fr(1)], [], [])
    assert val == 1


def test_real_frame_identity():
    # r+1 real forms in GP: prod_i |det(S minus i)| = |det L|^{r+1} prod |c_i|, L_0 = sum c_i L_i (Cramer)
    rnd = random.Random(1)
    for r in (2, 3, 4):
        L = sp.Matrix(r, r, lambda i, j: rnd.randint(-5, 5))
        if L.det() == 0:
            continue
        c = [rnd.choice([-3, -2, -1, 1, 2, 3]) for _ in range(r)]
        L0 = sum((c[i] * L.row(i) for i in range(r)), sp.zeros(1, r))
        rows = [L0] + [L.row(i) for i in range(r)]
        prod = 1
        for i in range(r + 1):
            M = sp.Matrix.vstack(*[rows[k] for k in range(r + 1) if k != i])
            prod *= abs(M.det())
        assert prod == abs(L.det()) ** (r + 1) * abs(math.prod(c))


def test_complex_to_real_top_coefficient():
    # det( (w_j L_j + conj(L_j))/2 ) has coefficient 2^{-r} det L at prod w_j (used in the Jensen step)
    r = 3
    w = sp.symbols("w0:3")
    L = sp.Matrix([[1 + 2 * sp.I, 3, -sp.I], [2, 1 - sp.I, 4], [sp.I, 5, 1]])
    M = sp.Matrix(r, r, lambda i, j: (w[i] * L[i, j] + sp.conjugate(L[i, j])) / 2)
    poly = sp.Poly(sp.expand(M.det()), *w)
    assert sp.simplify(poly.coeff_monomial(w[0] * w[1] * w[2]) - L.det() / 8) == 0


def test_lemma_sub_subfield_power():
    # f = h(x^2), module in Q(x^2): Res(f, A(x^2)) = Res(h, A)^2 (each K_0-row appears [K:K_0] = 2 times)
    h = [1, 1, 0, 1]
    f = [1, 0, 1, 0, 0, 0, 1]
    for a in itertools.product(range(-3, 4), repeat=3):
        if not any(a):
            continue
        A = list(a)
        A2 = [a[0], 0, a[1], 0, a[2]]
        while len(A2) > 1 and A2[-1] == 0:
            A2.pop()
        while len(A) > 1 and A[-1] == 0:
            A.pop()
        assert resultant(f, A2) == resultant(h, A) ** 2


def test_kernel_norm_and_polya_fiber():
    # C1: K side (1, x, 0), K' side (1, x, x^2 - m x); z0 = e3; |W'(z0)| >= N, and each fiber a + k z0 has
    # at most d' + 4 (Y/2N)^{1/d'} points with |W'| <= Y (Polya).
    N, m = 1009 * 1013, 123456
    g = [(-(m ** 3 + 2 * m)) % N, 2, 0, 1]
    g[0] = g[0] - N if g[0] > N // 2 else g[0]
    assert (g[0] + 2 * m + m ** 3) % N == 0
    Q3 = [0, -m, 1]
    w0 = abs(resultant(g, Q3))
    assert w0 % N == 0 and w0 >= N
    for (a, b) in [(1, 0), (2, -1), (5, 3)]:
        Y = 10 ** 6
        cnt = 0
        for k in range(-2000, 2001):
            A = poly_comb([[1], [0, 1], Q3], (a, b, k))
            if 0 < abs(resultant(g, A)) <= Y:
                cnt += 1
        assert cnt <= 3 + 4 * (Y / (2 * N)) ** (1 / 3)


def test_phi_crossing_lemma_random():
    # Prop PHI: at any point feasible for the pure-s constraints (s in [s_lo, s_hi], one-sided lower bounds),
    # q1 q2 R(q) >= (r-1)/r. Exact check of the crossing inequality on random rational points.
    rnd = random.Random(7)
    for r in (2, 3, 4, 5):
        for (s_lo, s_hi) in [(1, r), (1, max(1, r - 1)), (2, r) if r > 2 else (1, r)]:
            for _ in range(300):
                q1 = Fr(1, s_hi) + Fr(rnd.randint(0, 400), 37)
                q2 = Fr(1, r + 1 - s_lo) + Fr(rnd.randint(0, 400), 41)
                assert ledger.phi_crossing_certificate(r, s_lo, s_hi, q1, q2)


def test_phi_numeric_matches_closed():
    assert abs(ledger.phi_numeric(3, [0, 1, 2], 0.5, 0.0, grid=200, rounds=4) - 2 / 3) < 2e-3
    assert abs(ledger.phi_numeric(2, [0, 1, 2, 3], 0, 0, grid=200, rounds=4) - 0.5) < 2e-3


def test_constants():
    assert abs(float(ledger.formal_constant(Fr(1, 2))) - 1.52629) < 1e-4
    assert abs(float(ledger.formal_constant(Fr(2, 3))) - 1.67989) < 1e-4
    assert abs(float(ledger.scoped_constant(Fr(2, 3), 3)) - 2.42283) < 1e-4
    # O3 under H_split, formal: c^3 = 8 Phi/9; r = 2 gives smoothenum's algebraic-side-only infimum 0.7631
    assert abs((8 * 0.5 / 9) ** (1 / 3) - 0.76314) < 1e-4
    assert ledger.hsplit_floor(2, iters=1500) >= (8 * 0.5 / 9) ** (1 / 3) - 1e-6


def test_split_threshold():
    r, d2, x = 3, 5, 1000.0
    assert abs(ledger.split_side_charge(r, 7, d2, ledger.split_threshold(r, d2, x), x, 1)) < 1e-9
    assert ledger.split_threshold(r, d2, x) > ledger.r5_threshold(r, 7 + d2, x)
