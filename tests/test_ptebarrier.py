"""Tests for src/ptebarrier (Lemma R rigidity, Theorem W invariant, Prop S support lemma)."""
import math
import random

from src.ptebarrier.rigidity import (
    from_roots, mul, sub, add, trim, evalp, height, divmod_monic, floor_log2_delta,
    regime_ok, h_R, lemma_R_certificate, psi_invariant, view_a_small_value_set,
    prop_S_hypothesis, injection_bound, int_roots,
)
from src.ptebarrier.enum import split_polys, group_by_top
from src.ptebarrier.run_experiment import split_over_Z


def test_poly_basics():
    p = from_roots([1, 2, 3])
    assert p == [-6, 11, -6, 1]
    q, r = divmod_monic(mul(p, [1, 0, 1]), [1, 0, 1])
    assert r == [] and q == p
    assert evalp(p, 4) == 6
    assert int_roots(p) == {1, 2, 3}


def test_floor_log2_delta_exact():
    for B in range(1, 9):
        for K in range(1, 6):
            for h in (1, 2, 3, 7, 100):
                v = 2 ** (K - 1) * math.sqrt(K) * h * B ** (K - 1)
                L = floor_log2_delta(B, K, h)
                assert 2 ** L <= v * (1 + 1e-12) and v < 2 ** (L + 1) * (1 + 1e-12)


def test_h_R_monotone_and_consistent():
    for B in (1, 2, 5):
        for K in (1, 2, 3):
            for M in range(K + 1, 40):
                hr = h_R(B, M, K)
                if hr:
                    assert regime_ok(M, B, K, hr) and not regime_ok(M, B, K, hr + 1)
                else:
                    assert not regime_ok(M, B, K, 1)


def test_certificate_random_pairs():
    """Every step of the proof of Lemma R holds on random pairs (j = 0 and j = 1)."""
    rng = random.Random(1)
    for _ in range(400):
        B = rng.randint(1, 6)
        M = rng.randint(2, 9)
        A = [rng.randint(-B, B) for _ in range(M)]
        A2 = list(A)
        for _ in range(rng.randint(1, 3)):
            A2[rng.randrange(M)] = rng.randint(-B, B)
        Q = trim([rng.randint(-2, 2), rng.randint(-2, 2)]) or [1]
        Q2 = trim([rng.randint(-2, 2), rng.randint(-2, 2)]) or [1]
        for (q1, q2) in (([1], [1]), (Q, Q2)):
            c = lemma_R_certificate(A, A2, q1, q2, B=B)
            if c["distinct"]:
                assert all(c["steps"]), c
                assert c["consistent"], c


def test_codex_shared_defect_example_out_of_regime():
    # (X+1)(X+6)^2 - (X+2)^2(X+9) = 8X: a shared defect is possible at M = 3 (out of regime)
    P1, P2 = from_roots([-1, -6, -6]), from_roots([-2, -2, -9])
    assert sub(P1, P2) == [0, 8]
    c = lemma_R_certificate([-1, -6, -6], [-2, -2, -9], B=9)
    assert c["consistent"] and not c["regime"]


def test_e1_small_grid_no_violation():
    for B, Mmax in ((1, 14), (2, 10), (3, 8)):
        for M in range(2, Mmax + 1):
            ent = list(split_polys(B, M))
            assert len(ent) == math.comb(2 * B + M, M)
            for K in (1, 2, 3):
                if K >= M:
                    continue
                for mem in group_by_top(ent, K).values():
                    for i in range(len(mem)):
                        for k in range(i + 1, len(mem)):
                            h = max(abs(x - y) for x, y in zip(mem[i][1], mem[k][1]))
                            assert not regime_ok(M, B, K, h), (B, M, K, mem[i], mem[k])
                            if K == 1:
                                assert not set(mem[i][0]) & set(mem[k][0])
                assert max(len(m) for m in group_by_top(ent, K).values()) <= injection_bound(2 * B + 1, K)


def test_theorem_W_forced_structure_has_one_invariant():
    """Members hP = FQ + r of the forced shape P = Phi0*S, r = r0*S share h*P/r."""
    rng = random.Random(7)
    Phi0 = from_roots([3, -1, 4, 1, -5, 2, 2, 0])
    r0 = [2, -1]
    F = sub(Phi0, r0)
    invs = set()
    for _ in range(5):
        S = from_roots([rng.randint(-6, 6) for _ in range(2)])
        P = mul(Phi0, S)
        r = mul(r0, S)
        assert sub(P, r) == mul(F, S)  # exact identity P = F*S + r
        invs.add(str(psi_invariant([1], P, r)))
    assert len(invs) == 1


def test_prop_S_support_bound():
    rng = random.Random(3)
    checked = 0
    for _ in range(300):
        d = rng.randint(1, 4)
        m = rng.randint(-20, 20)
        Y = rng.randint(1, 6)
        F = [rng.randint(-10**9, 10**9) for _ in range(d)] + [1]
        if evalp(F, m) == 0:
            continue
        rho, H = rng.randint(0, 2), rng.randint(1, 5)
        if prop_S_hypothesis(F, m, Y, rho, H):
            S = view_a_small_value_set(F, m, Y, rho, H)
            assert len(S) <= d
            checked += 1
    assert checked > 50


def test_prop_S_on_constructed_member():
    # P split near m, F := P - c so P = F*1 + c: every root of P lies in S (|F(a)| = |c|)
    m, roots, c = 50, [47, 48, 49, 49], 3
    P = from_roots(roots)
    F = sub(P, [c])
    assert prop_S_hypothesis(F, m, 3, 0, c) is (abs(evalp(F, m)) > 1 * c * 1 * 5 * 3 ** 4)
    S = view_a_small_value_set(F, m, 3, 0, c)
    assert set(roots) <= set(S)


def test_dickson_example_split_member():
    F = [-7, -1, 9, 0, -6, 0, 1]
    p = list(F)
    p[1] += 1
    p[0] += 3
    assert split_over_Z(p) == [-2, -1, -1, 1, 1, 2]
    assert split_over_Z(F) is None
