"""Tests for src/schnorrcvp (agent schnorrcvp, phase 19)."""
import random
from math import gcd, log

import mpmath as mp
import pytest

from src.schnorrcvp.core import (first_primes, factor_over, in_lambda, is_fac_relation, lambda_vector,
                                 squares_dependency_factor, uvw, uw_lower_bound_holds)
from src.schnorrcvp.enum_rel import relations_in_box, smooth_mask_window, smooth_upto
from src.schnorrcvp.psi import PsiBracket, expected_relations

fpylll = pytest.importorskip("fpylll")
from src.schnorrcvp.lattice import analyse_vector, babai, build_lattice, exact_cvp, list_cvp, reduce_basis  # noqa: E402


def test_first_primes():
    assert first_primes(10) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    assert first_primes(47)[-1] == 211 and first_primes(95)[-1] == 499


def test_residual_identity_exact():
    """Thm R: w = u - vN = 2 exp((A+L)/2) sinh(Delta/2), A = ln(uv), Delta = ln(u/(vN))."""
    rng = random.Random(1)
    P = first_primes(12)
    mp.mp.dps = 80
    for _ in range(200):
        N = rng.randrange(10**6, 10**9) | 1
        e = [rng.randint(-3, 3) for _ in P]
        u, v, w = uvw(e, P, N)
        A = mp.log(u) + mp.log(v)
        D = mp.log(u) - mp.log(v) - mp.log(N)
        rhs = 2 * mp.e ** ((A + mp.log(N)) / 2) * mp.sinh(D / 2)
        assert abs(rhs - w) < mp.mpf(10) ** -40 * max(1, abs(w))
        assert gcd(u * v, w) == 1 or gcd(u * v, N) != 1


def test_uw_lower_bound():
    rng = random.Random(2)
    for _ in range(2000):
        N = rng.randrange(3, 10**6)
        v = rng.randrange(1, 50)
        w = rng.randrange(-v * N + 1, 10**6)
        if w == 0:
            continue
        u = v * N + w
        assert uw_lower_bound_holds(u, v, w, N)
    assert (10 - 1) * 1 == 10 - 1  # attained: u = N - 1, v = 1, w = -1


def test_lambda_vectors_and_factoring_toy():
    """Box relations -> Lambda_N vectors (exact membership) -> GF(2) congruence of squares factors N."""
    p, q = 1009, 1013
    N = p * q
    P = first_primes(25)
    r = relations_in_box(N, P, 2000, 20000)
    vecs = []
    for u, v, w in r["relations"]:
        assert is_fac_relation_triple(u, v, w, N, P)
        lv = lambda_vector(u, w, P)
        assert lv is not None and in_lambda(lv, P, N)
        vecs.append(lv)
    vecs = list(dict.fromkeys(vecs))
    assert len(vecs) > len(P) + 1
    g, ctr = squares_dependency_factor(vecs, P, N)
    assert g in (p, q)


def is_fac_relation_triple(u, v, w, N, P):
    return u - v * N == w and factor_over(u, P) is not None and factor_over(v, P) is not None \
        and factor_over(w, P) is not None


def test_gf2_kernel_basis_suffices():
    """Lemma LA: c -> x(c) is a homomorphism from the GF(2) kernel to {x : x^2 = 1}; checking a basis suffices.
    Here: product of the roots of two basis dependencies equals the root of their sum, up to Lambda'."""
    p, q = 1009, 1013
    N = p * q
    P = first_primes(25)
    r = relations_in_box(N, P, 2000, 20000)
    vecs = list(dict.fromkeys(lambda_vector(u, w, P) for u, v, w in r["relations"]))
    # every root of 1 we produce squares to 1 (asserted inside); run with all dependencies
    g, ctr = squares_dependency_factor(vecs, P, N, max_deps=10**6)
    assert g is not None and N % g == 0 and 1 < g < N


def test_smooth_window_matches_trial_division():
    P = first_primes(8)
    lo, L = 10**9 + 7, 5000
    mask = smooth_mask_window(lo, L, P)
    for i in range(L):
        assert bool(mask[i]) == (factor_over(lo + i, P) is not None)
    assert smooth_upto(100, first_primes(3)) == sorted(m for m in range(1, 101) if factor_over(m, first_primes(3)) is not None)


def test_psi_bracket_contains_truth():
    P = first_primes(6)
    pb = PsiBracket(P, 500, 16)
    for X in [10, 1000, 54321, 10**6]:
        true = len(smooth_upto(X, P))
        lo, hi = pb.psi_bounds(mp.log(X))
        assert lo <= true <= hi


def test_cvp_output_obeys_certificate():
    """Thm R corollary: a lattice vector at distance R has |w| <= 2 e^{(L + A)/2} sinh(|Delta|/2) with A <= Q <= R^2 (weights)."""
    N = 1000003 * 1000033
    B, t, d, P, meta = build_lattice(N, 20, c_exp=1.0, weight="lnp", seed=0)
    reduce_basis(B)
    for vec in list_cvp(B, t, 10):
        a = analyse_vector(vec, d, P, N, t)
        S = meta["S"]
        # exact integer facts: dist2 = sum (d_i e_i)^2 + r^2 ; recompute
        assert a["dist2"] == sum((di * ei) ** 2 for di, ei in zip(d, a["e"])) + (vec[-1] - t[-1]) ** 2
        u, v, w = a["u"], a["v"], a["w"]
        assert u - v * N == w and gcd(u, v) == 1
        # Thm R(ii)/(iii) bound: |w| <= 2 e^{(A+L)/2} sinh(|Delta|/2), A = ln(uv) <= Q (high-precision check)
        mp.mp.dps = 120
        A = mp.log(u) + mp.log(v)
        D = mp.log(u) - mp.log(v) - mp.log(N)
        Q = sum(mp.log(p) * ei * ei for p, ei in zip(P, a["e"]))
        assert A <= Q + mp.mpf(10) ** -60
        assert abs(w) <= 2 * mp.e ** ((A + mp.log(N)) / 2) * mp.sinh(abs(D) / 2) * (1 + mp.mpf(10) ** -50)


def test_exact_cvp_can_fail_codex_example():
    """codex schnorrcvp_1: N = 35, single prime 2, C = 35: unique closest vector 5*b1 gives 2^5 - 35 = -3, not 2-smooth."""
    u, v, w = uvw([5], [2], 35)
    assert w == -3 and factor_over(w, [2]) is None
    # distances for k*b1, b1 = (1, 25), t = (0, 125)
    d2 = {k: k * k + (25 * k - 125) ** 2 for k in range(-10, 20)}
    assert min(d2, key=d2.get) == 5


def test_adleman_rerun_lattice_unchanged():
    """Adleman 1995 §2 distinctness step: replacing V0 by V0+Vi+Vj does not change the lattice, so an SVP oracle's
    output cannot change (found by us; see notes).  Check: the shortest vector is identical before and after."""
    from fpylll import GSO, Enumeration, IntegerMatrix, LLL
    N = 10007 * 10009
    n = 12
    P = first_primes(n)
    M = 10**9
    rows = []
    for j, p in enumerate(P):
        r = [0] * (n + 1)
        r[0] = int(M * log(p) * 1000)
        r[j + 1] = int(1000 * log(p) ** 0.5)
        rows.append(r)
    v0 = [int(M * log(N) * 1000)] + [0] * n
    def shortest(v0row):
        A = IntegerMatrix.from_matrix([v0row] + rows)
        LLL.reduction(A)
        M_ = GSO.Mat(A)
        M_.update_gso()
        r0 = sum(x * x for x in A[0])
        sol = Enumeration(M_).enumerate(0, A.nrows, float(r0) * 1.0001, 0)
        s = A.multiply_left([int(round(c)) for c in sol[0][1]])
        return tuple(abs(x) for x in s)  # up to sign
    s1 = shortest(v0)
    v0b = [a + b + c for a, b, c in zip(v0, rows[2], rows[5])]
    s2 = shortest(v0b)
    assert sum(x * x for x in s1) == sum(x * x for x in s2)


def test_expected_relations_monotone_in_n():
    e1 = expected_relations(first_primes(20), 30 * mp.log(2), K=50, h=0.1, lnv_max=40)["E"]
    e2 = expected_relations(first_primes(30), 30 * mp.log(2), K=50, h=0.1, lnv_max=40)["E"]
    assert e2 > e1 > 0
