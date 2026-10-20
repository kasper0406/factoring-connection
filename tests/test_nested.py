"""Tests for agent nested (src/nested)."""
import random
from math import gcd, isqrt

from src.common.instances import balanced_semiprime
from src.nested import toy_e1 as T
from src.nested import barrier as BR
from src.thmA import numfield as nf


def _field(seed=5, bits=36):
    N = balanced_semiprime(bits, seed).N
    return N, T.sample_field(N, random.Random(seed))


def test_field_has_root_mod_N():
    N, A = _field()
    assert nf.hom_eval(A["f"], A["g0"], A["g1"]) == N * A["w"]


def test_box_points_complete_and_in_lattice():
    q, lam, V = 101, 37, 60
    pts = set(T.box_points(q, lam, V))
    brute = {(a, b) for a in range(-V, V + 1) for b in range(1, V + 1)
             if (a - lam * b) % q == 0 and gcd(a, b) == 1}
    assert pts == brute


def test_targets_divide_values():
    N, A = _field()
    ctx = T.Ctx(256, 4096)
    U = T.Uniform(A, ctx)
    rng = random.Random(1)
    for t in [rng.choice(U.alg) for _ in range(5)]:
        _, q, r = t
        for (a, b) in T.box_points(q, r, isqrt(q * 16))[:20]:
            assert T.F_val(A, a, b) % q == 0
    for q in U.rat[:5]:
        _, _, lam = T.rational_target(A, q)
        for (a, b) in T.box_points(q, lam, isqrt(q * 16))[:20]:
            assert T.G_val(A, a, b) % q == 0


def test_parent_children_are_prime_factors():
    N, A = _field(seed=7, bits=40)
    ctx = T.Ctx(512, 8192)
    found = 0
    for ell in range(1 << 18, (1 << 18) + 4000):
        if not T.arith.is_prime(ell):
            continue
        node = T.parent_node(A, ell, 8, ctx)
        if node is None:
            continue
        found += 1
        a, b = node["a"], node["b"]
        G, F = T.G_val(A, a, b), T.F_val(A, a, b)
        assert G % ell == 0
        for kind, p, r in node["children"]:
            assert ctx.B < p <= ctx.Y
            if kind == "r":
                assert (G // ell) % p == 0
            else:
                assert F % p == 0 and (a - r * b) % p == 0
        if found >= 3:
            break
    assert found >= 1


def test_d1_matches_bruteforce():
    N, A = _field(seed=3, bits=36)
    ctx = T.Ctx(256, 4096)
    U = T.Uniform(A, ctx)
    rng = random.Random(2)
    for t in [rng.choice(U.alg) for _ in range(4)] + [T.rational_target(A, rng.choice(U.rat)) for _ in range(4)]:
        kind, q, lam = t
        V = isqrt(q * 64)
        ok = False
        for (a, b) in T.box_points(q, lam, V):
            G, F = T.G_val(A, a, b), T.F_val(A, a, b)
            if G == 0 or F == 0:
                continue
            x, y = (abs(G) // q, abs(F)) if kind == "r" else (abs(G), abs(F) // q)
            if BR.is_smooth(x, ctx.B) and BR.is_smooth(y, ctx.B):
                ok = True
                break
        assert T.d1(A, t, 8, ctx, {}) == ok


def test_gamma_tail_matches_poisson():
    import math
    from fractions import Fraction as Fr
    # Pr[Gamma(1,1) >= V] = e^{-V};  Pr[Gamma(2,1) >= V] = e^{-V}(1+V)
    assert BR.gamma_upper_tail(1, Fr(3)) == 1
    assert BR.gamma_upper_tail(2, Fr(3)) == 4
    assert abs(BR.af_weight(2, Fr(100), Fr(1)) - math.exp(-math.log(100)) * (1 + math.log(100))) < 1e-12


def test_overdetermination_counts():
    from fractions import Fraction as Fr
    # lam < 1/2: one level of children already exceeds d + 3 parameters for d >= 5
    assert BR.overdetermined(9, Fr(1, 4), 1)
    # lam in (1/2, 1): one level fits, two levels do not
    assert not BR.overdetermined(9, Fr(3, 4), 1)
    assert BR.overdetermined(9, Fr(3, 4), 2)


def test_backtracking_model():
    from fractions import Fraction as Fr
    # constant p per level: overhead m*levels*log(1/p); polylog p: extra factor log t
    assert BR.backtrack_log_cost(Fr(10), 4, [Fr(1)] * 4) == 40
    assert BR.backtrack_log_cost(Fr(10), 4, [Fr(5)] * 4) == 200


def test_nr_constant_example():
    c = BR.nr_constant(4, 1, 0.5, 1)
    assert 24 < c < 25
    assert BR.ad_heuristic_ok(4, 0.5, 1, 0.5, 1, 1)
    assert not BR.ad_heuristic_ok(4, 0.3, 1, 0.5, 1, 1)
