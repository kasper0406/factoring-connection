import itertools
import random

from src.ptecons import acs, lemmas
from src.ptecons.polys import (all_elem_sym, centered, dickson, evaluate, from_roots, int_roots_split,
                               poly_mod, resultant)


def test_gap_r_exhaustive():
    for M in (3, 4, 5, 6):
        tot, viol = lemmas.gap_r_exhaustive(M, 4)
        assert tot > 0 and viol == 0


def test_lin_r_real_directions():
    # several real directions: e_2 or e_3 not affine
    for betas in ([[1], [1], [0], [0]], [[2], [2], [-1], [0]], [[1, 0], [0, 1], [1, 1], [0, 0]]):
        deg = lemmas.lin_family_affine_degree([0, 1, 2, 3], betas)
        assert deg[1] >= 2 or deg[2] >= 3
    # one moving root: every coefficient affine
    deg = lemmas.lin_family_affine_degree([0, 1, 2, 3], [[1], [0], [0], [0]])
    assert max(deg) <= 1


def test_lin_escape_needs_roots_of_unity():
    import sympy as sp
    T = sp.Symbol("T")
    assert sp.expand(sp.prod([1 + u * T for u in (1, sp.I, -1, -sp.I)])) == 1 - T ** 4


def test_ledger_and_units():
    rng = random.Random(1)
    F = [3, -1, 0, 2, 0, 1]
    for _ in range(20):
        roots = [rng.randint(-9, 9) for _ in range(rng.randint(1, 9))]
        assert lemmas.ledger_holds(F, roots)


def test_cunningham_bound():
    rng = random.Random(2)
    for _ in range(200):
        n = rng.randint(3, 9)
        roots = [rng.choice([-1, 0, 1, 0, rng.randint(-6, 6)]) for _ in range(rng.randint(1, 12))]
        rhs, lhs = lemmas.cunningham_bound(n, roots)
        assert rhs >= lhs
    # codex ptecons_2 counterexample to the uncorrected odd-n form: n = 3, roots 1, (-2)^84
    rhs, lhs = lemmas.cunningham_bound(3, [1] + [-2] * 84)
    assert rhs >= lhs


def test_galois_orbit_collapse():
    assert lemmas.galois_orbit_collapse(5, 7, 2)
    assert lemmas.galois_orbit_collapse(6, 3, 1)


def test_dickson_integer_split_classification():
    for n in range(3, 13):
        exp = [-2, -1, 1, 2] if n in (3, 6) else ([-2, 0, 2] if n == 4 else [])
        assert lemmas.dickson_integer_split(n, 10) == exp
    # the codex example
    c = list(dickson(6))
    c[0] -= 2
    assert int_roots_split(c) == [-2, -1, -1, 1, 1, 2]
    assert poly_mod(c, [-7, -1, 9, 0, -6, 0, 1]) == [3, 1]


def test_root_support_rank_disjoint_orbits():
    fam = [tuple(range(10 * i, 10 * i + 5)) for i in range(7)]
    rows, cols = lemmas.root_support_rank(fam)
    assert rows < cols


def test_quart_g_matches_bruteforce():
    N = 1000003 * 999983
    new = sorted(r[0] for r in acs.quart_g(N, 10, 2, 40, 10 ** 9))
    out = []
    for k in (1, 2):
        m = acs.iroot(k * N, 4) + 1
        R0 = m ** 4 - k * N
        E1 = (R0 + m ** 3 // 2) // m ** 3
        R1 = R0 - E1 * m ** 3
        E2 = -((R1 + m ** 2 // 2) // m ** 2)
        for c in itertools.combinations_with_replacement(range(-10, 11), 4):
            if sum(c) == E1 and abs(all_elem_sym(c)[2] - E2) <= 40:
                out.append(tuple(m - a for a in c))
    assert new == sorted(out)


def test_acs_relations_are_correct():
    N = 1000003 * 999983
    for gen in (acs.cos(N, 5, 50), acs.quart_s(N, 400, 30, 20, 50), acs.quart_rand(N, 50, 1, 50)):
        for pieces, u, v in gen:
            prod = 1
            for p in pieces:
                prod *= p
            assert prod == u and (u - v) % N == 0 and abs(v) <= N // 2


def test_operational_no_leakage():
    import inspect
    src = inspect.getsource(acs)
    for bad in ("factorint", "factor(", "isprime", "labels", "phi"):
        assert bad not in src


def test_star_gaussian():
    rng = random.Random(3)
    units = [1, -1, 1j, -1j]
    for _ in range(300):
        half = []
        for _ in range(rng.randint(1, 6)):
            c = rng.choice(units) if rng.random() < 0.7 else complex(rng.randint(-3, 3), rng.randint(-3, 3))
            half.append(c)
        cs = half + [z.conjugate() for z in half]
        Mp, s, bound = lemmas.star_gaussian_check(cs)
        assert Mp <= bound
    # mu_4 orbit: prod (Y - u) = Y^4 - 1, s = 0
    assert lemmas.star_gaussian_check([1, -1, 1j, -1j])[:2] == (4, 0)
