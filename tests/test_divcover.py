import itertools
import math
import random

import pytest
import sympy

from src.divcover.core import (GAP, brute_min_height, c_min_height, covers, uncovered,
                               residue_cover_prob, residue_cover_prob_np, shapes_of_rank,
                               incidence_params, hs_r_bound_holds, random_model_log_height)


def test_residues_match_elements():
    rng = random.Random(1)
    for _ in range(200):
        r = rng.randint(1, 3)
        Q = GAP(rng.randint(-50, 50), tuple(rng.randint(-30, 30) for _ in range(r)),
                tuple(rng.randint(1, 5) for _ in range(r)))
        d = rng.randint(1, 40)
        assert Q.residues(d) == {e % d for e in Q.elements()}


def test_known_prime_cover():
    # {2926, 2958, 2990}: 2*7*11*19, 2*3*17*29, 2*5*13*23
    Q = GAP(2926, (32,), (3,))
    assert covers(Q, 30, "prime")
    assert not covers(Q, 30, "full")
    assert uncovered(Q, 31, "prime") == [31]


@pytest.mark.parametrize("n", [8, 12, 16])
@pytest.mark.parametrize("mode", ["prime", "full"])
@pytest.mark.parametrize("shape", [(2,), (3,), (4,), (2, 2), (3, 2)])
def test_c_search_matches_bruteforce(n, mode, shape):
    Hmax = 300 if len(shape) == 1 else 120
    bf = brute_min_height(n, shape, Hmax, mode)[0]
    cr = c_min_height(n, shape, Hmax, mode)
    assert bf == cr["best"]
    if cr["best"] is not None:
        Q = GAP(cr["b"], cr["c"], shape)
        assert Q.height() == cr["best"] and covers(Q, n, mode)


def test_c_search_rank3_witness():
    res = c_min_height(14, (2, 2, 2), 200, "prime")
    Q = GAP(res["b"], res["c"], (2, 2, 2))
    assert covers(Q, 14, "prime") and Q.height() == res["best"]
    # exactness: brute force over rank-3 shape agrees
    assert brute_min_height(14, (2, 2, 2), res["best"], "prime")[0] == res["best"]


def test_residue_prob_np_matches():
    for p in [2, 3, 5, 7, 11, 13]:
        for sh in [(3,), (2, 2), (4, 2), (2, 2, 2)]:
            assert abs(residue_cover_prob(p, sh) - residue_cover_prob_np(p, sh)) < 1e-12


def test_residue_prob_rank1_closed_form():
    # rank 1, L < p: P_p = ((p-1) * L + 1) / p^2   (c != 0: L residues; c == 0: 1 residue)
    for p in [7, 11, 13]:
        for L in [2, 3, 5]:
            assert abs(residue_cover_prob(p, (L,)) - ((p - 1) * L + 1) / p ** 2) < 1e-12


def test_shapes():
    assert shapes_of_rank(1, 12) == [(12,)]
    assert set(shapes_of_rank(2, 12)) == {(6, 2), (4, 3)}
    assert set(shapes_of_rank(3, 12)) == {(3, 2, 2)}


def test_random_model_monotone_in_rank():
    h1 = random_model_log_height(40, shapes_of_rank(1, 8))
    h2 = random_model_log_height(40, shapes_of_rank(2, 8))
    assert h2 < h1


def test_hs_r_lemma_random_block_systems():
    """Property test of Prop. HS-r (v <= 1 + Delta*(mu*Delta - 1)) on random element sets."""
    rng = random.Random(7)
    P = list(sympy.primerange(11, 60))
    checked = 0
    for _ in range(3000):
        k = rng.randint(3, 12)
        els = []
        for _ in range(k):
            sub = rng.sample(P, rng.randint(1, len(P) - 1))
            els.append(math.prod(sub) * rng.randint(1, 5))
        v, D, mu, ok, proper = incidence_params(els, P)
        if ok and proper:
            checked += 1
        assert hs_r_bound_holds(els, P)
    assert checked > 50


def test_hs_r_reduces_to_he_sahai_for_ap():
    # rank-1 AP with pq > length and p,q not dividing c: mu <= 1
    P = [11, 13, 17]
    for b in range(1, 400):
        Q = GAP(b, (1,), (60,))
        v, D, mu, ok, proper = incidence_params(list(Q.elements()), P)
        assert mu <= 1
