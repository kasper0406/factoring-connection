"""Tests for src/lem85: exact checks of the elementary steps (E2)-(E5) of notes/lem85_results.md."""
import math

import pytest

from src.lem85 import identities as I


@pytest.mark.parametrize("x,y,ell,s", [(400, 7, 11, 6), (600, 5, 13, 10), (900, 7, 17, 12), (500, 3, 7, 9)])
def test_lift_average_identity(x, y, ell, s):
    for r in range(1, s + 1):
        if math.gcd(r, s) == 1:
            avg, rhs = I.lift_average(x, y, ell, s, r)
            assert avg == rhs


@pytest.mark.parametrize("x,y,r,q", [(500, 7, 3, 30), (700, 11, 5, 42), (800, 7, 7, 60), (600, 13, 15, 2 * 7 * 11 * 17)])
def test_mobius_identity(x, y, r, q):
    lhs, rhs = I.mobius_identity(x, y, r, q)
    assert lhs == rhs


@pytest.mark.parametrize("y,M1,z", [(5, 6, 400), (7, 10, 600), (3, 4, 300)])
def test_mk_decomposition_unique(y, M1, z):
    for n in I.smooth_list(z, y):
        decs = I.mk_decompositions(n, y, M1)
        if n > M1:
            assert len(decs) == 1, (n, decs)
        else:
            assert decs == []


def test_sum_inv_phi_bound():
    for X in [1, 2, 5, 10, 50, 200, 600]:
        assert I.sum_inv_phi(X) <= 2 * (1 + math.log(X))
