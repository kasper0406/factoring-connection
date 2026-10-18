"""Tests for src/sieve2d (two-variable coverage harness)."""
from itertools import product
from math import gcd, isqrt

import numpy as np
import pytest
import sympy

from src.dcltheorem.params import Design, weight_set
from src.sieve2d.coverage import (band_tables, covered_by_count, degenerate_mask,
                                  exhaustive_integer_pairs, exhaustive_prime_pairs,
                                  model_integer, model_prime, multiplicity, sn_bounds)
from src.sieve2d.sampling import in_region, sample_integer_pairs, sample_prime_pairs


def brute_M(x, y, band, V):
    return sum(1 for l in band for (A, B) in V if (B * x - A * y) % l == 0)


@pytest.mark.parametrize("H,band", [(3, [7, 11, 19]), (6, [5, 7, 11, 13]), (4, [43, 47, 59])])
def test_multiplicity_matches_brute(H, band):
    V = weight_set(H)
    tabs = band_tables(band, V)
    rng = np.random.default_rng(1)
    x = rng.integers(1, 5000, size=3000).astype(np.int64)
    y = rng.integers(1, 5000, size=3000).astype(np.int64)
    # include points divisible by band primes
    x[:50] = band[0] * rng.integers(1, 100, size=50)
    y[:25] = band[0] * rng.integers(1, 100, size=25)
    M = multiplicity(x, y, tabs)
    Nc = covered_by_count(x, y, tabs)
    for i in range(3000):
        xi, yi = int(x[i]), int(y[i])
        assert M[i] == brute_M(xi, yi, band, V)
        assert Nc[i] == sum(1 for l in band if any((B * xi - A * yi) % l == 0 for A, B in V))


@pytest.mark.parametrize("H,band", [(3, [7, 11]), (5, [13, 17, 19]), (7, [5, 11])])
def test_models_match_residue_enumeration(H, band):
    V = weight_set(H)
    tabs = band_tables(band, V)
    mp = model_prime(tabs)
    mi = model_integer(tabs)
    # brute: enumerate residues mod prod(band)
    Lp = 1
    for l in band:
        Lp *= l
    vals_p, vals_i = [], []
    for a in range(Lp):
        for b in range(Lp):
            m = brute_M(a, b, band, V)
            vals_i.append(m)
            if all(a % l and b % l for l in band):
                vals_p.append(m)
    for vals, mod in ((vals_p, mp), (vals_i, mi)):
        v = np.array(vals, dtype=float)
        assert abs((v == 0).mean() - mod["f_model"]) < 1e-12
        assert abs(v.mean() - mod["mu"]) < 1e-9
        assert abs(v.var() - mod["var"]) < 1e-9
        assert abs(((v - v.mean()) ** 4).mean() - mod["m4"]) < 1e-6
        assert abs(((v - v.mean()) ** 3).mean() - mod["m3"]) < 1e-6


@pytest.mark.parametrize("n", [16, 19, 20])
def test_exhaustive_pairs(n):
    p, q = exhaustive_prime_pairs(n)
    plo, phi, qlo, qhi = sn_bounds(n)
    brute = set()
    for x in sympy.primerange(2, phi + qhi):
        for y in sympy.primerange(x + 1, 2 * x + 1):
            if (1 << (n - 1)) <= x * y < (1 << n):
                brute.add((x, y))
    assert set(zip(p.tolist(), q.tolist())) == brute
    x, y = exhaustive_integer_pairs(n)
    bi = {(a, b) for a in range(1, phi + 2) for b in range(a + 1, 2 * a + 1)
          if (1 << (n - 1)) <= a * b < (1 << n)}
    assert set(zip(x.tolist(), y.tolist())) == bi


def test_samplers_in_region_and_prime():
    p, q = sample_prime_pairs(48, 200, seed=3)
    assert in_region(48, p, q).all()
    assert all(sympy.isprime(int(a)) and sympy.isprime(int(b)) for a, b in zip(p, q))
    x, y = sample_integer_pairs(70, 200, seed=3)
    assert in_region(70, x, y).all()
    p2, q2 = sample_prime_pairs(48, 200, seed=3)
    assert (p2 == p).all() and (q2 == q).all()     # seeded


def test_degenerate_mask():
    V = weight_set(4)
    x = np.array([6, 6, 7], dtype=np.int64)
    y = np.array([9, 10, 11], dtype=np.int64)       # 6:9 = 2:3 is degenerate
    assert degenerate_mask(x, y, V).tolist() == [True, False, False]


def test_prime_pairs_never_covered_by_both_zero():
    # for prime pairs above the band no l divides p or q: M_l <= 1 when classes distinct
    d = Design(n=48, H=4, Y=200)
    tabs = band_tables(d.bands[0], d.V)
    p, q = sample_prime_pairs(48, 500, seed=5)
    assert covered_by_count(p, q, tabs).tolist() == multiplicity(p, q, tabs).tolist()


def test_cell_record_qmodel_matches_brute():
    from src.sieve2d.experiments import cell_record, eval_cells
    band, V = [7, 11], weight_set(3)
    # one full period of residues: empirical == model exactly
    a, b = np.meshgrid(np.arange(77, dtype=np.int64), np.arange(77, dtype=np.int64))
    x, y = a.ravel() + 77, b.ravel() + 77
    M = eval_cells(x, y, band, [V])[0]
    rec = cell_record(0, "int", "t", 0, 3, 5, band, V, M)
    assert abs(rec["Qstat"] - rec["Qmodel"]) < 1e-9
    assert rec["uncovered"] == round(rec["E_model"])
