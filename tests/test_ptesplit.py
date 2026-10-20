import math

import numpy as np
from flint import fmpz_poly

from src.ptesplit import census as C
from src.ptesplit import nccheck as NC


def test_multisets_complete_and_sorted():
    for M, R in [(3, 5), (4, 7), (5, 3)]:
        a = C.multisets(M, R)
        assert a.shape == (C.count_multisets(M, R), M)
        assert np.all(np.diff(a.astype(int), axis=1) >= 0)
        assert len({tuple(r) for r in a}) == a.shape[0]


def test_census_matches_bruteforce():
    for M, R in [(4, 12), (5, 8), (6, 6)]:
        rows = {r["k"]: r for r in C.census_cell(M, R)}
        for k, row in rows.items():
            d = C.brute_fibres(M, R, M - k)
            sizes = [len(v) for v in d.values()]
            assert row["n_fibres"] == len(d)
            assert row["max_fibre"] == max(sizes)
            assert row["pairs"] == sum(s * (s - 1) // 2 for s in sizes)


def test_injection_bound_and_kleiman_small():
    rows = C.census_cell(4, 30)
    for r in rows:
        assert r["max_fibre"] <= r["injection_bound"]
        if r["k"] == 1:
            assert r["kleiman_violations"] == 0


def test_known_ideal_pte_size5():
    # [0,4,8,16,17] =_4 [1,2,10,14,18]
    A, Bm = [0, 4, 8, 16, 17], [1, 2, 10, 14, 18]
    assert all(sum(a ** i for a in A) == sum(b ** i for b in Bm) for i in range(1, 5))
    rows = {r["k"]: r for r in C.census_cell(5, 18)}
    assert rows[1]["n_nontrivial"] >= 1
    c = abs(math.prod(A) - math.prod(Bm))
    assert c % math.factorial(4) == 0


def test_nc_ledger_identity():
    F = fmpz_poly([-2, 0, 0, 1])  # X^3 - 2, irreducible
    row = NC.nc_row(F, 4, 4, 40)
    assert row["irreducible"] and row["violations"] == 0
    assert row["n_members"] >= row["n_distinct_r"] >= 1
