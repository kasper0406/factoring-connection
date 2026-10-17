from src.seed_checks import check_seed_lemmas, check_extended, kummer_support, carries
from math import comb


def test_brief_smoke():
    assert check_seed_lemmas(200) == 153


def test_extended_1000():
    st = check_extended(1000)
    assert st["composites"] == 831
    assert st["scalar_zero_composites"] == [561]


def test_kummer_examples():
    assert carries(1, 1, 2) == 1
    assert kummer_support(15) == {3, 5, 6, 9, 10, 12}
    for N in (45, 49, 105):
        assert kummer_support(N) == {k for k in range(1, N) if comb(N, k) % N}
