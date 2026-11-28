"""Compare bit-cube products, binomial divisibility and direct balanced factoring."""

from itertools import combinations
from math import isqrt

from sympy import primerange

if __package__:
    from .check_cofactor_zero_queries import direct_zero
    from .cofactor_zero_queries import semiprime_from_pairs
    from .submask_batch import _make_cube, submask_window_zero, submask_zero
else:
    from check_cofactor_zero_queries import direct_zero
    from cofactor_zero_queries import semiprime_from_pairs
    from submask_batch import _make_cube, submask_window_zero, submask_zero


def check_small():
    queries = 0
    for p, q in combinations(primerange(3, 51), 2):
        if q >= 2*p:
            continue
        n = p*q
        lo, hi = isqrt((n-1)//2)+1, isqrt(2*n)
        for mask in range(0, 1 << hi.bit_length(), 2):
            work = {}
            result = submask_zero(n, mask, stats=work)
            assert result == direct_zero(n, mask) == bool((p & mask) and (q & mask))
            assert submask_window_zero(n, mask) == result
            odds = range(lo | 1, hi+1, 2)
            assert work["avoiding"] == sum(not (value & mask) for value in odds)
            assert work["hitting"] == len(odds)-work["avoiding"]
            assert work["divisor"] in (None, p, q)
            queries += 1
        assert semiprime_from_pairs(n, submask_zero) in (p, q)
    print(f"{queries} mask comparisons for both evaluators against direct products")


def check_regimes():
    stats = {}
    assert not submask_zero(17*19, 8, stats=stats)
    assert stats["saturated_rows"] == 1 and stats["divisor"] in (17, 19)
    assert not submask_zero(77, 1 << 5, stats=stats)
    assert stats["side"] == "hitting" and stats["block_values"] == stats["edge_gcds"] == 0
    assert submask_zero(77, 0b1110, stats=stats)
    assert stats["side"] == "avoiding" and stats["block_values"] == stats["edge_gcds"] == 0

    n, mask = 1000003*1000033, 6
    assert not submask_zero(n, mask, stats=stats)
    assert stats["edge_gcds"] > 0 and stats["block_values"] > 0
    print(f"40-bit two-ratio query: {stats['avoiding']} candidate integers, "
          f"degree {stats['polynomial_degree']}, {stats['block_values']} blocks, "
          f"{stats['edge_gcds']} boundary gcds")


def check_benchmarks():
    # Positions for the negative control were chosen using the special prime.
    # The evaluator receives only the selected positions and the product.
    p, q = 2**127-1, 162*2**120+1
    n = p*q
    all_positions = (1 << (n.bit_length()-1))-2
    free = (1 << 121) | (1 << 125) | (1 << 127) | sum(1 << i for i in range(3, 35, 2))
    for omitted, expected in ((free, False), (free ^ (1 << 125) ^ (1 << 124), True)):
        stats = {}
        mask = all_positions ^ omitted
        assert submask_zero(n, mask, stats=stats) == expected
        assert expected == bool((p & mask) and (q & mask))
        assert stats["free"] == 19 and stats["avoiding"] == 262144
        assert stats["side"] == "avoiding" and stats["polynomial_degree"] == 512
        assert stats["block_values"] == 512 and stats["edge_gcds"] == 0
        assert stats["divisor"] == (None if expected else q)
        print(f"255-bit query {expected}: {stats['avoiding']} candidates, "
              f"degree {stats['polynomial_degree']}, {stats['block_values']} blocks")


def check_middle_windows():
    # An incidental hit in the relaxed coefficient interval must yield the actual bit.
    stats = {}
    assert submask_window_zero(47111*70667, 20676, stats=stats)
    assert stats["divisor"] == 70667 and stats["width"] == 47

    p, q = 2**127-1, 162*2**120+1
    n = p*q
    lower, upper = isqrt((n-1)//2)+1, isqrt(2*n)
    all_positions = (1 << (n.bit_length()-1))-2
    middle = ((1 << 40)-1) << 40
    for high, expected in ((125, False), (124, True)):
        omitted = middle | (1 << high) | (1 << 121) | (1 << 127)
        mask = all_positions ^ omitted
        avoiding = _make_cube(1, omitted, lower, upper).size
        assert avoiding == 4*(1 << 40)
        stats = {}
        assert submask_window_zero(n, mask, stats=stats) == expected
        assert expected == bool((p & mask) and (q & mask))
        assert stats["free"] == 43 and stats["grouped"] == 40
        assert stats["stride"] == 1 << 40 and stats["width"] == (1 << 40)-1
        assert stats["outside"] == stats["intervals"] == 4
        assert stats["maximum_dimension"] == 3
        assert stats["divisor"] == (None if expected else q)
        print(f"255-bit middle window {expected}: {avoiding} candidates, "
              f"{stats['intervals']} certified dimension-three affine calls")


if __name__ == "__main__":
    check_small()
    check_regimes()
    check_benchmarks()
    check_middle_windows()
