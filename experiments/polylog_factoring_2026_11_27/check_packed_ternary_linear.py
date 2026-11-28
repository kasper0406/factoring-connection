"""Compare signed CRT with unconstrained frame enumeration and native sphere counts."""

from collections import Counter
from itertools import combinations
from math import gcd

from sympy import primerange

if __package__:
    from .check_mixed_ternary_embeddings import _all_frames
    from .packed_ternary_decision import _packing_parameters, brute_sphere_count
    from .packed_ternary_linear import _centered_crt, recover_packed_linear_factor, select_linear_x
    from .short_jacobi_grams import first_vector_orbits, select_x
else:
    from check_mixed_ternary_embeddings import _all_frames
    from packed_ternary_decision import _packing_parameters, brute_sphere_count
    from packed_ternary_linear import _centered_crt, recover_packed_linear_factor, select_linear_x
    from short_jacobi_grams import first_vector_orbits, select_x


def check_selectors():
    cases = 0
    for p, q in combinations(primerange(3, 150), 2):
        n = p*q
        for limit in (7, 12):
            old = select_x(n, limit)
            x, factor = select_linear_x(n, limit)
            assert old == (None, None) or x or factor
            if factor:
                assert factor in (p, q)
            if x:
                assert 2 <= x <= limit and gcd(n, x) == 1
            cases += 1
    assert select_linear_x(15, 7) == (None, 5)
    assert select_linear_x(1121, 7) == (None, None)
    assert select_linear_x(1121, 9) == (9, None)
    return cases


def check_case(n, x):
    scale, norm = _packing_parameters(n, x)
    marked = set(first_vector_orbits(x))
    partners = {u: [] for u in marked}
    for u, v in _all_frames(n, x):
        if u in marked:
            partners[u].append(v)
    packed = [tuple(a - scale*b for a, b in zip(u, v))
              for u, vectors in partners.items() for v in vectors]
    assert all(sum(c*c for c in point) == norm for point in packed)

    counts = {}
    records = []
    for singleton_only in (False, True):
        queries = multiple = 0

        def read(m, modulus, residue):
            nonlocal queries, multiple
            assert m == norm
            if modulus not in counts:
                counts[modulus] = Counter(tuple(c % modulus for c in point) for point in packed)
            count = brute_sphere_count(m, modulus, residue)
            assert count == counts[modulus][residue]
            queries += 1
            multiple += count > 1
            return count % 2 if singleton_only else count > 0

        result = recover_packed_linear_factor(n, x, read, singleton_only=singleton_only)
        u, j = result.first_vector, result.index
        expected = {v[j] for v in partners[u] if gcd(*v) < n}
        assert expected == {result.coordinate}
        assert -n < result.coordinate < 0 and gcd(n, result.coordinate) == result.factor
        assert result.queries == queries and result.crt_modulus > 2*n
        records.append((queries, multiple, result))
    assert records[0][2].factor == records[1][2].factor
    if (n, x) == (15, 4):
        assert partners[(0, 1, 4)] == [(0, 0, -15)] and len(partners[(2, 2, 3)]) == 2
    if (n, x) == (15, 11):
        assert not partners[(4, 5, 9)]
    if (n, x) == (39, 7):
        assert not partners[(0, 5, 5)]
    if (n, x) == (395, 11):
        # The genuine signed coordinate and the public axis coordinate agree mod 11.
        assert records[0][2].coordinate == -55 and records[0][2].crt_modulus % 11 == 0
    if (n, x) == (395, 6):
        # Two distinct partners coincide mod 11. Parity removes exactly this prime.
        assert records[1][1] > 0 and records[1][2].crt_modulus % 11 != 0
        assert records[0][2].crt_modulus % 11 == 0
    return records


def main():
    primes = (11, 13)
    for t in range(-71, 72):
        assert _centered_crt([(p, t % p) for p in primes]) == (t, 143)
    print(f"{check_selectors()} public selector comparisons; 143 centered CRT lifts")
    for n, x in ((15, 4), (15, 11), (39, 7), (65, 4), (77, 2), (395, 6), (395, 11)):
        exact, parity = check_case(n, x)
        found = exact[2]
        print(f"n={n} x={x} u={found.first_vector} t={found.coordinate} gcd={found.factor} "
              f"queries={exact[0]}/{parity[0]} multiple={exact[1]}/{parity[1]}")


if __name__ == "__main__":
    main()
