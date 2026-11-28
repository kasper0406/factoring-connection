"""Check the exceptional dyadic lifts and a real cyclotomic splitting example."""

from math import prod
from sympy import Poly, symbols


def subgroup(a, m):
    result = {1}
    value = a % m
    while value != 1:
        result.add(value)
        value = value*a % m
    return result


def check_lifts():
    tested = 0
    for k in range(4, 12):
        m = 1 << k
        eps = 1+m//2
        assert pow(5, m//8, m) == eps
        exceptional = {1, m-1, eps, -eps % m}
        for a in range(1, m, 2):
            other = (a+m//2) % m
            assert a*eps % m == other
            assert (subgroup(a, m) != subgroup(other, m)) == (a in exceptional)
            for h in (1, 3, 5):
                assert bool((h*a) & (m//2)) != bool((h*other) & (m//2))
            tested += 1
    return tested


def check_fields():
    x = symbols("x")
    real_field = x**8-8*x**6+20*x**4-16*x*x+2
    pairs = [(67, 131), (179, 211)]
    assert [prod(pair) % 32 for pair in pairs] == [9, 9]
    for i, pair in enumerate(pairs):
        for h in (1, 3, 5):
            assert bool((h*pair[0]) & (h*pair[1]) & 16) == bool(i)
    patterns = {}
    for p in (67, 131, 179, 211, 97, 113):
        factors = Poly(real_field, x, modulus=p).factor_list()[1]
        patterns[p] = sorted(f.degree() for f, e in factors for _ in range(e))
    assert all(patterns[p] == [8] for pair in pairs for p in pair)
    assert patterns[97] == [1]*8 and patterns[113] == [2]*4
    return patterns


if __name__ == "__main__":
    print("Dyadic lift classifications:", check_lifts())
    print("Real cyclotomic splitting degrees:", check_fields())
