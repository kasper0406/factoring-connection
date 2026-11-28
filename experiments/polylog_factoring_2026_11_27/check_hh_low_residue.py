"""Independent checks of the order guarantee, normalization, collisions and crossover."""

from itertools import combinations
from math import gcd

from sympy import n_order, primerange

if __package__:
    from .hh_low_residue import (_bounded_order, hh_collision, hh_giants, hh_low_factor,
                                 hh_order_base, hh_parameters)
else:
    from hh_low_residue import (_bounded_order, hh_collision, hh_giants, hh_low_factor,
                                hh_order_base, hh_parameters)


def check_orders():
    counts = 0
    for n in (15, 21, 35, 55, 77, 91, 143, 323, 2021):
        for bound in (1, 2, 7, 16, 41, 100):
            d, beta = hh_order_base(n, bound)
            if d:
                assert 1 < d < n and n % d == 0 and beta is None
            else:
                assert n_order(beta, n) > bound
            for base in range(1, min(n, 80)):
                if gcd(n, base) == 1:
                    order = int(n_order(base, n))
                    assert _bounded_order(n, base, bound) == (order if order <= bound else None)
                    counts += 1
    print(f"{counts} bounded multiplicative orders matched independent orders")


def check_small():
    counts = 0
    for p, q in combinations(primerange(3, 160), 2):
        if q >= 2*p:
            continue
        n = p*q
        for m in (2, 4, 8, 16, 32, 64, 128):
            assert hh_low_factor(n, m, p % m) in (p,q)
            assert hh_low_factor(n, m, (p+2) % m) in (None,p,q)
            counts += 1
    print(f"{counts} small promised states split, with wrong-class controls checked")


def check_collisions():
    cases = (
        (1777, 1831, 8, 8, 53974, "collision"),
        (1187, 1913, 4, 2, 243246, "baby_gcd"),
        (1123, 1307, 8, 3, 800746, "integer_sum"),
        (937, 1249, 8, 3, 83077, "short_order"),
        (1279, 2273, 8, 4, 2775119, "two_part"),
        (2203, 2437, 16, 3, 1508647, "saturated"),
    )
    for p, q, m, k, beta, expected in cases:
        n = p*q
        r = p % m
        g = gcd(r-1, m)
        lam = hh_parameters(n, m, k)
        giants = list(hh_giants(n, m, r, k, lam))
        stride = m*m
        assert all(exp*g == a*n+b-j and (a*q-b) % g == 0 for a,b,j,exp in giants)
        assert any((a*q+b*p-j) % stride == 0 and 0 <= (a*q+b*p-j)//stride <= 2*lam
                   for a,b,j,exp in giants)
        stats = {}
        value = hh_collision(n, m, beta, g, lam, giants, stats=stats)
        assert stats["status"] == expected
        assert value in ((None,) if expected == "short_order" else (p,q))
        if expected == "two_part":
            assert n_order(beta, n) > m*(2*lam+1)
            assert n_order(pow(beta, stride//g, n), n) <= 2*lam+1


def check_medium():
    for p, q, m in ((850000013, 1345000021, 128), (17179869209, 25769803799, 128),
                    (17179869209, 25769803799, 2)):
        n = p*q
        stats = {}
        assert hh_low_factor(n, m, p % m, stats=stats) in (p,q)
        assert stats["status"] == "collision"
        print(f"{n.bit_length()} bits, modulus {m}, {stats['baby_count']} babies, "
              f"{stats['giant_count']} giants, recovered a divisor by product evaluation")


if __name__ == "__main__":
    check_orders()
    check_small()
    check_collisions()
    check_medium()
