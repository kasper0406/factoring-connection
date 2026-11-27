"""Check subgroup grouping and direct gcds on planted sparse quarter rules."""

from itertools import combinations
from math import gcd, lcm

from check_quarter_fourier import quarter


def subgroup_order(prime, exponents):
    order = prime-1
    degrees = list(dict.fromkeys([0, 1, *exponents]))
    denominator = 4*(len(exponents)+1)
    while True:
        for left, right in combinations(degrees, 2):
            common = gcd(left-right, order)
            if common < order and common**denominator*prime > order**denominator:
                order = common
                break
        else:
            return order


def direct_gcd(n, c, terms, h, k, offset=0):
    values = [coefficient*pow(c, exponent-1, n) % n for coefficient, exponent in terms]
    for mask in range(1 << len(values)):
        target = h*sum(value for i, value in enumerate(values) if mask & (1 << i))
        divisor = gcd(target-k+offset, n)
        if 1 < divisor < n:
            return divisor
    return None


def check(p, q):
    n = p*q
    hidden = 1+p*(-2*pow(p, -1, q) % q)
    e = 1+lcm((p-1)//2, (q-1)//2)
    terms = [(1, e), (2, 0)]
    correlation = roots = units = 0
    for c in range(1, n):
        if gcd(c, n) != 1:
            continue
        units += 1
        polynomial = sum(a*pow(c, k, n) for a, k in terms) % n
        correlation += quarter(n, polynomial)*quarter(n, c*hidden)
        if pow(c, e-1, n) in (hidden, n-hidden):
            assert direct_gcd(n, c, terms, 1, 1) in (p, q)
            roots += 1
    assert roots == units//2
    print(n, correlation, units, roots)


def check_intersection(offset=0):
    p, q = 101, 109
    n = p*q
    hidden = 1+p*(-2*pow(p, -1, q) % q)
    first = 2*q*pow(q, -1, p)
    second = -2*p*pow(p, -1, q)
    terms = [(-hidden-offset, 1), (first, p), (second, q)]
    exponents = [exponent for _, exponent in terms]
    local = [
        {i for i, e in enumerate(exponents) if (e-1) % subgroup_order(r, exponents) == 0}
        for r in (p, q)
    ]
    assert local == [{0, 1}, {0, 2}]
    for c in (1, 2, 5, 51, 1099):
        values = [a*pow(c, e-1, n) % n for a, e in terms]
        assert (sum(values[i] for i in local[0])+offset) % n == 1
        assert (sum(values[i] for i in local[1])+offset) % n == n-1
        assert gcd(sum(values[i] for i in local[0] & local[1])-1+offset, n) == q
        assert direct_gcd(n, c, terms, 1, 1, offset) in (p, q)


if __name__ == "__main__":
    for p, q in [(31, 43), (103, 107), (211, 239)]:
        check(p, q)
    check_intersection(2*pow(3,-1,101*109) % (101*109))
    check_intersection()
    print("The intersection test handled distinct local exponent sets.")
