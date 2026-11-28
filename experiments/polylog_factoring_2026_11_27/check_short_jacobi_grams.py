"""Exact squareclass controls and a held-out balanced-semiprime selector cohort."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import gcd, isqrt, prod
from random import Random

from sympy import factorint, isprime, jacobi_symbol, primerange

if __package__:
    from .short_jacobi_grams import first_vector_orbits, select_x
else:
    from short_jacobi_grams import first_vector_orbits, select_x


def mask(value, indices):
    row = int(value < 0)
    for prime, exponent in factorint(abs(value)).items():
        if exponent % 2:
            row ^= 1 << indices.setdefault(prime, len(indices)+1)
    return row


def reduce_row(row, basis):
    while row and row.bit_length() in basis:
        row ^= basis[row.bit_length()]
    return row


def check_triangle():
    indices, basis = {}, {}
    assert (-2)*(-5)*(-10) == -10**2
    triangular = 0
    for y in range(1, 129):
        value = -(y*y+1)
        row = reduce_row(mask(value, indices), basis)
        if row:
            basis[row.bit_length()] = row
        if y < 3:
            continue
        assert reduce_row(1, basis) == 0
        assert reduce_row(mask(2, indices), basis) == 0
        for ell in primerange(5, 2*y+2):
            if ell % 4 == 1:
                r = next(x for x in range(1, (ell+1)//2) if (x*x+1) % ell == 0)
                assert r <= y and 1 <= (r*r+1)//ell < ell
                assert reduce_row(mask(ell, indices), basis) == 0
                triangular += 1
    return triangular


def check_fixed_density():
    # For fixed y, PNT in arithmetic progressions gives the limiting prime
    # sign law. Check the exact finite law on units modulo a common conductor.
    values = [-1, -2, -5, -10, -17]
    modulus = 8*5*17
    indices, basis, units = {}, {}, []
    for v in values:
        row = reduce_row(mask(v, indices), basis)
        if row:
            basis[row.bit_length()] = row
    for a in range(1, modulus, 2):
        if gcd(a, modulus) == 1:
            units.append(tuple(int(jacobi_symbol(v, a)) for v in values))
    counts = Counter(units)
    patterns = 1 << len(basis)
    assert len(counts) == patterns and len(set(counts.values())) == 1
    assert patterns*sum(count*count for count in counts.values()) == len(units)**2
    return len(units), patterns


def mean_partial_sum():
    indices, basis = {}, {}
    result = Fraction(1)
    for x in range(121):
        row = reduce_row(mask(-(x*x+1), indices), basis)
        if row:
            basis[row.bit_length()] = row
        result += Fraction(1, 2**len(basis))
    assert len(basis) == 86
    return result


def local_sign(value, p):
    k = pow(value, (p-1)//2, p)
    return -1 if k == p-1 else k


def reference(p, q, limit):
    for x in range(limit+1):
        value = -(x*x+1)
        sign = local_sign(value, p)*local_sign(value, q)
        if sign == -1:
            return x, None
        if sign == 0 and gcd(value, p*q) < p*q:
            return None, gcd(value, p*q)
    return None, None


def check_selector():
    primes = list(primerange(3, 120))
    cases = 0
    for p, q in combinations(primes, 2):
        for limit in (0, 2, 12):
            assert select_x(p*q, limit) == reference(p, q, limit)
            cases += 1
        assert select_x(p*q, (p+q-2)//2) != (None, None)
    assert select_x(65, 4) == (None, 5)
    # Selected by finding a repeated 26-coordinate sign pattern among
    # consecutive primes, independently of all cohort seeds.
    p, q = 73477, 80557
    assert isprime(p) and isprime(q)
    assert select_x(p*q, 26) == (None, None)
    assert select_x(p*q, 27) == reference(p, q, 27) == (27, None)
    return cases


def check_orbits():
    for x in (*range(12), 27):
        norm = x*x+1
        b = isqrt(norm)
        expected = {v for v in product(range(-b, b+1), repeat=3)
                    if sum(a*a for a in v) == norm}
        shapes = list(first_vector_orbits(x))
        assert shapes == sorted(set(shapes))
        assert all(0 <= a <= b <= c for a, b, c in shapes)
        rebuilt = {tuple(s*a for s, a in zip(signs, v))
                   for shape in shapes for v in permutations(shape)
                   for signs in product((-1, 1), repeat=3)}
        assert rebuilt == expected


def check_conic():
    n, x = 35, 4
    u, v = (2, -3, -2), (-10, 30, 15)
    def dot(s, t):
        return sum(a*b for a, b in zip(s, t))
    a = x*x+1
    assert (dot(u, u), dot(u, v), dot(v, v)) == (a, -n*x, n*n)
    assert jacobi_symbol(-a, n) == -1
    target = (x, 1, 0)
    delta = tuple(i-j for i, j in zip(u, target))
    d = a-dot(u, target)
    w = tuple(d*j-k*dot(delta, v) for j, k in zip(v, delta))
    t, z = n*d+w[0], w[2]
    assert gcd(d, n) == 1 and w[1] == -x*t and z*z == t*(2*d*n-a*t)
    common = gcd(t, z)
    s, h = t//common, z//common
    norm = a*s*s+h*h
    assert (d, s, h, norm) == (12, 1, -2, 21)
    assert 2*d*n % norm == 0 and gcd(norm, n)*gcd(*v) == n


def draw_prime(rng, lower, upper):
    while True:
        p = rng.randrange(lower | 1, upper+1, 2)
        # isprime performs deterministic testing throughout this range.
        if isprime(p):
            return p


def cohort():
    seeds = (94183, 198017, 536339, 873463)
    for bits in (48, 64, 96, 128):
        lower = isqrt(2**(bits-1)-1)+1
        upper = 2**(bits//2)-1
        assert upper < 2**64 and upper < 2*lower
        visited = set()
        histogram = Counter()
        for seed in seeds:
            rng = Random(seed*256+bits)
            for _ in range(64):
                pair = []
                while len(pair) < 2:
                    p = draw_prime(rng, lower, upper)
                    if p not in visited:
                        visited.add(p)
                        pair.append(p)
                n = prod(pair)
                assert n.bit_length() == bits
                result = select_x(n, bits)
                assert result == reference(*pair, bits)
                x, factor = result
                assert factor is None
                histogram[x if x is not None else 'exhausted'] += 1
        assert len(visited) == 512 and sum(histogram.values()) == 256
        print(f"{bits}-bit, first x and frequency {dict(histogram.most_common())}")


if __name__ == '__main__':
    triangles = check_triangle()
    units, patterns = check_fixed_density()
    selectors = check_selector()
    check_orbits()
    check_conic()
    print(f"{triangles} triangular checks, {units} units in {patterns} equally sized sign classes, "
          f"{selectors} selector controls and 13 exact signed-orbit checks")
    print("Adversarial N=5919086689, first accepted x=27")
    mean = mean_partial_sum()
    print(f"Mean series through x=120: {mean} = {float(mean):.14f}")
    cohort()
