"""Compare the full ideal-norm product with a character-filtered trial product."""

from functools import cache
from itertools import combinations
from math import gcd, isqrt, prod

from sympy import divisors, factorint, kronecker_symbol, primerange


def fundamental_discriminant(a):
    """Discriminant of the field Q(sqrt(-a)), including nonsquarefree a."""
    squarefree = -prod(int(p) for p, exponent in factorint(a).items() if exponent % 2)
    return squarefree if squarefree % 4 == 1 else 4*squarefree


@cache
def character(delta, n):
    return int(kronecker_symbol(delta, n))


@cache
def ideal_count(delta, n):
    """Count index-n ideals independently from their rank-two lattice bases.

    Write omega^2 = trace * omega + constant. Each ideal has a unique
    basis (c*h, c*(z+omega)), c^2*h = n, 0 <= z < h. Closure under omega
    requires z*(z+trace) == constant mod h.
    """
    trace = delta % 4
    constant = (delta - trace*trace) // 4
    count = 0
    for c in range(1, isqrt(n)+1):
        if n % (c*c) == 0:
            h = n // (c*c)
            count += sum((z*(z+trace) - constant) % h == 0 for z in range(h))
    return count


def audit_coefficients(delta, limit):
    for n in range(1, limit+1):
        expected = sum(character(delta, int(d)) for d in divisors(n))
        assert ideal_count(delta, n) == expected
        if delta in (-4, -8):
            # Both are PIDs. Count all generators, then divide by the unit group.
            a = -delta//4
            representations = 0
            for j in range(-isqrt(n//a), isqrt(n//a)+1):
                rem = n-a*j*j
                k = isqrt(rem)
                if k*k == rem:
                    representations += 1 if k == 0 else 2
            units = 4 if delta == -4 else 2
            assert representations == units*expected


def evaluate_products(n, delta):
    """Input-only residues. The prime factors are used only by the fixture caller."""
    ideal_product = character_product = 1
    for i in range(1, isqrt(2*n)+1):
        ideal_product = (ideal_product * pow(i, ideal_count(delta, i), n)) % n
        if character(delta, i) == 1:
            character_product = (character_product * i) % n
    return gcd(ideal_product, n), gcd(character_product, n)


def check():
    params = (0, 1, 2, 3, 4, 7, 8, 12)
    primes = list(primerange(3, 100))
    deltas = {x: fundamental_discriminant(x*x+1) for x in params}
    assert deltas[2] == -20 and deltas[4] == -68 and deltas[7] == -8
    for delta in set(deltas.values()):
        audit_coefficients(delta, isqrt(2*primes[-1]**2))

    # Nonprincipal example. There are two ideals of norm 23, while 23
    # is absent from the principal form i^2+17*j^2 in Q(sqrt(-17)).
    p, q, delta = 19, 23, -68
    assert ideal_count(delta, q) == 2 and ideal_count(delta, p) == 0
    assert all(isqrt(q-17*j*j)**2 != q-17*j*j
               for j in range(isqrt(q//17)+1))
    assert evaluate_products(p*q, delta) == (q, q)

    cases = 0
    for p, q in combinations(primes, 2):
        if q > 2*p:
            continue
        n = p*q
        assert q <= isqrt(2*n) < 2*p
        for x, delta in deltas.items():
            if gcd(n, x*x+1) != 1 or character(delta, n) != -1:
                continue
            split = p if character(delta, p) == 1 else q
            assert evaluate_products(n, delta) == (split, split)
            cases += 1
    print(f"{cases} mixed balanced fixtures, {len(set(deltas.values()))} field discriminants")
    print("N=437: two nonprincipal ideals of norm 23; both full-menu products give gcd 23")


if __name__ == "__main__":
    check()
