"""Check public diagonal frequencies and factoring frequencies for CRT quarter parity."""

from collections import defaultdict
from itertools import product
from math import cos, gcd, pi, prod, sin
from random import Random


def sign_vectors(count):
    return [(1, *rest) for rest in product((1, -1), repeat=count-1) if -1 in rest]


def smooth_coefficients(n, limit):
    size = 2*(n//4)+1
    return {
        h: ((2*size-n)/n if h == 0 else
            (1-abs(h)/(limit+1))*2*sin(pi*size*h/n)/(n*sin(pi*h/n)))
        for h in range(-limit, limit+1)
    }


def vector_spectrum(n, count, limit):
    spectrum = {(0,)*count: 1.0}
    coefficients = smooth_coefficients(n, limit)
    for signs in sign_vectors(count):
        new = defaultdict(float)
        for vector, coefficient in spectrum.items():
            for h, weight in coefficients.items():
                target = tuple(value+h*sign for value, sign in zip(vector, signs))
                new[target] += coefficient*weight
        spectrum = new
    return spectrum


def check(factors, limit, rng):
    n = prod(factors)
    idempotents = [n//r*pow(n//r, -1, r) for r in factors]
    vectors = vector_spectrum(n, len(factors), limit)
    public, factorable = {}, {}
    for vector, weight in vectors.items():
        frequency = sum(t*x for t, x in zip(vector, idempotents)) % n
        if len(set(vector)) == 1:
            assert frequency == vector[0] % n
            public[frequency] = weight
        else:
            assert any(1 < gcd(frequency-t, n) < n for t in vector)
            factorable[frequency] = weight
    assert not (public.keys() & factorable.keys())
    roots = [sum(t*x for t, x in zip(vector, idempotents)) % n
             for vector in sign_vectors(len(factors))]
    coeff = smooth_coefficients(n, limit)
    for _ in range(20):
        c = rng.randrange(n)
        product_value = prod(sum(a*cos(2*pi*(c*h*u % n)/n) for h, a in coeff.items())
                             for u in roots)
        sum_value = sum(a*cos(2*pi*(c*f % n)/n) for part in (public, factorable)
                        for f, a in part.items())
        assert abs(product_value-sum_value) < 1e-10
    print(n, len(public), len(factorable))


if __name__ == "__main__":
    rng = Random(261129)
    check((101, 109, 127), 7, rng)
    check((103, 113, 131, 137), 2, rng)
    check((101**2, 109, 127), 7, rng)
