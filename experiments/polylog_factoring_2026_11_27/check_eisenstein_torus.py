"""Compare native cubic symbols to known-prime Euler values and enumerate torus fibers."""

from collections import Counter
from itertools import combinations, product
from math import gcd
from random import Random

from sympy import factorint, nextprime, primerange, sqrt_mod

if __package__:
    from .eisenstein_torus import (OMEGA, ONE, cubic_jacobi, mixed_probe, norm, rotation_labels,
                                   torus_mark, torus_power, torus_product, torus_sample)
else:
    from eisenstein_torus import (OMEGA, ONE, cubic_jacobi, mixed_probe, norm, rotation_labels,
                                  torus_mark, torus_power, torus_product, torus_sample)


def euler_symbol(r, z):
    """Known-prime reference from the two split residue fields or the inert residue field."""
    if r % 3 == 2:
        value = torus_power(r, z, (r*r-1)//3)
        return (ONE, OMEGA, (r-1, r-1)).index(value)
    root = int((-1+sqrt_mod(-3, r))*pow(2, -1, r) % r)
    labels = []
    for w in (root, root*root % r):
        value = pow((z[0]+z[1]*w) % r, (r-1)//3, r)
        labels.append((1, w, w*w % r).index(value))
    return sum(labels) % 3


def check_raw():
    tested = 0
    for n in range(2, 96):
        if n % 3 == 0:
            continue
        factors = factorint(n)
        for z in product(range(-4, 5), repeat=2):
            expected = None if gcd(norm(z), n) > 1 else sum(
                e*euler_symbol(p, z) for p, e in factors.items()) % 3
            assert cubic_jacobi(n, z) == expected, (n, z, expected)
            tested += 1
    return tested


def crt(x, y, p, q):
    return tuple(a+p*((b-a)*pow(p, -1, q) % q) for a, b in zip(x, y))


def local_torus(r):
    return [z for z in product(range(r), repeat=2) if norm(z) % r == 1]


def check_tori():
    tested = shifts = 0
    for p, q in combinations(primerange(5, 30), 2):
        n = p*q
        ep = 1 if p % 3 == 1 else -1
        eq = 1 if q % 3 == 1 else -1
        histogram = Counter()
        for x, y in product(local_torus(p), local_torus(q)):
            z = crt(x, y, p, q)
            lp, lq = torus_mark(p, x), torus_mark(q, y)
            label = cubic_jacobi(n, z)
            assert label == (-ep*lp-eq*lq) % 3
            assert (label == 0) == ((lp == lq) if n % 3 == 2 else ((lp+lq) % 3 == 0))
            histogram[label, lp == lq] += 1
            rotations = rotation_labels(n, z)
            for j, rotated_label in enumerate(rotations):
                rotated = torus_product(n, z, torus_power(n, OMEGA, j))
                assert cubic_jacobi(n, rotated) == rotated_label
                shifted = (1+rotated[0], rotated[1])
                if gcd(norm(shifted), n) == 1:
                    assert cubic_jacobi(n, shifted) == (-rotated_label) % 3
                    shifts += 1
            tested += 1
        for label in range(3):
            diagonal = histogram[label, True]
            total = diagonal+histogram[label, False]
            assert total == (p-ep)*(q-eq)//3
            if n % 3 == 1:
                assert 3*diagonal == total
            else:
                assert diagonal == (total if label == 0 else 0)
    return tested, shifts


def check_probes():
    rng = Random(89331)
    output = []
    # The large controls use primes below 2^64, where SymPy's primality test is deterministic.
    p = int(nextprime(1 << 59))
    while p % 9 != 2:
        p = int(nextprime(p))
    for q_class in (1, 7, 4):
        q = int(nextprime(3 << 58))
        while q % 9 != q_class:
            q = int(nextprime(q))
        n = p*q
        z, immediate = mixed_probe(n, rng)
        assert immediate is None and norm(z) % n == 1
        lp, lq = torus_mark(p, z), torus_mark(q, z)
        assert lp != lq and cubic_jacobi(n, z) in (1, 2)
        output.append((n, n % 9, z, cubic_jacobi(n, z), (lp, lq), torus_power(n, z, 3) == ONE))
    z, factor = torus_sample(35, Random(773))
    assert (factor in (5, 7)) if factor is not None else (norm(z) % 35 == 1)
    return output


if __name__ == "__main__":
    print("Native symbols against factored Euler symbols:", check_raw())
    print("Complete torus-element comparisons and unit-shift checks:", check_tori())
    for row in check_probes():
        print("N, N mod 9, public witness, cubic label, prime-check endpoints, witness has cube one:", row)
