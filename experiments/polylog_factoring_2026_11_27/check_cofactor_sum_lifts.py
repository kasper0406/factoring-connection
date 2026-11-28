"""Compare dyadic sum congruences, signed-search ranks and actual factors."""

from itertools import combinations
from math import gcd, isqrt
from random import Random

from sympy import isprime, nextprime, primerange

if __package__:
    from .cofactor_sum_lifts import signed_fermat_factor, state_sum_class
else:
    from cofactor_sum_lifts import signed_fermat_factor, state_sum_class


def check_precision():
    rng = Random(591)
    checked = 0
    for bits in range(1, 8):
        m = 1 << bits
        high = 16*m*m
        for _ in range(40):
            n = rng.randrange(1, high, 2)
            a = rng.randrange(1, m, 2)
            base = a+n*pow(a, -1, high)
            stride = high
            sums = set()
            for x in range(a, high, m):
                s = (x+n*pow(x, -1, high)) % high
                stride = gcd(stride, s-base)
                sums.add(s)
            assert state_sum_class(n, a, bits) == (base % stride, stride)
            if bits >= 2:
                period = 8*m*m
                assert any((s+period) % high in sums for s in sums)
            checked += 1
    return checked


def ranks(p, q):
    a, c = (p+q)//2, (q-p)//2
    return (1+(a-isqrt(p*q-1)-1)//(c & -c)**2,
            1+c//(a & -a)**2)


def check_semiprimes():
    checked = 0
    primes = list(primerange(3, 230))
    for p, q in combinations(primes, 2):
        best = min(ranks(p, q))
        for budget in {1, best, max(1, best-1)}:
            found = signed_fermat_factor(p*q, budget)
            assert (found in (p, q)) == (best <= budget)
            checked += 1
    return checked


def cohorts():
    rng = Random(74175)
    rows = []
    for prime_bits in (48, 96, 160):
        e = prime_bits//3
        stride = 1 << e
        for direction in (1, -1):
            p = int(nextprime((1 << (prime_bits-1))+rng.getrandbits(prime_bits-3)))
            # The gap is about 2^(3k/5). Choose q in the desired signed state,
            # with the exact v2 of q-p or q+p equal to e.
            gap = 1 << (3*prime_bits//5)
            start = p+gap
            q = start+(direction*p+stride-start) % (2*stride)
            while not isprime(q):
                q += 2*stride
            n = p*q
            assert p < q < 2*p and stride**4 < n
            plus, minus = ranks(p, q)
            assert (plus if direction == 1 else minus) == 1
            assert signed_fermat_factor(n, 1) == p
            ordinary = 1+(p+q)//2-(isqrt(n-1)+1)
            rows.append((prime_bits, direction, e, ordinary, plus, minus))
    background = []
    for bits in (16, 32, 64):
        hits = 0
        for _ in range(100):
            p, q = sorted(int(nextprime((1 << (bits-1))+rng.getrandbits(bits-2)))
                          for _ in range(2))
            while p == q:
                q = int(nextprime(q))
            found = signed_fermat_factor(p*q, 1)
            assert (found is not None) == (min(ranks(p, q)) == 1)
            hits += found is not None
        background.append((bits, hits))
    return rows, background


if __name__ == "__main__":
    print("Maximal state congruences verified:", check_precision())
    print("Exact semiprime budget predictions:", check_semiprimes())
    planted, background = cohorts()
    print("Planted controls (prime bits, sign, shared bits, ordinary rank, plus rank, minus rank):")
    for row in planted:
        print(row)
    print("Unplanted one-point hits among 100 per prime size:", background)
