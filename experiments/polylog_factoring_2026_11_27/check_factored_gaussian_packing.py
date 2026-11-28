"""Independent checks for the Gaussian prime pair construction."""

from math import gcd, isqrt, prod
from random import Random

from sympy import factorint, jacobi_symbol, nextprime

if __package__:
    from .factored_gaussian_packing import _enough_prime_vectors, _search_lines, certified_prime, factored_packing, line_octant
    from .check_fixed_field_packing import check_recovery, sphere_progression
else:
    from factored_gaussian_packing import _enough_prime_vectors, _search_lines, certified_prime, factored_packing, line_octant
    from check_fixed_field_packing import check_recovery, sphere_progression


def check_lines():
    cases = 0
    for n in (9, 15, 25, 35):
        for a in range(n):
            for b in range(n):
                if gcd(a, b, n) != 1:
                    continue
                height = n + 5
                expected = {
                    (x, y) for x in range(1, height + 1) for y in range(1, x)
                    if height**2 < 2 * (x*x + y*y) <= 2 * height**2
                    and gcd(x, y) == 1 and (a*y - b*x) % n == 0
                }
                actual = list(line_octant(n, a, b, height))
                assert len(actual) == len(expected) and set(actual) == expected
                cases += 1
    return cases


def check_challenge(n, x, seed, *, positive_scale=False):
    result = factored_packing(n, x, Random(seed), positive_scale=positive_scale)
    a = x*x + 1
    params = result.parameters
    (i, j), (k, l) = result.gaussian_primes
    c = n * result.t
    p, q = result.primes
    assert sorted((p, q)) == sorted((i*i + j*j, k*k + l*l))
    assert p != q and certified_prime(p) and certified_prime(q)
    assert gcd(p*q, 2*a*n*params.Q) == 1
    assert 0 < j < i and 0 < l < k and i % 2 == k % 2
    assert (i*k+j*l, abs(i*l-j*k)) == (params.B, c)
    assert c < params.B and result.t >= 2 and result.t % 2 == 0
    if positive_scale:
        assert c * 8 * (x+1) < params.B
    else:
        assert result.height == 1 + isqrt(64*n*n.bit_length()-1)
        assert params.T <= result.height**4
    assert params.M == (params.P*x + params.Q*n)**2 + params.P**2 == a * (p*q)**2
    # The only unrestricted factorization here concerns the short base norm A.
    factored = dict(factorint(a))
    factored[p] = factored[q] = 2
    assert result.norm_factorization == factored
    assert prod(prime**e for prime, e in factored.items()) == params.M
    print(f"N bits={n.bit_length()} x={x} prime bits={p.bit_length()},{q.bit_length()} "
          f"M bits={params.M.bit_length()} lines={result.lines} t={result.t}")
    return result


def check_retries():
    n, x, height = 15, 14, 62
    assert not _enough_prime_vectors(n, 2*(x*x+1)*n, 1, certified_prime)
    assert _enough_prime_vectors(n, 2*(x*x+1)*n, height, certified_prime)
    widened = _search_lines(n, x, Random(41027), certified_prime, 1, 1)
    assert widened.height > 1 and widened.height.bit_count() == 1

    class FirstBatchMisses(Random):
        def __init__(self):
            super().__init__(41027)
            self.misses = 2 * 16 * height.bit_length()**3

        def randrange(self, n):
            if self.misses:
                self.misses -= 1
                return 0
            return super().randrange(n)

    held = factored_packing(n, x, FirstBatchMisses())
    assert held.height == height and held.parameters.P < 0
    print(f"Height doubling from 1 reached {widened.height}; sufficient height {height} retained")


if __name__ == "__main__":
    print(f"Direct annulus comparisons: {check_lines()}")
    check_retries()
    negative = check_challenge(15, 14, 41027)
    assert negative.parameters.P < 0
    decisions = check_recovery(15, 14, 5, (negative.parameters.B, negative.t), primitive=True)
    print(f"  negative-P native recovery: {dict(decisions)}")
    for n, x in ((15, 0), (35, 0), (65, 4), (77, 1), (91, 0), (183, 0)):
        result = check_challenge(n, x, 41027)
        if n == 15:
            decisions = check_recovery(n, x, 3, (result.parameters.B, result.t), primitive=True)
            print(f"  N=15 input-only sphere decisions: {dict(decisions)}")
        if n == 77:
            result = check_challenge(n, x, 41027, positive_scale=True)
            params = result.parameters
            assert result.primes == (56269, 42089)
            modulus = params.Q * 17
            yes, no = (2555529, 1776889, 3503880), (4112809, 219609, 5061160)
            assert (sum(i*i for i in yes) - params.M) % modulus == 0
            assert (sum(i*i for i in no) - params.M) % modulus == 0
            w, = sphere_progression(params.M, yes, modulus)
            assert not list(sphere_progression(params.M, no, modulus))
            v = tuple((params.P*u-i) // params.Q for u, i in zip((1, 1, 0), w))
            assert v == (-14, -63, 42) and gcd(*v) == 7
            print(f"  N=77 exact sphere counts 1 and 0; factored radius data {result.primes}")
    for bits in (12, 18, 24, 32, 40):
        n = int(nextprime(2**bits + 412)) * int(nextprime(2**bits + 908))
        x = next(k for k in range(101) if jacobi_symbol(-(k*k+1), n) == -1)
        check_challenge(n, x, 41027)
    n = 32416187567 * 32416190071
    check_challenge(n, 2, 41027)
