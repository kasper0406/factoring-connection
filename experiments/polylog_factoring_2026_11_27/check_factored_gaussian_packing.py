"""Independent checks for the Gaussian prime pair construction."""

from math import gcd, prod
from random import Random

from sympy import factorint, jacobi_symbol, nextprime

if __package__:
    from .factored_gaussian_packing import certified_prime, factored_packing, line_octant
    from .check_fixed_field_packing import check_recovery, sphere_progression
else:
    from factored_gaussian_packing import certified_prime, factored_packing, line_octant
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


def check_challenge(n, x, seed):
    result = factored_packing(n, x, Random(seed))
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
    assert c * 8 * (x+1) < params.B and result.t >= 2
    assert params.M == (params.P*x + params.Q*n)**2 + params.P**2 == a * (p*q)**2
    # The only unrestricted factorization here concerns the short base norm A.
    factored = dict(factorint(a))
    factored[p] = factored[q] = 2
    assert result.norm_factorization == factored
    assert prod(prime**e for prime, e in factored.items()) == params.M
    print(f"N bits={n.bit_length()} x={x} prime bits={p.bit_length()},{q.bit_length()} "
          f"M bits={params.M.bit_length()} lines={result.lines} t={result.t}")
    return result


if __name__ == "__main__":
    print(f"Direct annulus comparisons: {check_lines()}")
    for n, x in ((15, 0), (35, 0), (65, 4), (77, 1), (91, 0), (183, 0)):
        result = check_challenge(n, x, 41027)
        if n == 15:
            decisions = check_recovery(n, x, 3, (result.parameters.B, result.t), primitive=True)
            print(f"  N=15 input-only sphere decisions: {dict(decisions)}")
        if n == 77:
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
