"""Check the batch, including all parameter choices across ordinary small prime cohorts."""

from itertools import combinations
from math import gcd, isqrt

from sympy import isprime, jacobi_symbol, primerange

if __package__:
    from .large_conic_batch import NormHit, batch_conic_factor
else:
    from large_conic_batch import NormHit, batch_conic_factor


def direct(n, x):
    norm, limit = x*x + 1, isqrt(8*n)
    for b in range(1, isqrt(limit // norm) + 1):
        for a in range(1 + isqrt(limit - norm*b*b)):
            divisor = gcd(a*a + norm*b*b, n)
            if 1 < divisor < n:
                return NormHit(divisor, a, b)
    return None


def first_multipliers(prime, parameters):
    """Smallest multiplier e<=3 represented at each x, grading with the known prime."""
    admissions = {}
    for x in parameters:
        norm = x*x + 1
        for e in (1, 2, 3):
            target = e*prime
            if any(isqrt(target - norm*b*b)**2 == target - norm*b*b
                   for b in range(1, 1 + isqrt(target // norm))):
                admissions[x] = e
                break
    return admissions


def cohort(bound, *, run_native=False):
    primes = list(primerange(bound, 2*bound))
    start = isqrt(bound - 1) + 1
    parameters = range(start, 2*start)
    tables = {p: first_multipliers(p, parameters) for p in primes}
    represented = sum(any(e <= 2 for e in row.values()) for row in tables.values())
    prime_x = sum(e <= 2 for row in tables.values() for e in row.values())
    pairs = actual_pairs = genuine_pairs = actual_settings = native_checks = 0
    for p, q in combinations(primes, 2):
        n, limit = p*q, isqrt(8*p*q)
        genuine = set()
        actual = set()
        for prime in (p, q):
            for x, e in tables[prime].items():
                if e <= 2:
                    genuine.add(x)
                if e*prime <= limit:
                    actual.add(x)
        assert genuine <= actual
        if run_native:
            for x in parameters:
                result = batch_conic_factor(n, x)
                assert result == direct(n, x) and bool(result) == (x in actual)
                native_checks += 1
        pairs += 1
        actual_pairs += bool(actual)
        genuine_pairs += bool(genuine)
        actual_settings += len(actual)
    return (bound, len(primes), represented, prime_x, pairs, genuine_pairs,
            actual_pairs, actual_settings, native_checks)


def main():
    n, x = 10412720415486989, 8192
    work = {}
    hit = batch_conic_factor(n, x, stats=work)
    assert hit == NormHit(83317541, 4026, 1)
    assert n // hit.factor == 124976329 and all(map(isprime, (hit.factor, n // hit.factor)))
    assert jacobi_symbol(-x*x-1, n) == -1
    assert hit.a**2 + (x*x+1)*hit.b**2 == hit.factor
    print(f"planted {n=} {x=} {hit=} {work=} t_scan={2*n // (x*x+1)}")
    assert batch_conic_factor(91, 4, stats=work) == direct(91, 4) == NormHit(7, 2, 1)
    assert work["saturated_blocks"] == 1
    print("saturated block at n=91, x=4 gives gcd(21,91)=7")
    print("R primes represented_primes represented_prime_x pairs genuine_pairs "
          "actual_pairs actual_pair_x checked_natively")
    for bound in (64, 256, 1024, 4096):
        print(*cohort(bound, run_native=bound == 64))


if __name__ == "__main__":
    main()
