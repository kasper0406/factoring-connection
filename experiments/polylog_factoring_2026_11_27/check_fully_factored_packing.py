"""Check the prime-conductor construction, complete factorizations, and scalar decisions."""

from collections import Counter
from math import gcd, isqrt, prod
from random import Random

from sympy import jacobi_symbol, nextprime, primerange

if __package__:
    from .factored_gaussian_packing import certified_prime
    from .fixed_field_packing import recover_factor
    from .fully_factored_packing import fully_factored_packing
    from .packed_ternary_decision import brute_sphere_count
else:
    from factored_gaussian_packing import certified_prime
    from fixed_field_packing import recover_factor
    from fully_factored_packing import fully_factored_packing
    from packed_ternary_decision import brute_sphere_count


def check_transcript(n, x, result):
    p = result.parameters
    b, q, a = p.B, result.auxiliary_prime, x*x+1
    for value, factors in ((p.Q, result.scale_factorization), (p.M, result.norm_factorization)):
        assert all(certified_prime(s) and e > 0 for s, e in factors.items())
        assert prod(s**e for s, e in factors.items()) == value
    assert b > 2*n and p.T == b*b+4*n*n*q*q and p.Q == 4*a*b*q
    assert (x*p.P+n*p.Q)**2+p.P**2 == p.M == a*p.T**2
    assert gcd(p.P, p.Q) == gcd(p.T, 2*a*n) == 1
    if p.P < 0:
        assert n*n*(2*q) + 2*x*n*b < a*b*b*(2*q)


def cohort():
    for bits in (16, 24, 32, 40):
        boxes, trials = Counter(), []
        for seed in range(8):
            rng = Random(813731*bits + 2441*seed)
            lo, hi = isqrt(2**(bits-1))+1, 2**(bits//2)
            primes = set()
            while len(primes) < 2:
                p = int(nextprime(rng.randrange(lo, hi)))
                if p < hi:
                    primes.add(p)
            n = prod(primes)
            assert n.bit_length() == bits
            x = next(k for k in range(32) if jacobi_symbol(-(k*k+1), n) == -1)
            result = fully_factored_packing(n, x, rng)
            check_transcript(n, x, result)
            boxes[result.box] += 1
            trials.append(result.prime_pairs)
        print(f"{bits}-bit inputs, eight seeds, boxes {dict(boxes)}, prime-pair trials {trials}")


def check_recovery():
    for n, x, factor in ((15, 0, 3), (65, 4, 5), (77, 1, 7)):
        result = fully_factored_packing(n, x, Random(90017+n))
        check_transcript(n, x, result)
        hits = Counter()

        def decision(norm, modulus, residue):
            p = result.parameters
            ell, remainder = divmod(modulus, p.Q)
            assert norm == p.M and remainder == 0 and certified_prime(ell)
            assert ell not in result.scale_factorization
            count = brute_sphere_count(norm, modulus, residue)
            hits[count] += 1
            return count > 0

        got = recover_factor(n, x, decision, gaussian=(result.parameters.B, 2*result.auxiliary_prime))
        assert got == factor and 0 in hits and 1 in hits
        print(f"N={n}, x={x}, T={result.parameters.T}, Q={result.parameters.Q}, "
              f"factor={got}, scalar queries by multiplicity {dict(sorted(hits.items()))}")


def wide_transcripts():
    for bits in (64, 96):
        bound = 2**(bits//2)
        n = int(nextprime(bound-20000))*int(nextprime(bound-10000))
        x = next(j for j in range(32) if jacobi_symbol(-(j*j+1), n) == -1)
        record = fully_factored_packing(n, x, Random(431747+bits))
        check_transcript(n, x, record)
        assert n.bit_length() == bits
        print(f"{bits}-bit input, {record.parameters.T.bit_length()}-bit certified conductor, "
              f"{record.parameters.M.bit_length()}-bit squared norm, {record.prime_pairs} prime pairs")


def exact_box():
    n, x = 251*257, 0
    pairs = hits = 0
    assert jacobi_symbol(-(x*x+1), n) == -1
    for b in primerange(2*n, 4*n+1):
        for q in (17, 19, 23, 29, 31):
            pairs += 1
            hits += certified_prime(b*b + 4*n*n*q*q)
    assert (hits, pairs) == (4157, 52995)
    print(f"All first-box prime pairs at N={n} had {hits} prime conductors in {pairs} trials")


if __name__ == "__main__":
    cohort()
    exact_box()
    wide_transcripts()
    check_recovery()
