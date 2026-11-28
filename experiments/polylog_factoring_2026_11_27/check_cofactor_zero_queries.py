"""Check whole-modulus Boolean decoders with direct products and large valuation controls."""

from functools import cache
from math import prod
from random import Random
from sympy import factorint, nextprime

if __package__:
    from .check_cofactor_bit_binomials import direct_ratio
    from .cofactor_zero_queries import factorization, largest_prime, semiprime_from_pairs, zero_from_complement
else:
    from check_cofactor_bit_binomials import direct_ratio
    from cofactor_zero_queries import factorization, largest_prime, semiprime_from_pairs, zero_from_complement


@cache
def _ratio(n, j):
    return direct_ratio(n >> j, n)


def direct_zero(n, mask):
    residue = 1
    while mask:
        bit = mask & -mask
        residue = residue * _ratio(n, bit.bit_length()-1) % n
        mask -= bit
    return pow(residue, n.bit_length(), n) == 0


def signature(n, p):
    value = 0
    while n >= p:
        n //= p
        value |= n
    return value


def verify():
    semiprimes = 0
    for n in range(1, 2501):
        factors = factorint(n)
        assert factorization(n, direct_zero) == factors
        if n % 2 and n > 1:
            assert largest_prime(n, direct_zero) == max(factors)
            if len(factors) == 2 and sum(factors.values()) == 2:
                calls = set()
                def pair_query(modulus, mask):
                    assert modulus == n and 1 <= mask.bit_count() <= 2
                    calls.add(mask)
                    return direct_zero(n, mask)
                assert semiprime_from_pairs(n, pair_query) in factors
                if max(factors) < 2*min(factors):
                    assert len(calls) <= 2*(n.bit_length()-2)
                semiprimes += 1
    print(f"2500 complete direct factorizations; {semiprimes} semiprimes from single/pair tests")

    queries_checked = 0
    for n in range(3, 192, 2):
        for packed in range(1 << (n.bit_length()-2)):
            mask = packed << 1
            assert zero_from_complement(n, mask) == direct_zero(n, mask)
            queries_checked += 1
    print(f"{queries_checked} complete-domain comparisons for the complement evaluator")

    rng = Random(862)
    primes = [int(nextprime(2**bits+rng.randrange(2**(bits-1)))) for bits in (19, 37, 61, 117)]
    large = {primes[0]: 3, primes[1]: 2, primes[2]: 1, primes[3]: 1}
    n = prod(p**e for p, e in large.items())
    # A specified dense query at a balanced input. These benchmark positions
    # were chosen with the test factors; verification itself is input-only.
    p, q = 162*2**120+1, 231*2**120+1
    dense_n = p*q
    all_bits = (1 << (dense_n.bit_length()-1))-2
    free = (1 << 121) | (1 << 125) | (1 << 127)
    assert p == 1+free
    work = {}
    assert not zero_from_complement(dense_n, all_bits ^ free, stats=work)
    assert work["unselected"] == 3 and work["candidates"] <= 8
    assert zero_from_complement(dense_n, all_bits)
    print(f"{dense_n.bit_length()}-bit dense queries, both confirmed with at most eight masks")

    queries = 0
    def planted(modulus, mask):
        nonlocal queries
        queries += 1
        return all(mask & signature(modulus, p) for p in large if modulus % p == 0)
    assert factorization(n, planted) == large
    print(f"{n.bit_length()}-bit pre-factored control, complete answer from {queries} booleans")

    transcript = []
    def record(modulus, mask):
        answer = direct_zero(modulus, mask)
        transcript.append((mask, answer))
        return answer
    assert largest_prime(3**3*5**2*7, record) == 7
    print(f"N=4725 largest prime 7, {len(transcript)} direct yes/no queries")


if __name__ == "__main__":
    verify()
