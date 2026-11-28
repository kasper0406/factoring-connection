"""Finite check that one aligned block counts subsets using roots of x²=1.

Each fixture has a single squarefree modulus for all positive target sums.
The core primes are successive primes above 32*k*L. We scan a prescribed
progression inside dyadic bands for the extra prime. This finite search
gives no complexity bound for constructing it. SymPy's isprime test is
deterministic below 2^64; the combined modulus is unrestricted.
"""

from collections import Counter
from itertools import product
from math import gcd, prod

from sympy import isprime, nextprime


def ceil_div(a, b):
    return -(-a // b)


def crt(residues, moduli):
    modulus = prod(moduli)
    result = 0
    for value, q in zip(residues, moduli, strict=True):
        other = modulus // q
        result += value * other * pow(other, -1, q)
    return result % modulus


def band_prime(core, residue, total):
    maximum_prime = (1 << 64) - 1
    tested = 0
    for h in range(core.bit_length() - 1, (core * maximum_prime).bit_length()):
        power = 1 << h
        lo = ceil_div(power, core)
        hi = min(maximum_prime, power * (4 * total + 1) // (4 * total * core))
        first = residue + max(0, ceil_div(lo - residue, core)) * core
        for candidate in range(first, hi + 1, core):
            tested += 1
            if candidate > 2 and isprime(candidate):
                return candidate, h, tested
    raise AssertionError("Fixture exhausted the deterministic primality range")


def verify(weights):
    total = sum(weights)
    scale = 1 << (16 * (total + 1) - 1).bit_length()
    cores = []
    coefficients = []
    q = 32 * len(weights) * scale
    for weight in weights:
        q = int(nextprime(q))
        k = ceil_div(q * weight, scale)
        assert q < 2 ** 64 and isprime(q) and 0 < k < q and gcd(k, q) == 1
        assert 0 <= scale * k - q * weight < scale
        cores.append(q)
        coefficients.append(k)

    core = prod(cores)
    residues = [pow(k * (core // q), -1, q) for k, q in zip(coefficients, cores)]
    prime, h, tested = band_prime(core, crt(residues, cores), total)
    moduli = cores + [prime]
    n = prod(moduli)
    assert len(set(moduli)) == len(moduli)
    assert (1 << h) < n and 4 * total * n <= (4 * total + 1) * (1 << h)
    assert coefficients == [pow(n // q, -1, q) for q in cores]
    width, remainder = divmod(1 << (h + 1), scale)
    assert not remainder and width > 0 and width & (width - 1) == 0

    block_counts = Counter()
    roots = set()
    for signs in product((-1, 1), repeat=len(moduli)):
        x = crt(signs, moduli)
        assert x not in roots and x * x % n == 1
        roots.add(x)
        target = sum(w for w, sign in zip(weights, signs) if sign == 1)
        if signs[-1] == -1 and target > 0:
            predicted = 2 * sum((n // q) * k for q, k, sign in zip(cores, coefficients, signs)
                                if sign == 1) - 1
            assert x == predicted and x // width == target
        elif signs[-1] == -1:
            assert x == n - 1
        elif target == total:
            assert x == 1
        else:
            assert x > n // 2
        if 1 <= x // width <= total:
            block_counts[x // width] += 1
    assert len(roots) == 2 ** len(moduli)
    wanted = Counter(sum(w * f for w, f in zip(weights, bits))
                     for bits in product((0, 1), repeat=len(weights)))
    for target in range(1, total + 1):
        lo, hi = target * width, (target + 1) * width
        assert 1 < lo < hi <= n // 2 and n - 1 >= hi and lo % width == 0
        assert block_counts[target] == wanted[target] == sum(lo <= x < hi for x in roots)
        half = n // 2
        assert not any(lo + half <= x < hi + half for x in roots)
        first = sum(lo <= x < lo + half for x in roots)
        second = sum(hi <= x < hi + half for x in roots)
        assert first - second == wanted[target]
    print(weights, cores, prime, h, tested, width,
          [block_counts[t] for t in range(1, total + 1)], sep="  ")


if __name__ == "__main__":
    print("weights  core_primes  anchor_prime  band_exponent  candidates_tested  block_width  counts_for_T=1..S")
    for instance in ((1,), (1, 1), (1, 1, 2), (1, 2, 3), (2, 2, 3), (1, 1, 1, 1)):
        verify(instance)
