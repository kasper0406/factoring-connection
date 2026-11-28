"""Check collisions, direct gcd hits, and prime support for positive binary substrings."""

from collections import Counter
from math import gcd, prod

from sympy import divisor_count, primerange

if __package__:
    from .check_digit_order_support import order_powers
    from .check_root_scale_order_hints import prime_pairs
else:
    from check_digit_order_support import order_powers
    from check_root_scale_order_hints import prime_pairs


def proper_words(n):
    """Return distinct positive words. The complete input and empty words are excluded."""
    words = set()
    bits = n.bit_length()
    for low in range(bits):
        suffix = n >> low
        for width in range(1, bits-low+1):
            if low or width < bits:
                word = suffix % (1 << width)
                if word:
                    words.add(word)
    return words


def check_fibres(radius):
    odds = range(radius | 1, 2*radius, 2)
    partners = range(radius, 2*radius)
    divisors = [0] + [int(divisor_count(i)) for i in range(1, 2*radius*radius+1)]
    prime_set = set(map(int, primerange(radius, 2*radius)))
    cases = direct_hits = off_diagonal = 0
    for top in range(1, (4*radius**2-1).bit_length()+1):
        h = 1 << top
        for low in range(top):
            b = 1 << low
            c = h//b
            distances = ((p*v) % h for p in odds for v in range(1, radius))
            close_pairs = sum(min(offset, h-offset) < b for offset in distances)
            divisor_majorant = sum(number for n, number in enumerate(divisors) if
                                   min(n % h, (-n) % h) < b)
            assert close_pairs <= divisor_majorant
            collisions = 0
            for p in odds:
                counts = Counter((p*q//b) % c for q in partners)
                collisions += sum(count*(count-1) for count in counts.values())
            assert collisions <= 2*radius*close_pairs
            off_diagonal += collisions

            # Classify each direct hit by its nonzero outside bits J=N-Bz.
            complements = Counter()
            for p in prime_set:
                for q in prime_set-{p}:
                    n = p*q
                    z = (n//b) % c
                    if z and z % p == 0 and n != b*z:
                        j = n-b*z
                        assert j > 0 and j % p == 0 and j % h < b
                        assert q == j//p + b*(z//p)
                        complements[j, p] += 1
            assert all(count <= (c-1)//radius for count in complements.values())
            direct_hits += sum(complements.values())
            cases += 1
    return cases, off_diagonal, direct_hits


def check_shifted_gcd_witnesses():
    hits = whole = 0
    primes = tuple(map(int, primerange(16, 32)))
    for p in primes:
        for q in primes:
            if p == q:
                continue
            n = p*q
            for a, b in ((1, -3), (2, 0), (2, 1), (3, 5), (6, -7)):
                value = a*n+b
                for low in range(value.bit_length()):
                    scale = 1 << low
                    for width in range(1, value.bit_length()-low+1):
                        word = (value >> low) % (1 << width)
                        if not word or word % p:
                            continue
                        outside = value-scale*word
                        if outside == b:
                            assert word % n == 0
                            whole += 1
                        else:
                            assert (outside-b) % p == 0
                            assert a*q == scale*(word//p)+(outside-b)//p
                            hits += word % n != 0
    assert hits and whole
    return hits, whole


def benchmark(bits, count=128):
    totals = [0]*4
    for p, q in prime_pairs(bits, count, seed_offset=4_000_037):
        n = p*q
        words = proper_words(n)
        assert n not in words and 0 not in words
        if bits == 16:
            binary = bin(n)[2:]
            substrings = {binary[i:j] for i in range(len(binary))
                          for j in range(i+1, len(binary)+1)}
            substrings.discard(binary)
            assert words == {int(s, 2) for s in substrings} - {0}
        best = [1, 1, 1]
        for r in (p, q):
            powers = order_powers(r)
            supported = {ell: val for ell, val in powers.items()
                         if any(z % ell == 0 for z in words)}
            best[1] = max(best[1], prod(sorted(supported.values(), reverse=True)[:2]))
            best[2] = max(best[2], prod(supported.values()))
            for z in words:
                best[0] = max(best[0], prod(v for ell, v in supported.items() if z % ell == 0))
        for i, mass in enumerate(best):
            totals[i] += mass**4 >= n
        totals[3] += any(1 < gcd(n, z) < n for z in words)
    return totals


if __name__ == '__main__':
    for r in (16, 32, 64):
        print("radius, blocks, off-diagonal collisions, direct hits", r, *check_fibres(r))
    print("shifted direct gcd witnesses and whole-modulus cases", *check_shifted_gcd_witnesses())
    print("prime bits | any word | best two pooled bases | all pooled bases | any direct gcd")
    for bits in (16, 32, 48, 64):
        print(bits, *benchmark(bits), flush=True)
