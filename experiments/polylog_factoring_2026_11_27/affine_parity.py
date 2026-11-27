"""Recover affine-isolated mixed roots using unweighted disk parity."""

from math import gcd, isqrt
from random import Random

from affine_isolation import make_batch, query_form
from fixed_window_split import reference_count


def auxiliary_primes(n):
    primes, product = [], 1
    candidate = 2
    while product < n:
        if all(candidate % d for d in range(2, isqrt(candidate)+1)):
            if gcd(candidate, n) == 1:
                primes.append(candidate)
                product *= candidate
        candidate += 1
    return primes


def reference_form(n):
    (leading, _, _), cutoff = query_form(n, 1, 0, 1)
    height = n*cutoff
    last = height + n+1 + (1-height) % (n+1)
    return (leading, 0, last), cutoff


def progression_queries(n, scale, offset, width, primes):
    for prime in primes:
        inverse = pow(prime, -1, n)
        group = []
        for r in range(min(prime, width)):
            size = (width-1-r)//prime+1
            a, b = scale*inverse % n, (offset-r)*inverse % n
            form, cutoff = query_form(n, a, b, size)
            correction = sum((b + sign*a) % n < size for sign in (1, -1))
            if size % 2 == 0:
                endpoint = (size-b)*pow(a, -1, n) % n
                correction += int((endpoint**2-1) % n == 0)
            group.append((r, form, cutoff, correction))
        yield prime, group


def decode(n, scale, offset, width, groups, disk_count, baseline):
    coordinate, modulus = 0, 1
    for prime, group in groups:
        odd = []
        for r, form, cutoff, correction in group:
            residue = (disk_count(n, form, cutoff) - baseline) % 4
            assert residue in (0, 2)
            if (residue//2-correction) % 2:
                odd.append(r)
        if len(odd) != 1:
            return None
        coordinate += modulus * ((odd[0]-coordinate) * pow(modulus, -1, prime) % prime)
        modulus *= prime
    if coordinate >= width:
        return None
    root = (coordinate-offset) * pow(scale, -1, n) % n
    if root in (1, n-1) or (root**2-1) % n:
        return None
    divisor = gcd(root-1, n)
    assert 1 < divisor < n
    return divisor


def find_factor(n, rng, disk_count, repeats=None):
    if repeats is None:
        repeats = min(16*(n.bit_length()+1), 65_536)
    divisor, batches = make_batch(n, rng, repeats)
    if divisor is not None:
        return divisor
    primes = auxiliary_primes(n)
    form, cutoff = reference_form(n)
    baseline = disk_count(n, form, cutoff)
    for scale, offset, width in batches:
        groups = progression_queries(n, scale, offset, width, primes)
        divisor = decode(n, scale, offset, width, groups, disk_count, baseline)
        if divisor is not None:
            return divisor
    return None


if __name__ == "__main__":
    rng = Random(261127)
    trials = singleton_trials = queries = 0

    def counted(n, form, cutoff):
        global queries
        queries += 1
        return reference_count(n, form, cutoff)

    for n in (15, 21, 45, 105, 315, 385, 945, 1155, 1575, 3465):
        roots = [r for r in range(2, n-1) if (r*r-1) % n == 0]
        primes = auxiliary_primes(n)
        reference, cutoff = reference_form(n)
        baseline = counted(n, reference, cutoff)
        widths = [1 << j for j in range(n.bit_length()) if (1 << j) <= n//4]
        for width in widths:
            for _ in range(30):
                scale = rng.randrange(1, n)
                while gcd(scale, n) != 1:
                    scale = rng.randrange(1, n)
                offset = rng.randrange(n)
                groups = list(progression_queries(n, scale, offset, width, primes))
                assert all(cap == cutoff for _, group in groups for _, _, cap, _ in group)
                divisor = decode(n, scale, offset, width, groups, counted, baseline)
                if sum((scale*r+offset) % n < width for r in roots) == 1:
                    assert divisor and 1 < divisor < n and n % divisor == 0
                    singleton_trials += 1
                trials += 1
    print(f"{trials} trials, {queries} signed counts, {singleton_trials} singleton decodes")
