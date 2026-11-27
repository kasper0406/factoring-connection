"""Check one-bit, small-selector count queries on synthetic CRT components."""

from collections import Counter
from itertools import product
from math import gcd, prod

from sympy import nextprime

from idempotent_count import PackedQuery, decode_radix, evaluate, make_query, radix_bits
from stationary_dual import prime_power_base


def selector_query(n, threshold):
    """Construct a packed query of count product(q + threshold)."""
    if threshold < 1:
        raise ValueError('Expected positive threshold')
    bits = (threshold - 1).bit_length()
    source = make_query(n, bits)
    if threshold == 1 << bits:
        return source

    length = len(source.congruences)
    names = list(source.names[:length])
    congruences = list(source.congruences)
    for index in range(bits):
        z = 2 + 2 * index
        state = length + index
        previous = state - 1
        names.append(f'less{index}')
        equation = [(1, (state,))]
        if (threshold >> index) & 1:
            equation.extend(((-1, ()), (1, (z,))))
            if index:
                equation.append((-1, (z, previous)))
        elif index:
            equation.extend(((-1, (previous,)), (1, (z, previous))))
        congruences.append(tuple(equation))

    length = len(congruences)
    bounds = [(0, n - 1)] * (length - 1) + [(1, 1)] + [(-4*n, 4*n)] * length
    base = 16 * n * n + 1
    terms = []
    for index, equation in enumerate(congruences):
        multiplier = base**index
        terms.extend((multiplier * coefficient, indices) for coefficient, indices in equation)
        terms.append((-n * multiplier, (length + index,)))
        names.append(f'h{index}')
    return PackedQuery(n, bits, base, tuple(names), tuple(bounds), tuple(congruences), tuple(terms))


def explicit_selector_count(n, threshold):
    """Verify the Boolean state equations by scanning small primal controls."""
    query = selector_query(n, threshold)
    bits = query.bits
    total = 0
    for s in range(n):
        if s * (s - 1) % n:
            continue
        for values in product(range(n), repeat=bits):
            w = [s * z % n for z in values]
            if any(z * (v - 1) % n for z, v in zip(values, w)):
                continue
            z_coordinates = [item for row in zip(values, w) for item in row]
            states = []
            if threshold != 1 << bits:
                previous = 0
                for index, z in enumerate(values):
                    previous = ((1 - z + z * previous) if (threshold >> index) & 1
                                else previous * (1 - z)) % n
                    states.append(previous)
                if states[-1] != 1:
                    continue
            for y in range(0, n, n // gcd(s, n)):
                point = query.lift((s, y, *z_coordinates, *states))
                assert point is not None and evaluate(query.terms, point) == 0
                assert all(lo <= value <= hi for value, (lo, hi) in zip(point, query.bounds))
                total += 1
    return total


def polynomial_from_roots(residues, prime):
    coefficients = [1]
    for component in residues:
        updated = [0] * (len(coefficients) + 1)
        for index, value in enumerate(coefficients):
            updated[index] = (updated[index] + component * value) % prime
            updated[index + 1] = (updated[index + 1] + value) % prime
        coefficients = updated
    return coefficients


def decode_zero_bits(n, predicate):
    """Recover components from bits saying whether each box count is zero mod a prime."""
    nbits = n.bit_length()
    bits = radix_bits(n)
    radix = 1 << bits
    degree = 0
    combined = 1
    coefficients = [0]
    prime = 1
    work = Counter()
    # At most nbits^3 primes can divide the component-difference product.
    for _ in range(nbits**3 + 2 * nbits + 1):
        prime = int(nextprime(prime))
        residues = [(-threshold) % prime for threshold in range(1, prime + 1)
                    if predicate(threshold, prime)]
        work['primes'] += 1
        work['bit_queries'] += prime
        work['largest_prime'] = prime
        if len(residues) < degree:
            continue
        if len(residues) > degree:
            degree = len(residues)
            combined = 1
            coefficients = [0] * (degree + 1)
        local = polynomial_from_roots(residues, prime)
        inverse = pow(combined, -1, prime)
        coefficients = [value + combined * ((wanted - value) * inverse % prime)
                        for wanted, value in zip(local, coefficients)]
        combined *= prime
        if combined <= (n << nbits) or coefficients[0] != n or max(coefficients) >= radix:
            continue
        observation = sum(value << (index * bits)
                          for index, value in enumerate(coefficients))
        try:
            components = decode_radix(n, observation)
        except ValueError:
            continue
        if all(prime_power_base(value) is not None for value in components):
            work['degree'] = degree
            return components, work
    raise ValueError('No verified component decomposition from the supplied bits')


if __name__ == '__main__':
    checked = 0
    for components in ((3,), (5,), (2, 3), (3, 5), (9, 5)):
        n = prod(components)
        for threshold in (1, 2, 3, 5, 6, 7):
            assert explicit_selector_count(n, threshold) == prod(q + threshold for q in components)
            checked += 1
    print(f'{checked} explicit bounded selector counts checked')

    for components in ((3, 5), (9, 5), (9, 25, 7), (503, 937),
                       (11, 13, 17, 19, 23, 29, 31, 37),
                       (3**25, 5**7, 7**8, 11**14, 13**5, 17**4, 19, 23)):
        n = prod(components)

        def synthetic_predicate(threshold, prime):
            return any((threshold + component) % prime == 0 for component in components)

        found, work = decode_zero_bits(n, synthetic_predicate)
        assert found == tuple(sorted(components))
        print(f'N has {n.bit_length()} bits, {len(found)} components, {dict(work)}')
