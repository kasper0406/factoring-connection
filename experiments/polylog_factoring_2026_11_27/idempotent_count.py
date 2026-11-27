"""Encode CRT prime-power components in one ordinary quadratic box count.

The query uses selector s, free variable y, and r pairs (z_j, w_j). Over
N, impose s(s-1)=0, sy=0, w_j-s*z_j=0, and z_j(w_j-1)=0. Over each
prime-power component q, s=0 gives q choices and s=1 gives A=2**r choices.
The CRT count is therefore P(A)=product(q+A). The same equations can be
enforced by complete cubic additive-character sums over multipliers, with
square-root normalization N**(2+2*r).

All primal coordinates lie in [0,N-1], and each congruence C=0 mod N
has a unique quotient h in [1-N,N-1] when it is true. Each |C|<N**2,
so |C-N*h|<2*N**2 for every candidate. The base B=4*N**2+1 makes one
packed equation sum(B**j*(C_j-N*h_j))=0 equivalent to all equations.

With r=0 the box has four coordinates and a rank-two quadratic part.
When there are exactly two components, P(1)-N-1 is their sum. For any
number of components use r=2*N.bit_length(): every coefficient of P is
at most 2**omega(N)*N < A, and the base-A digits of P(A) are exact.
Construction and decoding take polynomial bit time. Finding the count
is the oracle task. The reference counts below use power-of-N scans.
"""

from argparse import ArgumentParser
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, product
from json import dumps
from math import gcd, isqrt, prod
from sys import set_int_max_str_digits

from sympy import Poly, Symbol


def evaluate(terms, point):
    return sum(coefficient * prod(point[index] for index in indices)
               for coefficient, indices in terms)


@dataclass(frozen=True)
class PackedQuery:
    modulus: int
    bits: int
    base: int
    names: tuple
    bounds: tuple
    congruences: tuple
    terms: tuple

    def lift(self, primal):
        """Return the unique in-box point above these primal coordinates, if any."""
        length = len(self.congruences)
        if len(primal) != length or any(not low <= value <= high
                                        for value, (low, high) in zip(primal, self.bounds)):
            raise ValueError('Expected one in-box value per primal coordinate')
        quotients = []
        for equation in self.congruences:
            quotient, remainder = divmod(evaluate(equation, primal), self.modulus)
            if remainder:
                return None
            quotients.append(quotient)
        return (*primal, *quotients)


def make_query(modulus, bits):
    if modulus <= 1 or bits < 0:
        raise ValueError('Expected modulus greater than one and nonnegative bit count')
    base = 4 * modulus**2 + 1
    names = ['s', 'y']
    congruences = [((1, (0, 0)), (-1, (0,))), ((1, (0, 1)),)]
    for index in range(bits):
        z, w = len(names), len(names) + 1
        names.extend((f'z{index}', f'w{index}'))
        congruences.extend((((1, (w,)), (-1, (0, z))),
                            ((1, (z, w)), (-1, (z,)))))

    length = len(congruences)
    bounds = [(0, modulus - 1)] * length
    bounds.extend([(1 - modulus, modulus - 1)] * length)
    terms = []
    multiplier = 1
    for index, equation in enumerate(congruences):
        terms.extend((multiplier * coefficient, indices) for coefficient, indices in equation)
        terms.append((-multiplier * modulus, (length + index,)))
        names.append(f'h{index}')
        multiplier *= base
    return PackedQuery(modulus, bits, base, tuple(names), tuple(bounds),
                       tuple(congruences), tuple(terms))


def radix_bits(modulus):
    return 2 * modulus.bit_length()


def reference_count(modulus, bits):
    """Scan selector residues without a factorization, for small controls."""
    total = 0
    for s in range(modulus):
        if s * (s - 1) % modulus:
            continue
        count_y = gcd(s, modulus)
        count_bit = 1
        if bits:
            count_bit = sum(z * (s * z - 1) % modulus == 0 for z in range(modulus))
        total += count_y * count_bit**bits
    return total


def decode_two_components(modulus, observation):
    total = round(observation) - modulus - 1
    discriminant = total * total - 4 * modulus
    if total < 1 or discriminant < 0:
        raise ValueError('The count does not give two positive CRT components')
    difference = isqrt(discriminant)
    if difference * difference != discriminant or (total - difference) % 2:
        raise ValueError('The count does not give two integral CRT components')
    components = ((total - difference) // 2, (total + difference) // 2)
    if components[0] <= 1 or gcd(*components) != 1 or prod(components) != modulus:
        raise ValueError('The decoded values are not a two-component CRT decomposition')
    return components


def decode_radix(modulus, observation):
    radix = 1 << radix_bits(modulus)
    value = round(observation)
    digits = []
    while value > 0:
        value, last = divmod(value, radix)
        digits.append(last)
    if len(digits) < 2 or digits[-1] != 1 or digits[0] != modulus:
        raise ValueError('The count does not encode a monic polynomial of the expected form')

    polynomial = Poly.from_list(digits[::-1], Symbol('x'), domain='ZZ')
    roots = polynomial.ground_roots()
    if len(roots) != polynomial.degree() or any(root >= -1 for root in roots):
        raise ValueError('The count does not give distinct negative integral roots')
    components = tuple(sorted(-int(root) for root in roots))
    if prod(components) != modulus or any(gcd(a, b) != 1 for a, b in combinations(components, 2)):
        raise ValueError('The decoded values are not a CRT decomposition')
    return components


def check():
    literal_checks = 0
    for modulus, components, bits in ((3, (3,), 0), (5, (5,), 0), (9, (9,), 0),
                                      (15, (3, 5), 0), (3, (3,), 1)):
        query = make_query(modulus, bits)
        expected = prod(q + (1 << bits) for q in components)
        observed = 0
        for point in product(*(range(lo, hi + 1) for lo, hi in query.bounds)):
            if evaluate(query.terms, point) == 0:
                assert point == query.lift(point[:len(query.congruences)])
                observed += 1
        assert observed == expected
        literal_checks += 1

    scan_checks = 0
    fixtures = ((7,), (3, 5), (9, 5), (5, 7), (9, 7), (3, 5, 7),
                (9, 25, 7), (3, 5, 7, 11), (2, 3, 5), (8, 9, 11))
    for components in fixtures:
        modulus = prod(components)
        for bits in (0, 1, 2, radix_bits(modulus)):
            count = reference_count(modulus, bits)
            assert count == prod(q + (1 << bits) for q in components)
            if bits == 0 and len(components) == 2:
                assert decode_two_components(modulus, count) == tuple(sorted(components))
            if bits == radix_bits(modulus):
                assert decode_radix(modulus, count) == tuple(sorted(components))
            scan_checks += 1

    radix_checks = 0
    for components in ((101, 109), (101**2, 109), (3**90, 5**46),
                       (11, 13, 17, 19, 23, 29, 31, 37),
                       (3**25, 5**7, 7**8, 11**14, 13**5, 17**4, 19, 23)):
        modulus = prod(components)
        radix = 1 << radix_bits(modulus)
        for error in (Fraction(-49, 100), Fraction(49, 100)):
            if len(components) == 2:
                small = prod(q + 1 for q in components) + error
                assert decode_two_components(modulus, small) == tuple(sorted(components))
            large = prod(q + radix for q in components) + error
            assert decode_radix(modulus, large) == tuple(sorted(components))
            radix_checks += 1
    print(f'{literal_checks} literal box counts, {scan_checks} selector scans, '
          f'{radix_checks} synthetic noisy radix decodes passed')


def main():
    # The packed coefficients and counts can legitimately exceed 4300 digits.
    set_int_max_str_digits(0)
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('modulus', nargs='?', type=int)
    parser.add_argument('--two-components', action='store_true', help='Use the four-coordinate r=0 query')
    parser.add_argument('--count', type=Fraction, help='Decode a supplied integer or rational observation')
    parser.add_argument('--reference', action='store_true', help='Enumerate a small count directly')
    args = parser.parse_args()
    if args.modulus is None:
        if args.two_components or args.count is not None or args.reference:
            parser.error('Supply a modulus for this option')
        check()
        return
    if args.modulus <= 1 or (args.reference and args.count is not None):
        parser.error('Expected modulus greater than one and at most one count source')
    bits = 0 if args.two_components else radix_bits(args.modulus)
    observation = args.count
    if args.reference:
        observation = reference_count(args.modulus, bits)
        print(f'Count by selector scan = {observation}')
    if observation is not None:
        decoder = decode_two_components if args.two_components else decode_radix
        try:
            print('CRT components', *decoder(args.modulus, observation))
        except ValueError as error:
            parser.error(str(error))
        return

    query = make_query(args.modulus, bits)
    print(dumps({'coordinates': query.names, 'inclusive_bounds': query.bounds,
                 'equation_terms': query.terms, 'right_hand_side': 0, 'packing_base': query.base}))


if __name__ == '__main__':
    main()
