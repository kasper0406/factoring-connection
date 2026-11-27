"""Decode CRT components from small-prime residues of idempotent counts.

The oracle supplies P(2**r) modulo an odd prime ell, where
P(X)=product(X+q) over the complete prime-power components q of N.
With K=floor(log2(N)), the first ell has K distinct powers of two. We
interpolate (P(X)-N)/X at K points. Its leading coefficient is 1, so its
degree modulo ell reveals k=number of components without cancellation.
Then only the k-1 unknown coefficients of (P(X)-N-X**k)/X are needed,
using k-1 oracle calls per additional prime.

At most K(K+1)/2 primes divide product(2**j-1 for 1<=j<=K). All other
odd primes have enough distinct powers. Each coefficient is at most
2**k*N, so at most 2K+3 accepted primes suffice for CRT. Searching that
many initial primes needs only a polynomial range in K, and the answers
have O(log K) bits. The largest query uses r=K-1 and dimension 4K.
The CLI oracle is synthetic and uses supplied prime-power components.
"""

from argparse import ArgumentParser
from collections import Counter
from itertools import combinations, count
from math import gcd, isqrt, prod

from idempotent_count import decode_radix, radix_bits, reference_count


def evaluate_polynomial(coefficients, x, modulus):
    value = 0
    for coefficient in reversed(coefficients):
        value = (value * x + coefficient) % modulus
    return value


def interpolate(points, values, prime):
    """Return ascending coefficients at distinct points over F_prime."""
    result = [0] * len(points)
    basis = [1]
    for x, y in zip(points, values):
        residual = y - evaluate_polynomial(result, x, prime)
        inverse = pow(evaluate_polynomial(basis, x, prime), -1, prime)
        correction = residual * inverse % prime
        for index, coefficient in enumerate(basis):
            result[index] = (result[index] + correction * coefficient) % prime
        # Multiply the next interpolation basis by (X-x).
        basis = [(-x*basis[0]) % prime] + [
            (left - x*right) % prime for left, right in zip(basis, basis[1:])
        ] + [basis[-1]]
    return result


def decode_residues(modulus, oracle):
    """Recover prime-power components and charge every prime and oracle probe."""
    if modulus <= 1:
        raise ValueError('Expected modulus greater than one')
    max_degree = modulus.bit_length() - 1
    degree = None
    coefficients = []
    combined_modulus = 1
    work = Counter()
    for prime in count(3, 2):
        if any(prime % divisor == 0 for divisor in range(3, isqrt(prime) + 1, 2)):
            continue
        work['prime_candidates'] += 1
        size = max_degree if degree is None else degree - 1
        points = [pow(2, exponent, prime) for exponent in range(size)]
        if len(set(points)) != size:
            work['collision_primes'] += 1
            continue

        values = []
        for exponent, point in enumerate(points):
            residual = oracle(exponent, prime) - modulus
            if degree is not None:
                residual -= pow(point, degree, prime)
            values.append(residual * pow(point, -1, prime) % prime)
        work['oracle_calls'] += size
        work['accepted_primes'] += 1
        work['largest_answer_prime'] = prime
        work['largest_selector_bits'] = max(work['largest_selector_bits'], size - 1)
        work['largest_query_dimension'] = 4 * work['largest_selector_bits'] + 4
        local = interpolate(points, values, prime)

        if degree is None:
            while local and local[-1] == 0:
                local.pop()
            if not local or local[-1] != 1:
                raise ValueError('Initial residues did not encode a monic product')
            degree = len(local)
            work['degree_discovery_calls'] = size
            work['decoded_degree'] = degree
            coefficients = [0] * (degree - 1)
            local.pop()  # The leading coefficient is public after degree discovery.
            if degree == 1:
                return (modulus,), work

        inverse = pow(combined_modulus, -1, prime)
        for index, target in enumerate(local):
            correction = (target - coefficients[index]) * inverse % prime
            coefficients[index] += correction * combined_modulus
        combined_modulus *= prime
        if combined_modulus > (1 << degree) * modulus:
            break

    # Recreate the same radix integer accepted by the exact-count decoder.
    radix = 1 << radix_bits(modulus)
    integer = 1
    for coefficient in reversed(coefficients):
        integer = integer * radix + coefficient
    integer = integer * radix + modulus
    work['crt_modulus_bits'] = combined_modulus.bit_length()
    return decode_radix(modulus, integer), work


def synthetic_oracle(components):
    """The factors of these controls are known to the oracle, never the decoder."""
    def evaluate(selector_bits, prime):
        a = pow(2, selector_bits, prime)
        value = 1
        for component in components:
            value = value * (component + a) % prime
        return value
    return evaluate


def check():
    small = ((3,), (3, 5), (9, 5), (7, 9, 25), (3, 5, 7, 11), (8, 9, 11))
    large = ((101, 109), (101**2, 109), (3**90, 5**46),
             (11, 13, 17, 19, 23, 29, 31, 37),
             (3**25, 5**7, 7**8, 11**14, 13**5, 17**4, 19, 23))
    exact_scans = 0
    total_work = Counter()
    for components in small + large:
        modulus = prod(components)
        oracle = synthetic_oracle(components)
        scanned = {}

        def checked_oracle(bits, prime):
            answer = oracle(bits, prime)
            if components in small:
                if bits not in scanned:
                    scanned[bits] = reference_count(modulus, bits)
                assert answer == scanned[bits] % prime
            return answer

        recovered, work = decode_residues(modulus, checked_oracle)
        assert recovered == tuple(sorted(components))
        assert work['oracle_calls'] == (work['degree_discovery_calls']
                                        + (work['accepted_primes'] - 1) * (len(components) - 1))
        exact_scans += len(scanned)
        total_work.update({key: work[key] for key in ('oracle_calls', 'accepted_primes')})
    print(f'{len(small + large)} modular decodes, {exact_scans} selector scans, '
          f'{total_work["oracle_calls"]} oracle calls across {total_work["accepted_primes"]} primes')


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('components', nargs='*', type=int,
                        help='Complete prime-power components for a synthetic run')
    args = parser.parse_args()
    if not args.components:
        check()
        return
    if (any(q < 2 for q in args.components)
            or any(gcd(a, b) != 1 for a, b in combinations(args.components, 2))):
        parser.error('Supply pairwise coprime components greater than one')
    modulus = prod(args.components)
    recovered, work = decode_residues(modulus, synthetic_oracle(args.components))
    print(f'N={modulus}; recovered components {recovered}; work {dict(work)}')


if __name__ == '__main__':
    main()
