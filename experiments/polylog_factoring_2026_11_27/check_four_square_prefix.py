"""Verify scalar-prefix approximations against sphere enumeration and exact divisor sums."""

from fractions import Fraction
from math import isqrt, prod

if __package__:
    from .check_four_square_incidence import prime_powers, quadruples
    from .four_square_prefix import (
        approximate_closed_prefix, approximate_factored_prefix, factor_from_open_prefix,
    )
else:
    from check_four_square_incidence import prime_powers, quadruples
    from four_square_prefix import (
        approximate_closed_prefix, approximate_factored_prefix, factor_from_open_prefix,
    )


def divisor_weights(factors):
    terms = [(1, 1)]
    for p, exponent in factors:
        powers = [(1, 1)] + [
            (p**j, p**(2*j) + p**(2*j - 1)) for j in range(1, exponent + 1)
        ]
        terms = [(d * v, w * local) for d, w in terms for v, local in powers]
    return terms


def exact_count(terms, cutoff, include_zero):
    return 8 * sum(w * (cutoff // d + include_zero) for d, w in terms)


def verify_lower(lower, exact, epsilon):
    assert lower <= exact <= (1 + epsilon) * lower, (lower, exact, epsilon)


def check_spheres():
    points = cutoffs = estimates = 0
    for n in range(1, 46, 2):
        hist = [0] * n
        for a, b, c, d in quadruples(n * n):
            points += 1
            for t in range(n):
                if t * a % n == t * b % n == t * c % n == t * d % n == 0:
                    hist[t] += 1
        factors = prime_powers(n)
        terms = divisor_weights(factors)
        total = 0
        for B, count in enumerate(hist):
            total += count
            opened = total - hist[0]
            assert total == exact_count(terms, B, True)
            assert opened == exact_count(terms, B, False)
            for epsilon in (Fraction(1, 2), Fraction(1, 9)):
                verify_lower(approximate_closed_prefix(n, B, epsilon), total, epsilon)
                for zero, exact in ((True, total), (False, opened)):
                    lower = approximate_factored_prefix(iter(factors), B, epsilon, zero)
                    verify_lower(lower, exact, epsilon)
                estimates += 3
            cutoffs += 1
    return points, cutoffs, estimates


def check_large_estimates():
    epsilon = Fraction(1, 2)
    smooth = [(p, 1) for p in (3, 5, 7, 11, 13, 17, 19)]
    n = prod(p for p, _ in smooth)
    terms = divisor_weights(smooth)
    B = n // 2
    for zero in (False, True):
        exact = exact_count(terms, B, zero)
        trimmed = approximate_factored_prefix(dict(smooth), B, epsilon, zero)
        verify_lower(trimmed, exact, epsilon)
        assert trimmed < exact

    # The cutoff on trial division strips 3² but leaves both large prime
    # bases in each fixture, including the repeated large prime in the second.
    stripped_examples = []
    for factors in ([(3, 2), (1013, 1), (1019, 1)],
                    [(3, 2), (2**31 - 1, 2), (2**61 - 1, 1)]):
        n = prod(p**j for p, j in factors)
        terms = divisor_weights(factors)
        for B in (0, 1, n // 3 + 1, n - 1):
            exact = exact_count(terms, B, True)
            stripped = approximate_closed_prefix(n, B, epsilon)
            verify_lower(stripped, exact, epsilon)
            assert stripped < exact
            stripped_examples.append((n, B, stripped, exact))
    return trimmed, stripped_examples


def check_decoder():
    cases = [prime_powers(n) for n in range(9, 226, 2)
             if any(n % d == 0 for d in range(3, isqrt(n) + 1, 2))]
    cases += [[(3, 6)], [(1013, 4)], [(97, 1), (2**61 - 1, 1)],
              [(101, 2), (103, 1), (107, 1)], [(2**61 - 1, 2)]]
    calls = 0
    for factors in cases:
        n = prod(p**j for p, j in factors)
        least = min(p for p, _ in factors)
        terms = divisor_weights(factors)
        for sign in (-1, 0, 1):
            visited = []

            def noisy(input_n, B, epsilon):
                assert input_n == n and epsilon == Fraction(1, 4)
                visited.append(B)
                # The sign-zero oracle chooses the error direction closest to
                # the decision threshold at each B, including the first jump.
                direction = sign or (1 if B < least else -1)
                return exact_count(terms, B, False) * (1 + direction * epsilon)

            assert factor_from_open_prefix(n, noisy) == least
            assert len(visited) <= (isqrt(n) - 1).bit_length()
            calls += len(visited)

    factors = [(911, 2), (1009, 1)]
    n = prod(p**j for p, j in factors)
    assert factor_from_open_prefix(
        n, lambda _, B, eps: approximate_factored_prefix(factors, B, eps, False)
    ) == 911

    for prime in (1, 3, 17, 101):
        try:
            factor_from_open_prefix(prime, lambda _, B, eps: 8 * B)
        except ValueError:
            continue
        raise AssertionError(f"decoder returned a factor on {prime}")
    return len(cases), 3 * len(cases), calls


if __name__ == "__main__":
    points, cutoffs, estimates = check_spheres()
    trimmed, stripped_examples = check_large_estimates()
    composites, decodings, calls = check_decoder()
    print(f"{points} sphere points, {cutoffs} cutoffs, {estimates} estimates")
    print(f"trimmed closed count {trimmed}; {len(stripped_examples)} rough-cofactor estimates")
    print(f"{composites} odd composites, {decodings} noisy decodings, {calls} oracle calls")
