"""Find a separating power translate of a squarefree root set, or a proper factor.

The input is a composite integer. Coefficients of monic H run from low to high.
"""

from dataclasses import dataclass
from math import gcd, lcm
from random import Random, SystemRandom

if __package__:
    from .check_batch_witness import BatchWitness, prepare_witness
    from .check_correlated_batch import squarefree_support
    from .check_norm_one_exponent import exponent_trial
    from .check_translated_batch import multiply, power
else:
    from check_batch_witness import BatchWitness, prepare_witness
    from check_correlated_batch import squarefree_support
    from check_norm_one_exponent import exponent_trial
    from check_translated_batch import multiply, power


def tower_power(n, base_poly, upper_poly, element, exponent):
    """Power in B[V]/upper_poly, where B=(Z/n)[U]/base_poly and both are monic."""
    if exponent < 0:
        raise ValueError("expected a nonnegative tower exponent")
    g = len(base_poly) - 1
    k = len(upper_poly) - 1
    zero = (0,) * g
    one = (1,) + zero[1:]

    def product(left, right):
        coefficients = [[0] * g for _ in range(2 * k - 1)]
        for i, a in enumerate(left):
            for j, b in enumerate(right):
                for pos, term in enumerate(multiply(a, b, base_poly, n)):
                    coefficients[i+j][pos] += term
        for i in range(2*k-2, k-1, -1):
            for j in range(k):
                reduction = multiply(coefficients[i], upper_poly[j], base_poly, n)
                for pos, term in enumerate(reduction):
                    coefficients[i-k+j][pos] -= term
        return tuple(tuple(v % n for v in row) for row in coefficients[:k])

    value = (one,) + (zero,) * (k - 1)
    while exponent:
        if exponent & 1:
            value = product(value, element)
        exponent >>= 1
        if exponent:
            element = product(element, element)
    return value


def relative_trial(n, g, k, exponent, rng):
    """One relative-norm-order trial. The subfield has degree g, the tower degree g*k."""
    base_poly = tuple(rng.randrange(n) for _ in range(g)) + (1,)

    def base_element():
        return tuple(rng.randrange(n) for _ in range(g))

    upper_poly = tuple(base_element() for _ in range(k)) + ((1,) + (0,) * (g-1),)
    y = tuple(base_element() for _ in range(k))
    value = tower_power(n, base_poly, upper_poly, y, exponent)
    divisor = gcd(n, *(v for row in value[1:] for v in row))
    return divisor if 1 < divisor < n else None


@dataclass(frozen=True)
class ShiftResult:
    factor: int | None = None
    shift: int | None = None
    witness: BatchWitness | None = None


def factor_or_shift(n, h, m, rng=None):
    """For composite n return a factor or a power-basis decoder on the radical of h."""
    from sympy import integer_nthroot

    if n < 4 or m < 1 or len(h) < 2 or h[-1] != 1:
        raise ValueError("expected composite n, positive m, and nonconstant monic h")
    if n % 2 == 0:
        return ShiftResult(factor=2)
    divisor = gcd(n, m)
    if 1 < divisor < n:
        return ShiftResult(factor=divisor)

    for exponent in range(2, n.bit_length()):
        root, exact = integer_nthroot(n, exponent)
        if exact:
            return ShiftResult(factor=root)

    divisor, support = squarefree_support(n, h)
    if divisor is not None:
        return ShiftResult(factor=divisor)
    degree = len(support) - 1
    if rng is None:
        rng = SystemRandom()
    limit = 4 * n.bit_length() * degree**2
    small_bound = (4 * n.bit_length() * degree**3)**2
    next_odd = 3
    while True:
        c = rng.randrange(n)
        image = power((c, -1), m, support, n)
        witness = prepare_witness(n, support, image)
        if witness.determinant_gcd == 1:
            return ShiftResult(shift=c, witness=witness)
        if witness.determinant_gcd != n:
            return ShiftResult(factor=witness.determinant_gcd)
        # Interleave the small-prime alternative with later shift trials.
        if next_odd <= small_bound and next_odd < n:
            if n % next_odd == 0:
                return ShiftResult(factor=next_odd)
            next_odd += 2
        # A missing cofactor up to limit is sampled with probability at least 1/limit.
        # Keeping the exponent short avoids forming the lcm of the whole range.
        order_multiple = m * rng.randrange(1, limit + 1)
        for d in range(2, degree + 1):
            result = exponent_trial(n, d, order_multiple, rng)
            if result is not None:
                return ShiftResult(factor=result[0])
            for g in range(2, d):
                if d % g == 0:
                    divisor = relative_trial(n, g, d // g, order_multiple, rng)
                    if divisor is not None:
                        return ShiftResult(factor=divisor)


def check():
    n = 35
    value = tower_power(n, (-2, 0, 1), ((0,-1), (0,0), (1,0)),
                        ((1,0), (1,0)), 26)
    assert value == ((11,9), (0,25))
    assert gcd(n, *value[1]) == 5
    splits = decoders = 0
    cases = ((35, (-2,5,0,0,1), 2), (77, (-2,0,1), 492),
             (143, (-1,0,1), 5), (149*211, (0,1,2,2,1), 149*211),
             (225, (1,2,1), 37), (45, (0,0,1), 5), (105, (-2,0,1), 105))
    for n, h, m in cases:
        for seed in range(4):
            result = factor_or_shift(n, h, m, Random(seed))
            if result.factor is not None:
                assert 1 < result.factor < n and n % result.factor == 0
                splits += 1
            else:
                assert result.witness.determinant_gcd == 1
                assert result.witness.image == power((result.shift,-1), m, result.witness.polynomial,n)
                decoders += 1
    # An inert quadratic at both primes forces collisions at every shift.
    result = factor_or_shift(77, (1,0,1), lcm(7**2-1,11**2-1), Random(9))
    assert result.factor in (7,11)
    return splits, decoders


if __name__ == "__main__":
    print("divisors and decoder controls", check())
