"""Bach–Miller–Shallit coefficient gcds given a multiple of a hidden norm-one order."""

from collections import Counter
from math import gcd, lcm
from random import Random, SystemRandom

if __package__:
    from .check_translated_batch import power
else:
    from check_translated_batch import power


def scalar_membership_probe(n, h, y, exponent):
    """A field scalar has zero coefficients in every positive degree."""
    value = power(y, exponent, h, n)
    divisor = gcd(n, *value[1:])
    return (divisor, "coefficients") if 1 < divisor < n else None


def exponent_trial(n, degree, exponent, rng):
    a = rng.randrange(1, n)
    divisor = gcd(a, n)
    if divisor > 1:
        return divisor, "scalar unit check"
    twos = (exponent & -exponent).bit_length() - 1
    value = pow(a, exponent >> twos, n)
    for _ in range(twos + 1):
        divisor = gcd(value - 1, n)
        if 1 < divisor < n:
            return divisor, "scalar ladder"
        value = value * value % n

    h = tuple(rng.randrange(n) for _ in range(degree)) + (1,)
    y = tuple(rng.randrange(n) for _ in range(degree))
    return scalar_membership_probe(n, h, y, exponent)


def split_norm_one_exponent(n, degree, exponent, rng=None):
    """Find a proper divisor in expected O(degree) trials under the supplied promise.

    The odd input n has at least two distinct prime divisors. For some prime r | n,
    (r**degree - 1) // (r - 1) divides exponent. Prime powers in n are allowed.
    """
    if n < 6 or degree < 2 or exponent < 1:
        raise ValueError("expected composite n >= 6, degree >= 2, and a positive exponent")
    if n % 2 == 0:
        return 2
    if rng is None:
        rng = SystemRandom()
    while True:
        hit = exponent_trial(n, degree, exponent, rng)
        if hit is not None:
            return hit[0]


def check():
    # Only the control constructs exponents from the hidden input primes.
    fixtures = ((31, 43, 2), (37, 43, 3), (41, 53, 4), (43, 47, 5), (47, 59, 6))
    totals = Counter()
    for p, q, degree in fixtures:
        n = p*q
        order_p = (p**degree - 1) // (p-1)
        order_q = (q**degree - 1) // (q-1)
        modes = (("one", order_p), ("both", lcm(order_p, order_q)),
                 ("both and scalars", lcm(order_p, order_q, p-1, q-1)))
        for mode, exponent in modes:
            rng = Random(99173 + n + exponent)
            factor = split_norm_one_exponent(n, degree, exponent, rng)
            assert factor in (p, q)
            for _ in range(12):
                result = exponent_trial(n, degree, exponent, rng)
                if result:
                    factor, source = result
                    assert factor in (p, q)
                    totals[source] += 1
                else:
                    totals["miss"] += 1
            totals["full factorizations"] += 1
    assert totals["coefficients"] and totals["scalar ladder"]
    # Both cubic norm-one orders divide this odd exponent. The powered element
    # is the scalar 32 modulo 37, and has coefficients (24, 33, 39) modulo 43.
    exponent = lcm(37**2 + 37 + 1, 43**2 + 43 + 1)
    h, y = (1, 5, 0, 1), (1, 1, 0)
    assert power(y, exponent, h, 37*43) == (1142, 592, 555)
    assert scalar_membership_probe(37*43, h, y, exponent) == (37, "coefficients")
    assert split_norm_one_exponent(3*101, 5, (3**5-1)//2, Random(4)) in (3, 101)
    for n in (37*43*47, 37**2*43):
        factor = split_norm_one_exponent(n, 3, (37**3-1)//36, Random(n))
        assert 1 < factor < n and n % factor == 0
        totals["additional composite splits"] += 1
    # A relative norm over a quadratic subfield kills only the odd-degree coefficients.
    relative = power((1, 1), 26, (-2, 0, 0, 0, 1), 35)
    assert relative == (11, 0, 9, 25)
    assert gcd(35, relative[1], relative[3]) == 5
    return dict(totals)


if __name__ == "__main__":
    print(check())
