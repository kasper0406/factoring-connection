"""Check the collision bound for nonuniform independent CRT marginals."""

from fractions import Fraction
from random import Random

from check_power_sampling import square_splits


def check(p, q, rng):
    n = p * q
    inverse = pow(p, -1, q)
    successes = {
        (i, j)
        for i in range(1, p)
        for j in range(1, q)
        if square_splits(i + p * ((j - i) * inverse % q), n)
    }
    for sample in range(400):
        left = [rng.randrange(13) for _ in range(p - 1)]
        right = [rng.randrange(13) for _ in range(q - 1)]
        if sample % 3:
            left[rng.randrange(p - 1)] += 37
        if sample % 3 == 2:
            right[rng.randrange(q - 1)] += 71
        if not sum(left) or not sum(right):
            continue
        alpha = [Fraction(x, sum(left)) for x in left]
        beta = [Fraction(x, sum(right)) for x in right]
        s = sum(alpha[i - 1] * beta[j - 1] for i, j in successes)
        a = sum(x * x for x in alpha)
        b = sum(y * y for y in beta)
        collision = a + b - 2 * a * b
        assert s <= max(alpha) + max(beta)
        if a > Fraction(1, 2) and b > Fraction(1, 2):
            mode = (alpha.index(max(alpha)) + 1, beta.index(max(beta)) + 1)
            if mode in successes:
                assert s > Fraction(1, 4)
                continue
        assert 4 * collision >= s * s


if __name__ == "__main__":
    rng = Random(261127)
    for p, q in [(3, 5), (5, 7), (11, 13), (17, 23), (11, 37)]:
        check(p, q, rng)
    print("Checked 2,000 nonuniform marginal pairs across five moduli.")
