"""Recover the unique scalar candidate for a target when Y generates the algebra."""

from collections import Counter
from dataclasses import dataclass
from math import gcd

from sympy import Matrix

if __package__:
    from .check_translated_batch import multiply, norm, power
else:
    from check_translated_batch import multiply, norm, power


def evaluate(coefficients, argument, modulus):
    value = 0
    for coefficient in reversed(coefficients):
        value = (value * argument + coefficient) % modulus
    return value


@dataclass(frozen=True)
class BatchWitness:
    modulus: int
    polynomial: tuple
    image: tuple
    determinant: int
    determinant_gcd: int
    inverse_polynomial: tuple | None

    def decode(self, target):
        """Return (candidate, resultant gcd) for an invertible power basis."""
        if self.inverse_polynomial is None:
            raise ValueError("decoding requires a unit power-basis determinant")
        n = self.modulus
        candidate = evaluate(self.inverse_polynomial, target, n)
        divisor = gcd(n, evaluate(self.polynomial, candidate, n),
                      evaluate(self.image, candidate, n) - target)
        return candidate, divisor


def prepare_witness(n, h, image):
    """Use the power basis of Y=image in (Z/n)[X]/H, with low coefficients first.

    A proper determinant_gcd is already a factor. A gcd of n is a singular
    result with no factor furnished by this determinant. A gcd of 1 supplies
    L with L(Y)=X, enabling decode(target) for any scalar target.
    """
    h = tuple(h)
    if n < 2 or len(h) < 2 or h[-1] != 1:
        raise ValueError("expected n >= 2 and a nonconstant monic H")
    value = multiply((1,), tuple(image) or (0,), h, n)
    columns = [multiply((1,), (1,), h, n)]
    for _ in range(1, len(h) - 1):
        columns.append(multiply(columns[-1], value, h, n))
    basis = Matrix(columns).T
    determinant = int(basis.det(method="berkowitz")) % n
    divisor = gcd(n, determinant)
    inverse = None
    if divisor == 1:
        x = Matrix(multiply((0, 1), (1,), h, n))
        numerators = basis.adjugate(method="berkowitz") * x
        denominator_inverse = pow(determinant, -1, n)
        inverse = tuple(int(a) * denominator_inverse % n for a in numerators)
    return BatchWitness(n, h, value, determinant, divisor, inverse)


def check():
    totals = Counter()
    inputs = ((3, 1), (0, 0, 1), (1, 0, 1), (0, -1, 0, 1), (-2, 5, 0, 0, 1))
    for n in (4, 6, 8, 9, 12, 15, 25, 35, 45, 49):
        for h in inputs:
            images = ((0,), (0, 1), (2, -1, 1), (0, 2), power((-1, 1), 5, h, n))
            for y in images:
                witness = prepare_witness(n, h, y)
                if witness.determinant_gcd != 1:
                    mode = "singular modulo N" if witness.determinant_gcd == n else "determinant splits"
                    totals[mode] += 1
                    assert witness.inverse_polynomial is None
                    continue
                totals["unit determinants"] += 1
                # Check L(Y)=X directly in the quotient, including linear H.
                composition = (0,) * (len(h) - 1)
                for coefficient in reversed(witness.inverse_polynomial):
                    composition = multiply(composition, witness.image, h, n)
                    composition = ((composition[0] + coefficient) % n, *composition[1:])
                assert composition == multiply((0, 1), (1,), h, n)
                for target in range(n):
                    candidate, divisor = witness.decode(target)
                    difference = (witness.image[0] - target, *witness.image[1:])
                    assert divisor == gcd(n, norm(difference, h, n))
                    assert {(x % divisor) for x in range(divisor)
                            if evaluate(h, x, divisor) == 0
                            and evaluate(y, x, divisor) == target % divisor} == {candidate % divisor}
                    totals["scalar targets"] += 1

    # At 5, X^4=2 and X^2 lies in the quadratic extension with (X^2)^2=2.
    # X^25-X reduces to 3X, so H has no linear or quadratic factor at 5.
    n, h, y = 35, (-2, 5, 0, 0, 1), (0, 0, 1)
    witness = prepare_witness(n, h, y)
    assert (witness.determinant, witness.determinant_gcd) == (10, 5)
    assert power((0, 1), 25, h, 5) == (0, 4, 0, 0)
    assert all(norm((-u, 0, 1, 0), h, 5) == (u*u-2)**2 % 5 != 0 for u in range(5))
    assert norm((-1, 0, 1, 0), h, n) == 11
    assert gcd(11, n) == 1

    # Distinguish a singular determinant from a failed scalar probe.
    singular = prepare_witness(35, (1, 0, 1), (0,))
    assert singular.determinant_gcd == 35 and singular.inverse_polynomial is None
    return dict(totals)


if __name__ == "__main__":
    print(check())
