"""Native N-th matrix-power Jacobians and checks of the nilpotent precision threshold."""

from math import gcd, isqrt
from random import Random

from sympy import GF, Matrix, diag, eye, zeros
from sympy.polys.matrices import DomainMatrix


def matrix_power(matrix, exponent, modulus):
    result = eye(matrix.rows)
    base = matrix % modulus
    while exponent:
        if exponent & 1:
            result = (result * base) % modulus
        exponent >>= 1
        if exponent:
            base = (base * base) % modulus
    return result


def power_jacobian(n, matrix):
    """Jacobian of A -> A^n modulo n, in row-major matrix-entry coordinates."""
    d = matrix.rows
    columns = []
    for i in range(d*d):
        direction = zeros(d)
        direction[i] = 1
        block = matrix.row_join(direction).col_join(zeros(d).row_join(matrix))
        derivative = matrix_power(block, n, n)[:d, d:]
        columns.append(derivative.reshape(d*d, 1))
    return Matrix.hstack(*columns)


def shift_factor(n):
    """Evaluate (1+z)^n through degree floor(sqrt(n)), then gcd its positive coefficients.

    For n=p*q with distinct primes p<q this yields q. The explicit table
    has floor(sqrt(n))+1 entries, so this has power-size cost.
    """
    degree = isqrt(n)

    def multiply(left, right):
        result = [0] * min(degree + 1, len(left) + len(right) - 1)
        for i, a in enumerate(left):
            if a:
                for j in range(min(len(right), len(result)-i)):
                    result[i+j] = (result[i+j] + a*right[j]) % n
        return result

    base = [1, 1]
    result = [1]
    exponent = n
    while exponent:
        if exponent & 1:
            result = multiply(result, base)
        exponent >>= 1
        if exponent:
            base = multiply(base, base)
    return gcd(n, *result[1:])


def rank(matrix, prime):
    return DomainMatrix.from_Matrix(matrix).convert_to(GF(prime)).rank()


def commutator_matrix(matrix):
    d = matrix.rows
    columns = []
    for i in range(d*d):
        direction = zeros(d)
        direction[i] = 1
        columns.append((direction*matrix - matrix*direction).reshape(d*d, 1))
    return Matrix.hstack(*columns)


def check():
    a = Matrix([[0, -1], [1, 3]])
    jacobian = power_jacobian(91, a)
    assert matrix_power(a, 91, 91) == Matrix([[25, 13], [78, 77]])
    assert (rank(jacobian, 7), rank(jacobian, 13)) == (2, 0)
    assert gcd(91, *map(int, jacobian)) == 13

    rng = Random(91)
    tested = 0
    for p, q in ((7, 13), (11, 17), (19, 29)):
        for d in (2, 3):
            jordan = eye(d)
            for i in range(d-1):
                jordan[i, i+1] = 1
            examples = [zeros(d), eye(d), jordan, jordan-eye(d), diag(0, *([1]*(d-1)))]
            examples += [Matrix(d, d, [rng.randrange(p*q) for _ in range(d*d)]) for _ in range(4)]
            for a in examples:
                differential = power_jacobian(p*q, a)
                orbit = commutator_matrix(matrix_power(a, p*q, p*q))
                for r in (p, q):
                    expected = rank(orbit, r)
                    assert rank(differential, r) == rank(differential.row_join(orbit), r) == expected
                    tested += 1

    # A three-level jet with noncommuting perturbations and a repeated eigenvalue.
    # The powered lift actually moves. Build its conjugator from the two projectors.
    n = 11*17
    a = Matrix([[1, 1, 0], [0, 1, 0], [0, 0, 2]])
    b, c = (Matrix(3, 3, [rng.randrange(n) for _ in range(9)]) for _ in range(2))
    lift = diag(a, a, a)
    lift[3:6, :3], lift[6:, 3:6], lift[6:, :3] = b, b, c
    power = matrix_power(lift, n, n)
    base = diag(*([matrix_power(a, n, n)] * 3))
    for r in (11, 17):
        second = pow(2, n, r)
        scale = pow(1-second, -1, r)
        lifted_projector = ((power - second*eye(9))*scale) % r
        base_projector = ((base - second*eye(9))*scale) % r
        transform = (lifted_projector*base_projector + (eye(9)-lifted_projector)*(eye(9)-base_projector)) % r
        assert power % r != base % r and rank(transform, r) == 9
        assert transform[:3, :3] == eye(3) and (power*transform-transform*base) % r == zeros(9)

    # x^2+1 is irreducible over both fields. Its divided difference is X+Y,
    # which gives the conjugator directly over the original ring.
    a = Matrix([[0, -1], [1, 0]])
    base = diag(a, a, a)
    lift = base.copy()
    b = Matrix([[2, 3], [5, 1]])
    lift[2:4, :2], lift[4:, 2:4] = b, b
    power = matrix_power(lift, 77, 77)
    for r in (7, 11):
        transform = ((power + base)*(-base)*pow(2, -1, r)) % r
        assert power % r != base % r and rank(transform, r) == 6
        assert transform[:2, :2] == eye(2) and (power*transform-transform*base) % r == zeros(6)

    # At ds>p take a single size-d nilpotent Jordan block, joined across s
    # coefficient layers by its wraparound entry epsilon. This is a ds-shift.
    for p, q, d, s in ((5, 11, 3, 2), (7, 13, 4, 2), (11, 17, 3, 4), (11, 17, 4, 4)):
        size = d*s
        shift = zeros(size)
        for i in range(size-1):
            shift[i+1, i] = 1
        difference = matrix_power(eye(size)+shift, p*q, p*q) - eye(size)
        assert (rank(difference, p), rank(difference, q)) == (size-p, 0)
        assert gcd(p*q, *map(int, difference)) == q
    boundary_jordan = Matrix([[1, 0, 0], [1, 1, 0], [0, 1, 1]])
    assert matrix_power(boundary_jordan, 55, 5) == eye(3)
    assert rank(power_jacobian(55, boundary_jordan), 5) > 0

    primes = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
    pairs = [(p, q) for i, p in enumerate(primes) for q in primes[i+1:]]
    assert all(shift_factor(p*q) == q for p, q in pairs)
    return tested, len(pairs)


if __name__ == '__main__':
    jacobians, shifts = check()
    print(f"N=91 gave factor 13. Checked {jacobians} local Jacobian images, four actual jet conjugators, "
          f"four sharp thresholds and {shifts} public shift factorizations.")
