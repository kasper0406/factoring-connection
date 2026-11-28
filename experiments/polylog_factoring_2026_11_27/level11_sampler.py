"""Sample representations by x²+xy+3y²+z²+zw+3w² without factoring the norm."""

from fractions import Fraction
from itertools import product
from math import gcd, isqrt
from random import SystemRandom

from sympy import Matrix, eye
from sympy.matrices.normalforms import hermite_normal_form

if __package__:
    from .four_square_incidence import _other_prime_powers
    from .four_square_sampler import gaussian_norm_pair
else:
    from four_square_incidence import _other_prime_powers
    from four_square_sampler import gaussian_norm_pair


def norm(point):
    x, y, z, w = point
    return x*x + x*y + 3*y*y + z*z + z*w + 3*w*w


def _primitive_pair(n, rng):
    while True:
        u, v = rng.randrange(n), rng.randrange(n)
        if gcd(u, v, n) == 1:
            return u, v


def _coordinate_inverse(n, omega, j):
    columns = [m.reshape(4, 1) for m in (eye(2), omega, j, omega*j)]
    basis = Matrix.hstack(*columns).applyfunc(lambda x: x % n)
    # An adjugate avoids pivoting on a nonunit over a composite modulus.
    return (basis.adjugate()*pow(int(basis.det()), -1, n)).applyfunc(lambda x: x % n)


def _inverse_splitting(n, rng):
    """Matrix entries to order coordinates over Z/n for n coprime to 11."""
    if n < 1 or n % 11 == 0:
        raise ValueError("the splitting needs a positive modulus coprime to 11")
    two = n & -n
    odd = n//two
    if odd == 1:
        odd_inverse = eye(4)
    else:
        a, b = gaussian_norm_pair(odd, -11, rng)
        j = Matrix([[0, 1], [-1, 0]])
        omega = (eye(2) + Matrix([[a, b], [b, -a]]))*pow(2, -1, odd)
        odd_inverse = _coordinate_inverse(odd, omega, j)
    if two == 1:
        return odd_inverse

    # Over Z/2^e, multiplication by omega and conjugation act on (1, omega).
    # Norm(s+omega) = -1 when s²+s+4 = 0. The derivative 2s+1 is odd.
    s, precision = 0, 2
    while precision < two:
        if (s*s + s + 4) & precision:
            s += precision
        precision *= 2
    omega = Matrix([[0, -3], [1, 1]])
    conjugation = Matrix([[1, 1], [0, -1]])
    even_inverse = _coordinate_inverse(two, omega, (s*eye(2) + omega)*conjugation)
    if odd == 1:
        return even_inverse
    difference = (odd_inverse - even_inverse)*pow(two, -1, odd)
    return even_inverse + two*difference.applyfunc(lambda x: x % odd)


def _remove_eleven(n):
    if n < 1:
        raise ValueError("the norm must be positive")
    exponent = 0
    while n % 11 == 0:
        exponent += 1
        n //= 11
    return n, exponent


def _restore_eleven(point, exponent):
    # Multiplication by sqrt(-11)=2*omega-1; its square is the scalar -11.
    if exponent % 2:
        x, y, z, w = point
        point = (-x-6*y, 2*x+y, -z-6*w, 2*z+w)
    scale = (-11)**(exponent//2)
    return tuple(scale*x for x in point)


def _reduced_basis(basis):
    """Four-dimensional LLL with exact rounding even for large input integers."""
    rows = [[int(x) for x in row] for row in basis.tolist()]
    i = 1
    while i < len(rows):
        orthogonal = []
        lengths = []
        for row in rows[:i]:
            vector = list(map(Fraction, row))
            for earlier, length in zip(orthogonal, lengths):
                mu = sum(x*y for x, y in zip(row, earlier))/length
                vector = [x - mu*y for x, y in zip(vector, earlier)]
            orthogonal.append(vector)
            lengths.append(sum(x*x for x in vector))
        for j in reversed(range(i)):
            mu = sum(x*y for x, y in zip(rows[i], orthogonal[j]))/lengths[j]
            shift = round(mu)
            if shift:
                rows[i] = [x - shift*y for x, y in zip(rows[i], rows[j])]
        # Compare the projection away from rows 0,...,i−2 with the previous length.
        projected = sum(x*x for x in rows[i]) - sum(
            sum(x*y for x, y in zip(rows[i], vector))**2/length
            for vector, length in zip(orthogonal[:-1], lengths[:-1])
        )
        if 4*projected >= 3*lengths[-1]:
            i += 1
        else:
            rows[i], rows[i-1] = rows[i-1], rows[i]
            i = max(1, i-1)
    return Matrix(rows)


def _line_basis(n, inverse, line):
    """Column HNF for the right ideal of a primitive projective line."""
    u, v = line
    # Matrices whose columns lie on (u,v) span this ideal modulo n.
    lifts = (inverse * Matrix([[u, 0], [0, u], [v, 0], [0, v]])).applyfunc(lambda x: x % n)
    return hermite_normal_form((n*eye(4)).row_join(lifts))


def _shell_points(columns, n):
    """List every point of norm n in a full-rank column lattice."""
    basis = _reduced_basis(columns.T)
    dual = basis.inv()
    # norm(x) >= ||x||²/2. Apply Cauchy–Schwarz to each dual basis vector.
    bounds = [isqrt(int(2*n*sum(v*v for v in dual.col(j)))) for j in range(4)]
    rows = [list(map(int, row)) for row in basis.tolist()]
    points = []
    for coeffs in product(*(range(-r, r+1) for r in bounds)):
        point = tuple(sum(c*x for c, x in zip(coeffs, col)) for col in zip(*rows))
        if norm(point) == n:
            points.append(point)
    return points


def _line_points(n, inverse, line):
    """Enumerate the norm-n generators of the right ideal of a projective line."""
    points = _shell_points(_line_basis(n, inverse, line), n)
    if len(points) not in (0, 4):
        raise ArithmeticError("a principal line ideal must have four norm-n generators")
    return points


def _primitive_trial(n, inverse, rng):
    points = _line_points(n, inverse, _primitive_pair(n, rng))
    return points[rng.randrange(len(points))] if points else None


def principal_line_predicate(n, rng=None):
    """Set up a test for lines in P¹(Z/n), using one coordinate splitting.

    Requires 11 coprime to n. The returned function takes a primitive pair
    (u, v) and reports whether its line ideal has a generator of norm n.
    """
    if rng is None:
        rng = SystemRandom()
    inverse = _inverse_splitting(n, rng)

    def principal(u, v):
        if gcd(u, v, n) != 1:
            raise ValueError("the coordinates must generate the unit ideal modulo n")
        return bool(_line_points(n, inverse, (u % n, v % n)))

    return principal


def sample_primitive(n, rng=None):
    """Uniform signed ordered primitive point of norm n. The shell is empty if 121|n."""
    if rng is None:
        rng = SystemRandom()
    n, exponent = _remove_eleven(n)
    if exponent >= 2:
        raise ValueError("the primitive shell is empty when the norm is divisible by 121")
    inverse = _inverse_splitting(n, rng)
    while True:
        point = _primitive_trial(n, inverse, rng)
        if point is not None:
            return _restore_eleven(point, exponent)


def _propose_content(n, rng):
    """Each possible content g gets unnormalised proposal probability proportional to 1/g²."""
    limit = isqrt(n)
    highest = limit.bit_length() - 1
    ticket = rng.randrange(1, 1 << (highest + 1))
    start = 1 << (highest - (ticket.bit_length() - 1))
    g = rng.randrange(start, 2*start)
    if g <= limit and n % (g*g) == 0 and rng.randrange(g*g) < start*start:
        return g
    return None


def sample(n, rng=None):
    """Uniform signed ordered point of any positive norm n, including all contents."""
    if rng is None:
        rng = SystemRandom()
    n, exponent = _remove_eleven(n)
    inverse = _inverse_splitting(n, rng)
    while True:
        g = _propose_content(n, rng)
        if g is None:
            continue
        m = n//(g*g)
        other = _other_prime_powers(n, m)
        u, _ = _primitive_pair(other, rng)
        if gcd(u, other) != 1:
            continue
        # Only one ideal trial: its acceptance supplies the primitive theta weight.
        point = _primitive_trial(m, inverse, rng)
        if point is not None:
            return _restore_eleven(tuple(g*x for x in point), exponent)
