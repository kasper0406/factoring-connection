"""Compute Res(H(X), (X-c)^m-u) modulo N without factoring N."""

from functools import cache
from math import gcd, prod

from sympy import Matrix, Poly, Symbol


def multiply(a, b, h, n):
    """Multiply residues in (Z/n)[X]/H. Coefficients run from low to high."""
    degree = len(h) - 1
    coefficients = [0] * (len(a) + len(b) - 1)
    for i, left in enumerate(a):
        for j, right in enumerate(b):
            coefficients[i + j] += left * right
    for k in range(len(coefficients) - 1, degree - 1, -1):
        lead = coefficients[k] % n
        for i in range(degree):
            coefficients[k - degree + i] -= lead * h[i]
    return tuple(v % n for v in (coefficients + [0] * degree)[:degree])


def power(value, exponent, h, n):
    """Raise a quotient-algebra element to a nonnegative integer power."""
    if exponent < 0:
        raise ValueError("expected a nonnegative exponent")
    base = multiply((1,), value, h, n)
    result = multiply((1,), (1,), h, n)
    while exponent > 0:
        if exponent & 1:
            result = multiply(result, base, h, n)
        exponent >>= 1
        if exponent:
            base = multiply(base, base, h, n)
    return result


def norm(value, h, n):
    degree = len(h) - 1
    columns = [value]
    for _ in range(1, degree):
        columns.append(multiply(columns[-1], (0, 1), h, n))
    # The transpose of the multiplication matrix has the same determinant.
    return int(Matrix(columns).det(method="berkowitz")) % n


def translated_resultant(n, h, c, m, u):
    """H is monic and nonconstant. Binary powering costs O(log m) products."""
    if n < 2 or m < 1 or len(h) < 2 or h[-1] != 1:
        raise ValueError("expected n >= 2, m >= 1 and a nonconstant monic H")
    value = power((-c, 1), m, h, n)
    return norm(((value[0] - u) % n, *value[1:]), h, n)


@cache
def field_hits(r, h, m, u):
    """Independently check geometric intersections by polynomial gcd over F_r."""
    x = Symbol("x")
    polynomial = Poly.from_list(h[::-1], x, modulus=r)
    shifts = set()
    for c in range(r):
        binomial = Poly(x - c, x, modulus=r)**m - u
        if polynomial.gcd(binomial).degree() > 0:
            shifts.add(c)
    assert len(shifts) <= (len(h) - 1) * gcd(m, r - 1)
    return shifts


def check():
    polynomials = ((-1, 1), (1, 0, 1), (-2, 1, 0, 1), (0, 1, 2, 2, 1))
    trials = splits = 0
    for p, q in ((5, 7), (7, 11), (11, 13)):
        n = p * q
        for h in polynomials:
            for m in (1, 2, 3, 6, 11):
                for u in (-1, 1, 2):
                    at_p, at_q = field_hits(p, h, m, u), field_hits(q, h, m, u)
                    for c in range(n):
                        actual = gcd(translated_resultant(n, h, c, m, u), n)
                        expected = (p if c % p in at_p else 1) * (q if c % q in at_q else 1)
                        assert actual == expected
                        trials += 1
                        splits += 1 < actual < n
    # (i-1)^4 = -4. The common root at 7 lies in its quadratic extension.
    assert translated_resultant(35, (1, 0, 1), 1, 4, 3) == 14
    assert all((x * x + 1) % 7 for x in range(7))

    # Over the integers this cubic has roots -3, 2, 4.
    n, m, c, u = 31439, 2**511 + 101, 5, 7
    actual = translated_resultant(n, (24, -10, -3, 1), c, m, u)
    assert actual == prod(pow(root - c, m, n) - u for root in (-3, 2, 4)) % n

    # A public seed yields the trace of a primitive fifteenth root at both primes.
    t = 3
    for _ in range(4):
        t = pow(t, n, n)
    assert t == 10994 and (t**4 - t**3 - 4 * t * t + 4 * t + 1) % n == 0
    trace_four = ((t * t - 2)**2 - 2) % n
    seed = (3 + trace_four) % n
    assert seed == 16683 and gcd(seed * seed + seed + 1, n) == 211
    # Ψ_15(X-3) reduced modulo X(X+1)(X²+X+1).
    assert norm((61, -108, 57, -15), (0, 1, 2, 2, 1), n) == 12449
    return trials, splits, field_hits.cache_info().currsize, actual


if __name__ == "__main__":
    trials, splits, field_bounds, stress_value = check()
    print(f"{trials} translations, {splits} proper gcds, {field_bounds} field bounds")
    print(f"512-bit exponent residue {stress_value}")
