"""Enumerate roots of x² = 1 mod n in a half-open interval by certified blocks."""

from math import isqrt

from sympy import Matrix, Poly, symbols

from experiments.subthird.gemini.run_checks import lll_reduce


def _block_roots(n, start, stop, shift=0, target=1):
    center = (start + stop - 1) // 2
    radius = stop - 1 - center
    if radius == 0:
        if ((center + shift)**2 - target) % n == 0:
            yield center
        return

    x = symbols("x")
    offset = center + shift
    quadratic = Poly(x**2 + (2*offset % n)*x + (offset**2 - target) % n, x)
    for m in range(2, n.bit_length() + 1):
        rows = []
        power = Poly(1, x)
        scales = [radius**i for i in range(2*m)]
        for j in range(m):
            coeffs = [n**(m-j) * int(power.nth(i)) for i in range(2*j+1)]
            for degree_shift in (0, 1):
                row = [0] * (2*m)
                for i, coefficient in enumerate(coeffs, degree_shift):
                    row[i] = coefficient * scales[i]
                rows.append(row)
            power *= quadratic
        try:
            short = Matrix(rows).lll()[0, :]
        except (AssertionError, OverflowError):
            short = lll_reduce(rows)[0]
        if sum(map(abs, short)) >= n**m:
            continue
        coefficients = [int(c)//s for c, s in zip(short, scales)]
        polynomial = Poly.from_list(coefficients[::-1], x)
        for root in polynomial.ground_roots():
            if root.is_Integer:
                candidate = center + int(root)
                if start <= candidate < stop and ((candidate + shift)**2 - target) % n == 0:
                    yield candidate
        return
    raise ArithmeticError("no certified row at the proven dimension")


def window_roots(n, start, stop):
    """Yield all roots in [start, stop), within the canonical residue range."""
    if n < 3 or n % 2 == 0 or not 0 <= start <= stop <= n:
        raise ValueError("expected odd n >= 3 and 0 <= start <= stop <= n")
    width = 2*(isqrt(n)//16) + 1
    for left in range(start, stop, width):
        yield from _block_roots(n, left, min(left + width, stop))


def progression_roots(n, first, step, count):
    """Yield canonical roots in the first count elements of a unit-step progression."""
    if n < 3 or n % 2 == 0 or not 0 <= count <= n:
        raise ValueError("expected odd n >= 3 and 0 <= count <= n")
    inverse = pow(step, -1, n)
    shift = first * inverse % n
    target = inverse**2 % n
    width = 2*(isqrt(n)//16) + 1
    for left in range(0, count, width):
        for index in _block_roots(n, left, min(left + width, count), shift, target):
            yield (first + step*index) % n
