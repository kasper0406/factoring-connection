"""Enumerate roots of x² = 1 mod n in a half-open interval by certified blocks."""

from math import isqrt

from sympy import Matrix, Poly, symbols

from experiments.subthird.gemini.run_checks import lll_reduce


def _block_roots(n, start, stop):
    center = (start + stop - 1) // 2
    radius = stop - 1 - center
    if radius == 0:
        if center * center % n == 1:
            yield center
        return

    x = symbols("x")
    quadratic = Poly(x**2 + (2*center % n)*x + (center**2 - 1) % n, x)
    for m in range(2, n.bit_length() + 1):
        rows = []
        power = Poly(1, x)
        scales = [radius**i for i in range(2*m)]
        for j in range(m):
            coeffs = [n**(m-j) * int(power.nth(i)) for i in range(2*j+1)]
            for shift in (0, 1):
                row = [0] * (2*m)
                for i, coefficient in enumerate(coeffs, shift):
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
                if start <= candidate < stop and candidate**2 % n == 1:
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
