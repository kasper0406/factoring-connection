"""List circle orbits from one supplied point, using quadratic auxiliary curves."""

from math import isqrt

from sympy import Matrix, Poly, symbols

from experiments.subthird.gemini.run_checks import lll_reduce


def circle_orbits(real, imag):
    n = real**2 + imag**2
    if not n:
        return {(0, 0)}
    blocks = 4 << ((3*n.bit_length() + 9) // 10)
    v = symbols("v")
    result = set()
    for offset in range(-blocks, blocks):
        h = Poly(v + offset, v)
        u = blocks**2 + h**2
        cosine = blocks**2 - h**2
        sine = 2*blocks*h
        x = real*cosine - imag*sine
        y = imag*cosine + real*sine
        polys = [u*u, x*u, y*u, x*x-y*y, 2*x*y]
        basis = Matrix([[p.nth(i) for i in range(5)] for p in polys])
        try:
            short = list(basis.lll()[0, :])
        except AssertionError:
            short = lll_reduce([[int(c) for c in row] for row in basis.tolist()])[0]
        if sum(map(abs, short)) >= blocks**4:
            raise ArithmeticError("auxiliary polynomial exceeds the arc bound")
        equation = Poly.from_list(short[::-1], v)
        for root in equation.ground_roots():
            if 0 <= root <= 1:
                denominator = u.eval(root)
                xx, yy = x.eval(root)/denominator, y.eval(root)/denominator
                if xx.is_Integer and yy.is_Integer:
                    result.add(tuple(sorted((abs(int(xx)), abs(int(yy))))))
    return result

def orbits_from_norm(n):
    """List circle orbits without a starting point or a factorization."""
    if n < 1:
        raise ValueError("expected a positive norm")
    blocks = 8 << ((3*n.bit_length() + 9) // 10)
    scale = 16 * (2*blocks + 1)**5
    radius_numerator = isqrt(n * scale**2)
    v = symbols("v")
    result = set()
    for offset in range((blocks + 1) // 2):
        h = Poly(v + offset, v)
        u = blocks**2 + h**2
        cosine = blocks**2 - h**2
        sine = 2*blocks*h
        polys = [scale*u*u, radius_numerator*cosine*u, radius_numerator*sine*u,
                 scale*n*(cosine*cosine-sine*sine), scale*n*2*cosine*sine]
        basis = Matrix([[p.nth(i) for i in range(5)] for p in polys])
        try:
            reduced, transform = basis.lll_transform()
            short, weights = reduced[0, :], transform[0, :]
        except AssertionError:
            short = Matrix([lll_reduce(basis.tolist())[0]])
            weights = short * basis.inv()
        a, b, c, d, e = weights
        # The surrogate error is at most (|b|+|c|)/scale.
        if sum(map(abs, short)) + (abs(b)+abs(c))*blocks**4 >= scale*blocks**4:
            raise ArithmeticError("radius-transfer certificate failed")
        even = a + b*v + d*(2*v*v-n)
        odd = c + 2*e*v
        equation = Poly(even**2 - (n-v*v)*odd**2, v)
        for x in equation.ground_roots():
            if x.is_Integer and 0 <= x*x <= n:
                y = isqrt(n - int(x)**2)
                if x*x + y*y == n:
                    result.add(tuple(sorted((abs(int(x)), y))))
    return result
