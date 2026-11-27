"""List circle orbits from one supplied point, using quadratic auxiliary curves."""

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
