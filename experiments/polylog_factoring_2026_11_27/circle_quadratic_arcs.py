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

def orbits_from_norm(n, degree=2):
    """List all circle orbits using certified harmonic arcs of the given degree."""
    if n < 1 or degree < 1:
        raise ValueError("expected a positive norm and a positive degree")
    variable = symbols("v")
    v = Poly(variable, variable)
    dimension = 2*degree + 1
    blocks = 2
    # Store only the arcs still awaiting a certificate at this resolution.
    offsets = [0]
    result = set()
    while offsets:
        scale = 32*dimension*(2*blocks + 1)**(2*degree)
        denominator = blocks**(2*degree)
        radial = [isqrt(scale**2 * n**j) for j in range(1, degree + 1)]
        children = []
        for offset in offsets:
            h = v + offset
            u = blocks**2 + h*h
            cosine, sine = blocks**2 - h*h, 2*blocks*h
            re, im = v.one, v.zero
            polys = [scale*u**degree]
            for j, multiplier in enumerate(radial, 1):
                re, im = cosine*re - sine*im, sine*re + cosine*im
                factor = multiplier*u**(degree-j)
                polys.extend((factor*re, factor*im))
            basis = Matrix([[p.nth(i) for i in range(dimension)] for p in polys])
            try:
                reduced, transform = basis.lll_transform()
                short, weights = reduced[0, :], transform[0, :]
            except AssertionError:
                rows = [[int(c) for c in row] for row in basis.tolist()]
                short = Matrix([lll_reduce(rows)[0]])
                weights = short * basis.inv()
            error_weight = sum(abs(weights[2*j-1]) + abs(weights[2*j])
                               for j, multiplier in enumerate(radial, 1)
                               if multiplier**2 != scale**2 * n**j)
            if sum(map(abs, short)) + error_weight*denominator >= scale*denominator:
                children.extend((2*offset, 2*offset+1))
                continue
            # Represent (x+iy)^j as re(x) + iy*im(x), with y^2=n-x^2.
            re, im = v.one, v.zero
            y_squared = n-v*v
            even, odd = weights[0]*v.one, v.zero
            for j in range(1, degree+1):
                re, im = v*re - y_squared*im, re + v*im
                even += weights[2*j-1]*re
                odd += weights[2*j]*im
            equation = even**2 - y_squared*odd**2
            for x in equation.ground_roots():
                if x.is_Integer and x*x <= n:
                    y = isqrt(n-int(x)**2)
                    if x*x + y*y == n:
                        result.add(tuple(sorted((abs(int(x)), y))))
        offsets = children
        blocks *= 2
    return result
