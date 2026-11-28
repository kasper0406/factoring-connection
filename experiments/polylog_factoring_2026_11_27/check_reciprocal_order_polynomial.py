"""Test the partner-prime certificate for a fixed polynomial order hint.

Let N=r*s, g=gcd(e,P(r)), and d=deg P, with P(0)=±1. The monic
F(X)=P(0)*X^d P(N/X) has F(s)=P(0)*s^d P(r). Thus the unknown divisor
B=s^d*g of M=e*N^d divides F(s). For fixed P, balanced primes, and
|P(r)|/g=N^o(1), B=N^{d-o(1)}. If e=N^alpha, the usual degree-d
univariate Coppersmith bound M^((log_M B)^2/d) is N^{d/(d+alpha)-o(1)}.
It exceeds sqrt(2*N) by a fixed power whenever alpha<=d-epsilon.

This finite check uses specified lattice dimensions and reports the divisor
lower bound actually certified by its reduced row. It invokes the same LLL
and integer-root isolation as affine_prime_window. No factor enters the reader.
"""

from math import gcd, isqrt

from sympy import Poly, Symbol, integer_nthroot

if __package__:
    from .affine_prime_window import _integer_roots, lll_reduce
else:
    from affine_prime_window import _integer_roots, lll_reduce


def reciprocal_roots(n, e, coefficients, *, dimension, multiplicity):
    """All balanced-window s with gcd(e*n^d,F(s)) >= the returned lower bound.

    The coefficients of P are supplied in increasing degree. The lower bound
    is independent of the roots and incorporates the actual row norm.
    """
    degree = len(coefficients) - 1
    if n < 6 or e < 1 or degree < 1 or not coefficients[-1] or abs(coefficients[0]) != 1:
        raise ValueError("expected n>=6, a positive exponent, and positive degree with P(0)=±1")
    if multiplicity < 1 or dimension < degree*multiplicity:
        raise ValueError("require positive multiplicity and dimension at least degree*multiplicity")
    modulus = e*n**degree
    lo, hi = isqrt((n-1)//2) + 1, isqrt(2*n)
    center = (lo+hi)//2
    radius = max(1, hi-center)
    x = Symbol("x")
    reciprocal = Poly.from_list(
        [coefficients[0]*c*n**i for i, c in enumerate(coefficients)], x, domain="ZZ"
    )
    shifted = reciprocal.shift(center).trunc(modulus)
    powers = [shifted**i for i in range(multiplicity+1)]
    scales = [radius**i for i in range(dimension)]
    rows = []
    for i in range(dimension):
        power = min(i//degree, multiplicity)
        shift = i-degree*power
        # Each row polynomial X^shift M^(m-power) F(center+X)^power
        # evaluates to a multiple of b^m whenever b divides M and F(s).
        factor = modulus**(multiplicity-power)
        row = [0]*dimension
        for j, coefficient in enumerate(reversed(powers[power].all_coeffs()), shift):
            row[j] = int(coefficient)*factor*scales[j]
        rows.append(row)
    reduced = lll_reduce(rows)
    row = min((list(map(int, v)) for v in reduced), key=lambda v: sum(map(abs, v)))
    assert any(row) and all(a % scale == 0 for a, scale in zip(row, scales))
    lower = int(integer_nthroot(sum(map(abs, row)), multiplicity)[0]) + 1
    # For b>=lower, the evaluated auxiliary polynomial has |value|<b^m.
    aux = Poly.from_list([row[i]//scales[i] for i in reversed(range(dimension))], x, domain="ZZ")
    roots = []
    for offset in _integer_roots(aux, radius):
        candidate = center+offset
        if lo <= candidate <= hi and gcd(modulus, int(reciprocal.eval(candidate))) >= lower:
            roots.append(candidate)
    return tuple(sorted(roots)), lower


def main():
    # Independently check completeness of the returned list, including composite
    # roots and both signs of the constant coefficient.
    checks = simultaneous = 0
    for n, e in ((15, 403), (35, 1767), (48, 3139), (55, 31), (77, 399), (143, 427)):
        for coefficients in ((-1, 1), (1, 1, 1), (-1, 0, 1)):
            roots, lower = reciprocal_roots(n, e, coefficients, dimension=7, multiplicity=2)
            degree = len(coefficients)-1
            expected = []
            for s in range(isqrt((n-1)//2)+1, isqrt(2*n)+1):
                value = sum(c*n**i*s**(degree-i) for i, c in enumerate(coefficients))
                if gcd(e*n**degree, value) >= lower:
                    expected.append(s)
            assert list(roots) == expected
            checks += 1
            simultaneous += len(roots) >= 2
    assert simultaneous

    r, s = 17179869209, 22333829989
    n = r*s
    order = r*r+r+1
    e = (order//13)*241
    assert order//gcd(e, order) == 13 and n < e < n*n
    roots, lower = reciprocal_roots(n, e, (1, 1, 1), dimension=10, multiplicity=3)
    witness = s*s*gcd(e, order)
    divisors = tuple(x for x in roots if 1 < x < n and n % x == 0)
    assert witness >= lower and s in divisors
    print(f"{checks} exhaustive window checks passed, {simultaneous} with simultaneous roots")
    print(f"n={n} e={e} certified_lower={lower} witness={witness} divisors={divisors}")


if __name__ == "__main__":
    main()
