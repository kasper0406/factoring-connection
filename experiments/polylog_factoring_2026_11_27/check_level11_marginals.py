"""Check the two CRT marginals of the discriminant-11 principal ideal bit.

Rows are right-ideal classes, with P the four-unit principal class and A
the six-unit class. The Brandt matrix for a prime r away from 11 has row
sum r+1, stationary weights (3/5, 2/5) from the inverse unit orders, and
second eigenvalue a_r from the weight-two form eta(z)^2 eta(11z)^2.
Thus B_PP=(3(r+1)+2a_r)/5 and B_AP=(3(r+1)-3a_r)/5.

Fixing a p-line gives I_p, with p-primary quotient O/I_p. At q its
localization is O_q. Its q-subideals are precisely the intersections
I_p ∩ I_q as the q-line varies. Left multiplication by a class
equivalence bijects those subideals with the corresponding Brandt row.
Consequently E[C_pq | L_p]=3/5 + a_q/(q+1) (C_p-3/5).
Putting lambda=a_q/(q+1), the centered projection has variance
lambda² Var(C_p). Its covariance with any f(L_p) has magnitude at
most |lambda| sd(f)/2 by Cauchy-Schwarz and Var(C_p) <= 1/4.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt
from random import Random

from sympy import Matrix

if __package__:
    from .level11_sampler import _inverse_splitting, principal_line_predicate
else:
    from level11_sampler import _inverse_splitting, principal_line_predicate


def eta_coefficients(limit):
    """Coefficients of t prod_{k>=1}(1-t^k)^2(1-t^(11k))^2."""
    coefficients = [0, 1] + [0] * (limit - 1)
    for k in range(1, limit):
        for stride in (k, k, 11*k, 11*k):
            for i in range(limit, stride - 1, -1):
                coefficients[i] -= coefficients[i - stride]
    return coefficients


def norm_points(n):
    """All generators, by 4(x²+xy+3y²)=(2x+y)²+11y²."""
    shells = [[] for _ in range(n + 1)]
    for y in range(-isqrt(4*n//11), isqrt(4*n//11) + 1):
        limit = isqrt(4*n - 11*y*y)
        for doubled in range(-limit + ((limit+y) % 2), limit + 1, 2):
            shells[(doubled*doubled + 11*y*y)//4].append(((doubled-y)//2, y))
    for i, shell in enumerate(shells):
        for x, y in shell:
            for z, w in shells[n-i]:
                yield (x, y, z, w)


def image_line(point, rows, r):
    """Return slope v/u of the image line, with r denoting the vertical line."""
    a, b, c, d = (sum(t*x for t, x in zip(row, point)) % r for row in rows)
    assert (a*d - b*c) % r == 0 and any((a, b, c, d))
    u, v = (a, c) if a or c else (b, d)
    return v*pow(u, -1, r) % r if u else r


def population_variance(values):
    mean = sum(values, Fraction()) / len(values)
    return sum((v - mean)**2 for v in values) / len(values)


def check_pair(p, q, coefficients):
    n = p*q
    # Use one composite splitting and reduce its matrix at *both* primes.
    # Independent prime splittings would generally permute the line labels.
    inverse = _inverse_splitting(n, Random(n))
    forward = inverse.adjugate()*pow(int(inverse.det()), -1, n)
    rows = [[int(x) % n for x in row] for row in forward.tolist()]
    p_lines = Counter(image_line(v, rows, p) for v in norm_points(p))
    q_lines = Counter(image_line(v, rows, q) for v in norm_points(q))
    fibers = Counter((image_line(v, rows, p), image_line(v, rows, q))
                     for v in norm_points(n))
    assert all(count == 4 for group in (p_lines, q_lines, fibers) for count in group.values())

    density = Fraction(3, 5)
    projected = []
    for r, s, lines, axis in [(p, q, p_lines, 0), (q, p, q_lines, 1)]:
        slope = Fraction(coefficients[s], s + 1)
        marginal = [Fraction(sum(cell[axis] == u for cell in fibers), s + 1)
                    for u in range(r + 1)]
        bits = [int(u in lines) for u in range(r + 1)]
        assert sum(bits) == density*(r + 1) + (1-density)*coefficients[r]
        assert marginal == [density + slope*(bit - density) for bit in bits]
        projected.append(population_variance(marginal))
        assert projected[-1] == slope**2 * population_variance(bits)
    mean = Fraction(len(fibers), (p+1)*(q+1))
    assert mean == density + (1-density)*Fraction(coefficients[p]*coefficients[q], (p+1)*(q+1))
    variance = mean*(1-mean)
    assert sum(projected) <= variance  # Centered one-line functions are orthogonal.
    if (p, q) == (19, 29):
        matrix = Matrix(p+1, q+1, lambda i, j: int((i, j) in fibers))
        assert matrix.to_DM().rank() == 20

    row_counts = Counter(i for i, _ in fibers)
    col_counts = Counter(j for _, j in fibers)
    # On a principal line the determinant gcd reaches p or q exactly when
    # one of the two local lines agrees with the target.
    for i in range(p+1):
        for j in range(q+1):
            row_probability = density + Fraction(coefficients[q], q+1)*(int(i in p_lines)-density)
            col_probability = density + Fraction(coefficients[p], p+1)*(int(j in q_lines)-density)
            double_hit = int((i, j) in fibers)
            predicted = (row_probability/(p+1) + col_probability/(q+1)
                         - Fraction(2*double_hit, (p+1)*(q+1))) / mean
            observed = Fraction(row_counts[i]+col_counts[j]-2*double_hit, len(fibers))
            assert observed == predicted

    if (p, q) == (3, 5):
        def crt(u, v):
            return (u + p*((v-u)*pow(p, -1, q) % q)) % n

        principal = principal_line_predicate(n, Random(n))
        principal_pairs = []
        for i in range(p + 1):
            for j in range(q + 1):
                pu, pv = (1, i) if i < p else (0, 1)
                qu, qv = (1, j) if j < q else (0, 1)
                pair = crt(pu, qu), crt(pv, qv)
                assert principal(*pair) == ((i, j) in fibers)
                if (i, j) in fibers:
                    principal_pairs.append((i, j, *pair))
        for i, j, u, v in principal_pairs:
            for k, l, a, b in principal_pairs:
                assert gcd(n, u*b-v*a) == (p if i == k else 1)*(q if j == l else 1)
    return len(fibers), (p+1)*(q+1), projected, variance - sum(projected)


if __name__ == "__main__":
    primes = (3, 5, 7, 13, 17, 19, 29)
    coefficients = eta_coefficients(max(primes))
    pairs = cells = principals = 0
    for p, q in combinations(primes, 2):
        hits, total, projections, interaction = check_pair(p, q, coefficients)
        pairs += 1
        cells += total
        principals += hits
        if (p, q) in ((3, 19), (19, 29)):
            print(f"{p}*{q}: principal {hits}/{total}, projected variances "
                  f"{projections[0]}, {projections[1]}, mixed variance {interaction}")
    print(f"{pairs} semiprimes, {cells} CRT cells, {principals} principal cells")
