"""Check marked norm pooling against primitive disks and actual fine Gram marks."""

from itertools import product
from math import gcd, isqrt

if __package__:
    from .check_mixed_ternary_embeddings import _all_frames
    from .marked_conic_slopes import _slope_line, _transport, marked_slope_partners
    from .pooled_marked_slopes import pooled_slope_factor, slope_form
else:
    from check_mixed_ternary_embeddings import _all_frames
    from marked_conic_slopes import _slope_line, _transport, marked_slope_partners
    from pooled_marked_slopes import pooled_slope_factor, slope_form


def check_disks():
    cases = needed = 0
    for n, A, D in [(77, 1, 1), (91, 17, 1), (143, 2, 2), (437, 17, 1)]:
        limit = 2*D*isqrt(2*n)
        assert limit < n
        norms = [(k, j, k*k+A*j*j) for k in range(-isqrt(limit), isqrt(limit)+1)
                 for j in range(isqrt(limit//A)+1) if 0 < k*k+A*j*j <= limit]
        for ell in (1, 2, 3, 4, 5, 8, 15):
            if gcd(n, ell) > 1:
                continue
            for left, right in product(range(ell), repeat=2):
                if gcd(ell, left, right) != 1:
                    continue
                a, b, c = slope_form(A, ell, (left, right))
                assert abs(b) <= a <= c and 4*a*c-b*b == 4*A*ell*ell
                line_hits = [(k, j, gcd(n, value)) for k, j, value in norms
                             if (k*right-j*left) % ell == 0 and gcd(n, value) > 1]
                factor = pooled_slope_factor(n, A, D, ell, (left, right))
                if any(gcd(k, j) == 1 for k, j, _ in line_hits):
                    assert factor is not None
                    needed += 1
                if factor is not None:
                    assert 1 < factor < n and n % factor == 0
                    assert any(d == factor for _, _, d in line_hits)
                cases += 1
    return cases, needed


def _mark(n, x, u, z, ell):
    den, transform = _transport(u, x, ell)
    first, _, last = transform(z)
    row = first+den*n
    line = _slope_line((2*den*n-(x*x+1)*row, last, row), den, {ell: 1})
    return den, line


def check_frames():
    cases, transported = 0, 0
    for n, x, ell in [(15, 1, 7), (77, 1, 13), (437, 4, 11)]:
        for u, v in _all_frames(n, x):
            if gcd(n, *v) in (1, n):
                continue
            z = tuple(t % ell for t in v)
            den, (left, right, modulus) = _mark(n, x, u, z, ell)
            factor = pooled_slope_factor(n, x*x+1, den, modulus, (left, right))
            assert 1 < factor < n and n % factor == 0
            cases += 1
            transported += den > 1
    # A proper gcd alone does not answer the mark.
    n, x, ell, u, z = 15, 1, 7, (1,1,0), (3,3,2)
    den, (left,right,modulus) = _mark(n, x, u, z, ell)
    assert pooled_slope_factor(n, 2, den, modulus, (left,right)) == 3
    assert not marked_slope_partners(n, x, u, z, {ell: 1})
    return cases, transported


def check_regimes():
    work = {}
    assert pooled_slope_factor(77, 10, 2, 1, (0, 0), stats=work) in (7, 11)
    assert work['saturated_blocks'] == 1
    assert pooled_slope_factor(77, 2, 77, 13, (1, 0), stats=work) == 7
    assert work['method'] == 'small input'
    assert pooled_slope_factor(77, 2, 7, 13, (1, 0), stats=work) == 7
    assert work['method'] == 'transport gcd'
    assert pooled_slope_factor(77, 2, 1, 97, (1, 0), stats=work) is None
    assert work['method'] == 'single direction' and work['tiles'] == 0
    assert pooled_slope_factor(77, 2, 1, 997, (31, 17), stats=work) is None
    assert work['method'] == 'empty ellipse'


def check_planted():
    p, q, k, j = 4343699, 6515557, 1959, 503
    n = p*q
    assert k*k+2*j*j == p and q % 8 == 5 and gcd(k,j) == 1
    v = (-n+2*q*j*j, -2*q*j*j, 2*q*k*j)
    assert sum(t*t for t in v) == n*n and v[0]+v[1] == -n and gcd(n, *v) == q
    for ell in (101, 1009, 10007, 10000019):
        den, (left,right,modulus) = _mark(n, 1, (1,1,0), tuple(t % ell for t in v), ell)
        assert den == 1 and modulus == ell and (k*right-j*left) % ell == 0
        work = {}
        assert pooled_slope_factor(n, 2, den, modulus, (left,right), stats=work) == p
        print(f'45-bit marked input, modulus {ell}: {work["method"]}, '
              f'{work["quadratic_factors"]} quadratic factors, {work["block_values"]} '
              f'block values, {work["scalar_gcds"]} scalar refinements')


if __name__ == '__main__':
    cases, needed = check_disks()
    accepted, transported = check_frames()
    print(f'{cases} line generators checked; all {needed} primitive norm hits covered')
    print(f'{accepted} informative frame marks found factors, {transported} at transport denominator >1')
    check_regimes()
    check_planted()
