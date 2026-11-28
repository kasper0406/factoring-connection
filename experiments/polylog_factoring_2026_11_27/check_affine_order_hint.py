"""Verify the eN divisor and CRT intersection identities and the input-only decoder."""

from itertools import combinations, product
from math import gcd, isqrt

from sympy import isprime, primerange

if __package__:
    from .affine_order_hint import factor_from_affine_order_hint
else:
    from affine_order_hint import factor_from_affine_order_hint


def check_identities():
    count = 0
    for p, q in combinations(primerange(5, 50), 2):
        if q > 2*p:
            continue
        n = p*q
        for r, partner in ((p, q), (q, p)):
            for a, b, sign, e in product((1, 2, -3), (1, -2), (-1, 1), range(1, 41)):
                g = gcd(e, a*r+sign*b)
                offset = sign*b*partner
                assert abs(offset) <= abs(b)*(isqrt(2*n-1)+1)
                multiple = gcd(e*n, a*n+offset)
                assert multiple % (partner*g) == 0
                assert gcd(multiple, n) == partner
                count += 1
    return count


def check_hint(r, partner, e, sign, *, a=1, b=1, above=False):
    n = r*partner
    assert isprime(r) and isprime(partner) and max(r, partner) <= 2*min(r, partner)
    assert gcd(a*b*e, n) == 1
    if above:
        assert e**4 > n**3
    stats = {}
    split = factor_from_affine_order_hint(n, e, a=a, b=b, stats=stats)
    g = gcd(e, a*r+sign*b)
    assert stats['method'] == 'lattice' and stats['lower'] <= partner*g
    assert split in (r, partner) and 0 < stats['offsets_tested'] < stats['dimension']
    return n.bit_length(), e.bit_length(), abs(a*r+sign*b)//g, split, stats['dimension'], stats['offsets_tested']


def run():
    print('Signed divisor and intersection identities:', check_identities())
    p, q = 1099511627873, 1649267441959
    n = p*q
    high = [check_hint(p,q,(p+sign)*2**21,sign,above=True) for sign in (-1,1)]
    high.append(check_hint(q,p,(q-1)*2**21,-1,above=True))
    # Public replay of a planted hint with a genuinely missing cofactor.
    r, partner = 4722366482869645213801, 7083549724304467820953
    e = (r-1)//2 * (2**38+1)
    high.append(check_hint(r,partner,e,-1,above=True))
    assert high[-1][2] == 2
    print('Hints above N^(3/4): (N bits, e bits, missing cofactor, factor, dimension, offsets)')
    for record in high:
        print(record)

    count = 0
    for a, b, sign in product((-3, 3), (-2, 2), (-1, 1)):
        check_hint(p,q,abs(a*p+sign*b),sign,a=a,b=b)
        count += 1
    print('Certified signed rational-shift cases:', count)

    stats = {}
    assert factor_from_affine_order_hint(n, n-1, stats=stats) is None
    assert stats['method'] == 'no certificate' and stats['dimension'] == 0
    assert factor_from_affine_order_hint(n, 11, stats=stats) is None
    assert stats['method'] == 'small exponent' and stats['dimension'] == 0
    assert factor_from_affine_order_hint(n, 1) is None
    assert factor_from_affine_order_hint(n, n) is None
    assert factor_from_affine_order_hint(n, 2*p) == p
    assert factor_from_affine_order_hint(n, 11, b=q) == q
    assert factor_from_affine_order_hint(n, 11, a=p) == p
    # Both hidden orders miss this exponent's chosen certified threshold.
    assert factor_from_affine_order_hint(n, p+7, stats=stats) is None
    assert stats['method'] == 'lattice' and stats['offsets_tested'] >= 1
    print('Eight range, gcd and no-split controls passed.')


if __name__ == '__main__':
    run()
