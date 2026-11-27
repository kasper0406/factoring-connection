"""Check adjacent-binomial carry witnesses by direct gcd and factor-assisted labels."""

from math import comb, gcd, isqrt
from itertools import combinations
from sympy import primerange

from check_lehman_window_baseline import cohorts


def indices(n):
    product = 1
    for d in primerange(5, 100*n.bit_length()+100):
        for j in range(1, d//2+1):
            yield d, j, (j-1)*n//d, j*n//d
        product *= d
        if product > (isqrt(n)+2)**5:
            return
    raise AssertionError('Prime range did not reach stopping product')


def first_hit(p,q,direct):
    n = p*q
    for d,j,k,m in indices(n):
        if d >= p:
            return None  # public trial division through d covers this case
        drop_p = (j*q)%d < ((j-1)*q)%d
        drop_q = (j*p)%d < ((j-1)*p)%d
        if direct:
            g = gcd(comb(m,k),n)
            assert g == (p if drop_p else 1)*(q if drop_q else 1)
        if drop_p != drop_q:
            return d,j


if __name__ == '__main__':
    small = 0
    for p,q in combinations(primerange(13,100),2):
        if q<2*p:
            first_hit(p,q,True)
            small += 1
    print('direct gcd inputs',small)
    for bits,planted,cases in cohorts():
        if not planted and bits in (32,48):
            hits = [first_hit(p,q,False) for p,q in cases]
            assert all(hits)
            print(bits,len(hits),max(d for d,j in hits))
