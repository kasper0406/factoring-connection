"""Check lower-digit isolation and the upper-digit frequency of the count."""

from fractions import Fraction
from math import gcd, prod
from random import Random

import numpy as np


def sign(a, mask, flip):
    return 1-2*((((a & mask).bit_count() & 1)^flip))


def weight(n, a, offset, width, mask, flip):
    return sign(a,mask,flip) if (a+offset)%n < width else 0


def input_blocks(n, c):
    square = n*n
    for x in range(1,n):
        if gcd(x,n)!=1:
            continue
        inverse = pow(x,-1,square)
        upper, lower = divmod((x+c*c*inverse)%square,n)
        yield x, lower, upper, 2*c*inverse%n, gcd(x*x-c*c,n)


def blocks(n, c):
    for _, lower, upper, velocity, d in input_blocks(n,c):
        yield lower,upper,velocity,d


def corrected(n,c,j,items,offset,width,mask,flip,select_input=False):
    total = 0
    half = (n-1)//2
    for lower,upper,velocity,d in items:
        w = weight(n,lower,offset,width,mask,flip)
        if w:
            residue = (upper+velocity*j)%d
            total += w*d*((half-residue)//d+1)
    for root in (-1,1):
        upper,lower = divmod((2*root*(c+n*j))%(n*n),n)
        if upper <= half:
            key = (root*c)%n if select_input else lower
            total -= n*weight(n,key,offset,width,mask,flip)
    return total


def direct(n,c,j,offset,width,mask,flip,select_input=False):
    square = n*n
    total = 0
    for x in range(1,square):
        if gcd(x,n)==1:
            upper,lower = divmod((x+(c+n*j)**2*pow(x,-1,square))%square,n)
            if upper <= (n-1)//2:
                key = x%n if select_input else lower
                total += weight(n,key,offset,width,mask,flip)
    for r in (-1,1):
        upper,lower = divmod((2*r*(c+n*j))%square,n)
        if upper <= (n-1)//2:
            key = (r*c)%n if select_input else lower
            total -= n*weight(n,key,offset,width,mask,flip)
    return total


def check(primes,rng):
    n = prod(primes)
    roots = [u for u in range(2,n-1) if u*u%n==1]
    width = max(1,n//(4*(len(roots)+2)))
    trials = 0
    while True:
        trials += 1
        c,offset = rng.randrange(1,n),rng.randrange(n)
        if gcd(c,n)!=1:
            continue
        selected = [u for u in roots if (2*c*u+offset)%n < width]
        if len(selected)==1:
            u = selected[0]
            break
    mask,flip = rng.randrange(1<<n.bit_length()),rng.randrange(2)
    items = list(blocks(n,c))
    totals = [0]*n
    for lower,_,_,d in items:
        totals[lower] += d
    assert Fraction(sum(v*v for v in totals),n*n) == prod(
        2+Fraction(2*(p-3),p*p) for p in primes
    )
    selected_items = [item for item in items if (item[0]+offset)%n < width]
    counts = [corrected(n,c,j,selected_items,offset,width,mask,flip) for j in range(n)]
    if n<=105:
        assert all(counts[j]==direct(n,c,j,offset,width,mask,flip)
                   for j in (0,1,n//2,n-1))
    values = (np.array(counts)-counts[0])/n
    lift = (u-n*((u*u-1)//n)*pow(2*u,-1,n))%(n*n)
    private = []
    for j in range(n):
        upper,lower = divmod((2*(c+n*j)*lift)%(n*n),n)
        private.append(weight(n,lower,offset,width,mask,flip)*int(upper<=(n-1)//2))
    rest = np.fft.fft(values-private)/n
    assert max(abs(rest[k]) for k in range(n) if gcd(k,n)==1)<1e-10
    peak = abs(np.fft.fft(values)[2*u%n]/n)
    divisor = gcd(2*u-2,n)
    assert peak>.3 and 1<divisor<n
    print(n,width,trials, f"{peak:.6f}",divisor)


if __name__ == "__main__":
    rng=Random(261133)
    for primes in [(3,5),(3,5,7),(3,5,7,11),(3,5,7,11,13)]:
        check(primes,rng)
