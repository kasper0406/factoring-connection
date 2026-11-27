"""Check input-digit isolation, including inputs with repeated prime factors."""

from fractions import Fraction
from math import gcd, prod
from random import Random

import numpy as np

from check_square_digit_isolation import corrected, direct, input_blocks, weight


def check(components, rng):
    n = prod(p**e for p,e in components)
    mixed = [u for u in range(2,n-1) if u*u%n==1]
    width = max(1,n//(4*(len(mixed)+2)))
    while True:
        c,offset = rng.randrange(1,n),rng.randrange(n)
        if gcd(c,n)!=1:
            continue
        roots = [u for u in mixed if (c*u+offset)%n < width]
        if len(roots)==1:
            u=roots[0]
            break
    mask,flip = rng.randrange(1<<n.bit_length()),rng.randrange(2)
    all_blocks = list(input_blocks(n,c))
    assert Fraction(sum(d*d for _,_,_,_,d in all_blocks),n*n) == prod(
        2+Fraction(2,p)-Fraction(p+3,p**(e+1)) for p,e in components
    )
    selected = [(x,upper,velocity,d) for x,_,upper,velocity,d in all_blocks
                if (x+offset)%n < width]
    counts = [corrected(n,c,j,selected,offset,width,mask,flip,select_input=True)
              for j in range(n)]
    if n<=225:
        assert all(counts[j]==direct(n,c,j,offset,width,mask,flip,select_input=True)
                   for j in (0,1,n//2,n-1))
    values = (np.array(counts)-counts[0])/n
    lift = (u-n*((u*u-1)//n)*pow(2*u,-1,n))%(n*n)
    sign = weight(n,c*u%n,offset,width,mask,flip)
    private = [sign*int(((2*(c+n*j)*lift)%(n*n))//n <= (n-1)//2) for j in range(n)]
    remainder = np.fft.fft(values-private)/n
    assert max(abs(remainder[k]) for k in range(n) if gcd(k,n)==1)<1e-10
    peak = abs(np.fft.fft(values)[2*u%n]/n)
    divisor = gcd(2*u-2,n)
    assert peak>.3 and 1<divisor<n
    print(n,width,f"{peak:.6f}",divisor)


if __name__ == "__main__":
    rng=Random(261134)
    for components in [((3,1),(5,1)),((3,2),(5,2)),((3,3),(5,2),(7,1)),
                       ((3,2),(5,1),(7,2)),((3,2),(5,2),(7,1),(11,1))]:
        check(components,rng)
