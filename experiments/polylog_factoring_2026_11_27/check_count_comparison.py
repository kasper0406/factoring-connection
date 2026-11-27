"""Check the one-bit gap comparator on isolated root fixtures."""

from math import gcd, prod
from random import Random

import numpy as np

from check_square_digit_isolation import input_blocks


def check(components,rng):
    n = prod(components)
    mixed = [u for u in range(2,n-1) if u*u%n==1]
    width = max(1,n//(8*(len(mixed)+2)))
    while True:
        c,offset = rng.randrange(1,n),rng.randrange(n)
        if gcd(c,n)!=1:
            continue
        selected = [u for u in mixed if (c*u+offset)%n < width]
        if len(selected)==1:
            u = selected[0]
            break
    a,b = rng.randrange(n),rng.randrange(n)
    digits = np.arange(n,dtype=np.int64)
    residual, stationary = np.zeros(n,dtype=np.int64), np.zeros(n,dtype=np.int64)
    half = (n-1)//2
    for x,_,upper,slope,d in input_blocks(n,c):
        if d==1 or x in (c,n-c) or (x+offset)%n>=width:
            continue
        sign = 1 if (a*x+b)%n<=half else -1
        count = (half-(upper+slope*digits)%d)//d
        change = sign*d*(count-count[0])
        if d==n:
            stationary += change
        else:
            residual += change
    assert set(abs(stationary))=={0,n}
    true_bit = abs(stationary)>0
    magnitude = abs(stationary+residual)
    approximate = np.array([rng.getrandbits(1) for _ in range(n)])
    approximate[4*magnitude<=n]=0
    approximate[4*magnitude>=3*n]=1
    spectrum = np.fft.fft(approximate)/n
    coefficient = abs(spectrum[2*u%n])
    bad = np.mean(approximate!=true_bit)
    assert coefficient>.2
    assert 1<gcd(2*u-2,n)<n
    print(n,width,round(float(np.linalg.norm(residual/n)/np.sqrt(n)),5),
          round(float(bad),5),round(float(coefficient),6))


if __name__ == "__main__":
    rng = Random(261135)
    for components in [(3,5),(9,25),(9,5,49),(101,109),(9,25,7,11)]:
        check(components,rng)
