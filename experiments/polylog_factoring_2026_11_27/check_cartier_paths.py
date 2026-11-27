"""Compare composite-index curve coefficients with truncated permutation paths."""

from itertools import combinations
from math import comb, gcd

from sympy import primerange


def coefficient(n, d, i, j):
    # [x^(n*i-j)] (1+x^d)^((n-1)/2)
    k, rem = divmod(n*i-j, d)
    if rem or not 0 <= k <= (n-1)//2:
        return 0
    return comb((n-1)//2, k)


def check():
    hits = {5:0, 7:0, 11:0}
    total = any_hit = 0
    for p,q in combinations(primerange(13,100),2):
        if q >= 2*p:
            continue
        n = p*q
        total += 1
        split = False
        for d in hits:
            g = (d-1)//2
            found = False
            for i in range(1,g+1):
                for j in range(1,g+1):
                    value = coefficient(n,d,i,j)
                    for prime, partner in ((p,q),(q,p)):
                        intermediate = partner*i % d
                        path = (1 <= intermediate <= g and prime*intermediate % d == j)
                        assert (value % prime != 0) == path
                    found |= 1 < gcd(value,n) < n
            hits[d] += found
            split |= found
        any_hit += split
    print(total,hits,any_hit)


if __name__ == '__main__':
    check()
