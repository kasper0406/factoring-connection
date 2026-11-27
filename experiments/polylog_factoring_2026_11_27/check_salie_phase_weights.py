"""Check phase coefficients and their quadratic-residuosity predicate."""

from math import gcd

from sympy import jacobi_symbol


def table(n, c, units, character):
    b = c*c*pow(4,-1,n) % n
    coefficients = [0]*n
    for m in units:
        coefficients[(m+b*pow(m,-1,n)) % n] += character[m]
    return coefficients


def check(p, q):
    n = p*q
    u = 1+p*(-2*pow(p,-1,q) % q)
    units = [c for c in range(n) if gcd(c,n) == 1]
    chars = [int(jacobi_symbol(c,n)) for c in range(n)]
    tested = positive = inert = 0
    for c in units[:20]:
        coefficients = table(n,c,units,chars)
        for t,w in enumerate(coefficients):
            expected = sum(chars[(t+x)%n] for x in (c,-c,c*u,-c*u))
            assert w == expected
            tested += 1
            d = (t*t-c*c)%n
            if chars[d] != 1:
                continue
            square = pow(d,(p-1)//2,p) == pow(d,(q-1)//2,q) == 1
            assert w == (4*chars[(t+c)%n] if square else 0)
            positive += square
            inert += not square
    print(n, tested, positive, inert)
    return tested


if __name__ == "__main__":
    total = sum(check(p,q) for p,q in [(11,13),(31,37),(101,109)])
    print(f"{total} phase coefficients checked")
