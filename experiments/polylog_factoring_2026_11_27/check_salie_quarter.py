"""Check the Salié normalization and single-term sampling variance."""

from cmath import exp, pi
from math import cos, gcd, sqrt

from sympy import jacobi_symbol


def phase(n, value):
    return exp(2j*pi*(value % n)/n)


def check(p, q):
    n = p*q
    hidden = 1+p*(-2*pow(p,-1,q) % q)
    gauss = sqrt(n)*(1 if n % 4 == 1 else 1j)
    inv4 = pow(4, -1, n)
    terms = [(m, pow(m,-1,n), int(jacobi_symbol(m,n)))
             for m in range(1,n) if gcd(m,n) == 1]
    min_variance = float("inf")
    scales = 0
    for c in range(1,min(n,36)):
        if gcd(c,n) != 1:
            continue
        b = c*c*inv4 % n
        summands = [symbol*phase(n,m+b*inverse)/gauss for m,inverse,symbol in terms]
        total = sum(summands)
        expected = 2*cos(2*pi*c/n)+2*cos(2*pi*(c*hidden%n)/n)
        assert abs(total-expected) < 1e-8
        best_real_variance = sum(abs(term.real) for term in summands)**2 - total.real**2
        min_variance = min(min_variance,best_real_variance/n)
        scales += 1
    print(n, scales, f"{min_variance:.4f}")


if __name__ == "__main__":
    for p,q in [(11,13),(31,37),(101,109),(9,25),(27,49)]:
        check(p,q)
