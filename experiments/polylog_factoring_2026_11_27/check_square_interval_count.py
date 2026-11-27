"""Check the unweighted hyperbola interval and its hidden unit frequency."""

from cmath import exp, pi
from math import gcd, prod

from check_quarter_fourier import exact_spectrum


def inside(n, residue):
    square = n*n
    residue %= square
    return min(residue,square-residue) <= square//4


def direct_count(n, c, units):
    return sum(inside(n,x+c*c*inverse) for x,inverse in units)


def analyze(prime_powers):
    n = prod(prime_powers)
    square = n*n
    p,q = prime_powers
    u = 1+p*(-2*pow(p,-1,q)%q)
    lift = (u-n*((u*u-1)//n)*pow(2*u,-1,n))%square
    units = [(x,pow(x,-1,square)) for x in range(square) if gcd(x,n)==1]
    counts, mixed = [], []
    for j in range(n):
        c = 1+n*j
        counts.append(direct_count(n,c,units)-2*n*inside(n,2*c))
        mixed.append(2*int(inside(n,2*c*lift)))
    result = [(value-counts[0])/n for value in counts]
    remainder = [value-private for value,private in zip(result,mixed)]
    hidden_freq = 2*u%n
    spectrum = exact_spectrum(result)
    remainder_spectrum = exact_spectrum(remainder)
    assert max(abs(value) for k,value in enumerate(remainder_spectrum) if gcd(k,n)==1) < 1e-9
    assert abs(spectrum[hidden_freq]) > .55
    assert 1 < gcd(hidden_freq-2,n) < n
    print(n, n, round(abs(spectrum[hidden_freq]),6), max(map(abs,result)))


if __name__ == "__main__":
    for factors in [(3,5),(5,7),(11,13),(9,25)]:
        analyze(factors)
