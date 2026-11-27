"""Check stationary blocks and the hidden upper-digit frequencies."""

from cmath import exp, pi
from math import gcd

from check_quarter_fourier import exact_spectrum


def phase(n, value):
    return exp(2j*pi*(value % n)/n)


def direct(n, c, units):
    square = n*n
    return sum(phase(square, m+c*c*inverse) for m,inverse in units)/n


def check(p, q):
    n = p*q
    square = n*n
    mixed = 1+p*(-2*pow(p,-1,q) % q)
    lift = mixed - n*((mixed*mixed-1)//n)*pow(2*mixed,-1,n)
    assert lift*lift % square == 1
    units = [(m,pow(m,-1,square)) for m in range(1,square) if gcd(m,n) == 1]
    values = []
    for j in range(n):
        c = 1+n*j
        normalized = direct(n,c,units)
        public = phase(square,2*c)+phase(square,-2*c)
        private = phase(square,2*c*lift)+phase(square,-2*c*lift)
        assert abs(normalized-public-private) < 1e-8
        values.append(normalized-public)
    actual = {frequency for frequency, coeff in enumerate(exact_spectrum(values))
              if abs(coeff) > 1e-6}
    assert actual == {2*mixed % n, -2*mixed % n}
    assert all(1 < gcd(frequency-2,n) < n for frequency in actual)
    print(n, n, sorted(actual))


if __name__ == "__main__":
    for p,q in [(3,5),(3,7),(5,7),(9,5),(9,7)]:
        check(p,q)
