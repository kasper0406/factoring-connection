"""Run the coarse decoder on public, directly summed Salié values."""

from fractions import Fraction
from math import acos, cos, gcd, pi, sqrt

from sympy import jacobi_symbol

from check_disagreement_decoder import recover_factor
from check_salie_quarter import phase


def public_rate_oracle(n):
    gauss = sqrt(n) * (1 if n % 4 == 1 else 1j)
    terms = [(m, pow(m, -1, n), int(jacobi_symbol(m, n)))
             for m in range(1, n) if gcd(m, n) == 1]
    inv4 = pow(4, -1, n)

    def rate(h):
        if gcd(h, n) != 1:
            return Fraction(2 * min(h, n-h), n)
        b = h*h*inv4 % n
        total = sum(symbol * phase(n, m+b*inverse)
                    for m, inverse, symbol in terms) / gauss
        coarse = round(total.real * 10000) / 10000
        hidden_cos = coarse / 2 - cos(2*pi*h/n)
        estimate = acos(max(-1, min(1, hidden_cos))) / pi
        return Fraction(round(100 * estimate), 100)

    return rate


if __name__ == "__main__":
    for modulus in (15, 77, 143, 323, 1147, 2021, 10403):
        factor = recover_factor(modulus, public_rate_oracle(modulus))
        assert 1 < factor < modulus and modulus % factor == 0
        print(modulus, factor)
