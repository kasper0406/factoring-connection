"""Run the coarse decoder on public, directly summed Salié values."""

from fractions import Fraction
from math import acos, cos, gcd, pi, sqrt

from sympy import jacobi_symbol

from check_disagreement_decoder import factor_from_phases, recover_factor
from check_salie_quarter import phase


def public_sum_oracle(n):
    gauss = sqrt(n) * (1 if n % 4 == 1 else 1j)
    terms = [(m, pow(m, -1, n), int(jacobi_symbol(m, n)))
             for m in range(1, n) if gcd(m, n) == 1]
    inv4 = pow(4, -1, n)

    def normalized_sum(h):
        b = h*h*inv4 % n
        total = sum(symbol * phase(n, m+b*inverse)
                    for m, inverse, symbol in terms) / gauss
        return total.real

    return normalized_sum


def public_rate_oracle(n):
    normalized_sum = public_sum_oracle(n)

    def rate(h):
        if gcd(h, n) != 1:
            return Fraction(2 * min(h, n-h), n)
        coarse = round(normalized_sum(h) * 10000) / 10000
        hidden_cos = coarse / 2 - cos(2*pi*h/n)
        estimate = acos(max(-1, min(1, hidden_cos))) / pi
        return Fraction(round(100 * estimate), 100)

    return rate


def recover_from_coarse_sums(n, estimate):
    def cosine(h):
        h %= n
        public = cos(2*pi*h/n)
        if gcd(h, n) != 1:
            return public
        return max(-1, min(1, estimate(h)/2 - public))

    shifts = [1 << j for j in range(n.bit_length() + 1)]
    values = [cosine(h) for h in shifts]
    anchor = next(h for h, v in zip(shifts, values) if abs(v) <= 0.81)
    phases = []
    for h, value in zip(shifts, values):
        angle = acos(value)
        if cosine(h-anchor) < cosine(h+anchor):
            angle = -angle
        phases.append(Fraction(round(4096*angle/(2*pi)), 4096) % 1)
    return factor_from_phases(n, phases)


if __name__ == "__main__":
    for modulus in (15, 77, 143, 323, 1147, 2021, 10403):
        factor = recover_factor(modulus, public_rate_oracle(modulus))
        assert 1 < factor < modulus and modulus % factor == 0
        normalized_sum = public_sum_oracle(modulus)

        def coarse_sum(h):
            return round(normalized_sum(h) * 2.5) / 2.5

        coarse_factor = recover_from_coarse_sums(modulus, coarse_sum)
        assert 1 < coarse_factor < modulus and modulus % coarse_factor == 0
        print(modulus, factor, coarse_factor)
