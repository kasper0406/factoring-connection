"""Evaluate the damped fourth harmonic theta series at a rational phase."""

from fractions import Fraction
from math import isqrt

from mpmath import mp

from circle_tangent_arcs import tangent_direction


def theta_sample(n, phase, bits=40):
    """Use heat time 1/n and target absolute error 2**(-bits)."""
    if n < 1 or bits < 1:
        raise ValueError("n and bits must be positive")
    phase = Fraction(phase) % 1
    a, q = tangent_direction(phase.numerator, phase.denominator, isqrt(n))
    offset = phase - Fraction(a, q)
    with mp.workprec(96 + 12 * n.bit_length() + 2 * bits):
        w = 1 / mp.mpf(n) - 2j * mp.pi * offset.numerator / offset.denominator
        decay = mp.pi**2 / (q*q*w)
        if q % 2:
            modulus, parity = q, None
            twist = -pow(4*a, -1, q)
            scale = q if q % 4 == 1 else -q
        elif q % 4 == 0:
            modulus, parity = q, 0
            twist = -pow(a, -1, q)
            scale = 2j*q if a % 4 == 1 else -2j*q
        else:
            modulus, parity = q//2, 1
            twist = -pow(8*a, -1, modulus)
            scale = 4*modulus if modulus % 4 == 1 else -4*modulus
        radius = int(mp.ceil(mp.sqrt((bits + 8*n.bit_length() + 96) /
                                     decay.real)))
        moments = [mp.mpc(0) for _ in range(3)]
        for h in range(radius + 1):
            if parity is not None and h % 2 != parity:
                continue
            index = h//2 if q % 4 == 0 else h
            residue = (twist * index**2) % modulus
            term = mp.exp(-decay*h*h + 2j*mp.pi*residue/modulus)
            if h:
                term *= 2
            for j in range(3):
                moments[j] += term * h**(2*j)
        s0, s2, s4 = moments
        return scale * mp.pi**5 * (s0*s4 - 3*s2**2) / (2*q**6*w**5)
