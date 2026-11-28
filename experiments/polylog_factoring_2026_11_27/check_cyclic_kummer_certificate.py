"""Check a cyclotomic Kummer splitting certificate from an unfactored odd input."""

from dataclasses import dataclass
from math import gcd, prod
from sympy import Poly, isprime, primerange, symbols


@dataclass(frozen=True)
class Certificate:
    order: int
    scale: int
    shift: int
    norm: int


def certificate(n, order, auxiliary_prime, root_lift):
    """The auxiliary prime and root of exact order modulo its square are public."""
    if n <= 1 or n % 2 == 0 or order < 3 or not isprime(order):
        raise ValueError("expected odd N and odd prime order")
    if (not isprime(auxiliary_prime) or (auxiliary_prime-1) % order
            or gcd(n, order*auxiliary_prime) != 1):
        raise ValueError("expected an unramified input and a split auxiliary prime")
    square = auxiliary_prime**2
    if pow(root_lift, order, square) != 1 or root_lift % auxiliary_prime == 1:
        raise ValueError("expected a primitive root of the stated order")
    scale = order**3*n
    r = (root_lift + pow(scale,-1,square)*(1-auxiliary_prime) + square//2) % square - square//2
    h = scale*r-1
    # Homogeneous cyclotomic polynomial Phi_order(h, scale).
    norm = sum(h**j*scale**(order-1-j) for j in range(order))
    return Certificate(order, scale, r, norm)


def character_terms(cert, factors, value):
    if prod(t**e for t,e in factors.items()) != cert.norm or gcd(cert.norm,value) != 1:
        raise ValueError("expected full norm factors and a coprime argument")
    terms = []
    for t, exponent in factors.items():
        if not isprime(t) or (t-1) % cert.order:
            raise ValueError("invalid prime factor of the cyclotomic norm")
        zeta = (cert.shift - pow(cert.scale, -1, t)) % t
        target = pow(value, (t-1)//cert.order, t)
        power = next(j for j in range(cert.order) if pow(zeta,j,t) == target)
        terms.append(exponent*power % cert.order)
    return tuple(terms)


def _power(base, exponent, polynomial):
    result = base.one
    while exponent:
        if exponent % 2:
            result = (result*base).rem(polynomial)
        base = (base*base).rem(polynomial)
        exponent //= 2
    return result


def check_field_character(cert, factors, prime):
    ell, scale, r = cert.order, cert.scale, cert.shift
    z = symbols('z')
    def poly(expression):
        return Poly(expression, z, modulus=prime)
    modulus = poly(sum(z**i for i in range(ell)))
    degree = next(f for f in range(1,ell) if pow(prime,f,ell) == 1)
    alpha = poly(1)
    for a in range(1,ell):
        conjugate = poly(1+scale*(z**a-r))
        alpha = (alpha*_power(conjugate, pow(a,-1,ell), modulus)).rem(modulus)
    value = _power(alpha,(prime**degree-1)//ell,modulus)
    exponent = sum(character_terms(cert,factors,prime)) * degree % ell
    assert value == _power(poly(z),exponent,modulus)


def check():
    n = 91
    cert = certificate(n, 5, 11, 9)
    factors = {11: 1, 61: 1, 24955108781231: 1}
    assert cert == Certificate(5, 11375, -1, 16744877992206001)
    assert cert.norm % 121 and gcd(cert.norm, cert.scale) == 1
    assert all(isprime(t) for t in factors)
    assert 7*13 == n
    assert character_terms(cert,factors,7) == (4,2,4)
    assert character_terms(cert,factors,13) == (2,0,3)
    assert character_terms(cert,factors,2) == (2,3,3)
    checks = 0
    for prime in primerange(2,350):
        if gcd(prime,cert.order*cert.norm) == 1:
            check_field_character(cert,factors,prime)
            checks += 1
    print(f"N={n}, cyclic degree {cert.order}, conductor {cert.norm}, "
          f"{checks} Frobenius comparisons")
    print("Local exponents at 7 and 13:", character_terms(cert,factors,7),
          character_terms(cert,factors,13))


if __name__ == "__main__":
    check()
