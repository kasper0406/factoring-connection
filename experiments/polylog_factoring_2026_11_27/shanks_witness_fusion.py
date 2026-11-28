"""Fuse cubic splitting witnesses to a short rational parameter or a factor."""

from fractions import Fraction
from functools import lru_cache
from math import gcd, isqrt, lcm

from sympy import Poly, primerange, symbols

from check_shanks_split_seeds import split_seed
from shanks_conductor_relations import (
    binary_norm, eisenstein_product as mul, exponent, rational_trivial,
    shortest_rational_lift,
)


ONE = (1,0)
OMEGA = (0,1)
ROOTS = (ONE, OMEGA, (-1,-1))


def conjugate(z):
    a, b = z
    return a-b, -b


def div(a, b):
    norm = mul(b,conjugate(b))[0]
    return tuple(Fraction(v,norm) for v in mul(a,conjugate(b)))


def gamma(beta):
    return div(beta,conjugate(beta))


def power(z, k):
    value = ONE
    while k:
        if k % 2:
            value = mul(z,value)
        z = mul(z,z)
        k //= 2
    return value


@lru_cache(maxsize=512)
def norm_one_cube_root(z):
    """Use T = h + 1/h and T^3-3T = z+1/z to avoid factoring large norms."""
    if z == ONE or z == (-1,0):
        return z
    trace = 2*z[0]-z[1]
    variable = symbols('v')
    polynomial = Poly(variable**3-3*variable-trace,variable,domain='QQ')
    for T in sorted(polynomial.ground_roots()):
        T = Fraction(T)
        if T*T != 1:
            h = div((z[0]+T,z[1]),(T*T-1,0))
            if power(h,3) == z:
                return h
    return None


def oriented_generator(prime, root):
    radius = 2*isqrt(prime)+2
    for a in range(-radius,radius+1):
        for b in range(-radius,radius+1):
            if a*a-a*b+b*b == prime and (a+b*root) % prime == 0:
                return a,b
    raise AssertionError("split prime has no generator in the norm disk")


def short_beta(pairs, vector, limit):
    beta = ONE
    for p in primerange(7,limit+1):
        if p % 3 != 1:
            continue
        valuations = [exponent(binary_norm(pair),p) for pair in pairs]
        if not any(valuations):
            continue
        root = next(r for r in range(p) if (r*r+r+1) % p == 0)
        signed = sum(e*v*(1 if (a-3*b*root) % p == 0 else -1)
                     for (a,b),e,v in zip(pairs,vector,valuations)) % 3
        if signed:
            pi = oriented_generator(p,root)
            beta = mul(beta,pi if signed == 1 else conjugate(pi))
    return beta


def coefficients_mod(z, n):
    den = lcm(*(Fraction(c).denominator for c in z))
    common = gcd(den,n)
    if common != 1:
        raise ValueError(f"Eisenstein denominator gcd {common}")
    return tuple(int(c*den)*pow(den,-1,n) % n for c in z)


def fuse(n, pairs, seeds, vector, limit):
    """Given valid witnessed pairs and a relation, return (factor, short_pair, root)."""
    if not all(e in (0, 1, 2) for e in vector):
        raise ValueError("expected relation coefficients 0, 1 or 2")
    product, witness = ONE, ONE
    for (a,b), x, e in zip(pairs,seeds,vector,strict=True):
        if (b*(x**3-3*x-1)-a*x*(x+1)) % n or gcd(x*(x+1)*(x*x+x+1),n) != 1:
            raise ValueError("invalid modular root certificate")
        product = mul(product,power(gamma((a,-3*b)),e))
        value = coefficients_mod(gamma((x,-1)),n)
        witness = coefficients_mod(mul(witness,power(value,e)),n)
    beta0 = short_beta(pairs,vector,limit)
    for rotation in ROOTS:
        beta = mul(beta0,rotation)
        h = norm_one_cube_root(div(product,gamma(beta)))
        if h is None:
            continue
        z = coefficients_mod(mul(witness,conjugate(h)),n)
        if beta[1] == 0:
            for root in ROOTS:
                common = gcd(n,z[0]-root[0],z[1]-root[1])
                if 1 < common < n:
                    return common,None,None
            assert z in [coefficients_mod(root,n) for root in ROOTS]
            return None,None,None

        # x=(z*omega_bar-omega)/(z-1). Any nonzero component gcd already splits n.
        denominator = (z[0]-1,z[1])
        den_norm = mul(denominator,conjugate(denominator))[0]
        common = gcd(den_norm,n)
        if 1 < common < n:
            return common,None,None
        if common == n:
            for root in ROOTS[1:]:
                rotated = coefficients_mod(mul(z,root),n)
                denominator = (rotated[0]-1,rotated[1])
                common = gcd(mul(denominator,conjugate(denominator))[0],n)
                if 1 < common < n:
                    return common,None,None
                if common == 1:
                    z = rotated
                    break
            else:
                raise AssertionError("three diagonal rotations cover the two local components")
        a,b = 3*beta[0],-beta[1]
        content = gcd(a,b)*(1 if b>0 else -1)
        pair = (a//content,b//content)
        numerator = mul(z,conjugate(OMEGA))
        x, other = coefficients_mod(div((numerator[0],numerator[1]-1),denominator),n)
        assert other == 0
        assert (pair[1]*(x**3-3*x-1)-pair[0]*x*(x+1)) % n == 0
        return None,pair,x
    raise ValueError("combination has ramification above the cutoff")


def check():
    for n,seeds,vector in [(1247,(63,66),(1,1)),(2813,(78,111),(2,1))]:
        pairs = [(split_seed(n,x)[1],1) for x in seeds]
        factor,pair,x = fuse(n,pairs,seeds,vector,7)
        assert factor is None and binary_norm(pair) <= 9*7
        print(f"N={n}, fused seeds {seeds} into parameter {pair} with public root {x}")

    trivial = splits = diagonal = 0
    for n in (713,1247,2813):
        for seed in range(n):
            _, t = split_seed(n,seed)
            if t is None:
                continue
            divisor,pair = shortest_rational_lift(n,t)
            if divisor or not rational_trivial(pair):
                continue
            trivial += 1
            factor,short,x = fuse(n,[pair],[seed],[1],3)
            assert short is None and x is None
            if factor:
                assert 1 < factor < n and n % factor == 0
                splits += 1
            else:
                diagonal += 1
    assert (trivial,splits,diagonal) == (135,90,45)
    print(f"Previously trivial fields: {splits}/{trivial} further gcd splits, {diagonal} diagonal roots")


if __name__ == "__main__":
    check()
