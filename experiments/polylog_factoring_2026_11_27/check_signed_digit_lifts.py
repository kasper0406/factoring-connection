"""Check signed factors of an unmodified binary digit polynomial and its dilations."""

from sympy import Poly, Symbol, isprime

if __package__:
    from .carry_parity import _divide, _gcd, binary_product, factor_from_carry_parity
    from .carry_polynomial import carry_factor
else:
    from carry_parity import _divide, _gcd, binary_product, factor_from_carry_parity
    from carry_polynomial import carry_factor


def dilate(bits, step):
    return sum(1 << (i*step) for i in range(bits.bit_length()) if bits >> i & 1)


def binary_derivative(bits):
    return sum(1 << (i-1) for i in range(1, bits.bit_length(), 2) if bits >> i & 1)


def main():
    x = Symbol('X')
    f = Poly(x**10+x**7-x**6+x**5+x**2-x+1, x)
    g = Poly(x**10+x**9+x**8-x**3+x+1, x)
    n = 2006801
    source = Poly.from_list(list(map(int, bin(n)[2:])), x)
    assert f*g == source and set(source.factor_list()[1]) == {(f, 1), (g, 1)}
    for factor in (f, g):
        fourth = factor.compose(Poly(x**4, x))
        assert fourth.degree() == 40 and fourth.factor_list() == (1, [(fourth, 1)])
    assert f.eval(2) == 1123 and g.eval(2) == 1787 and isprime(1123) and isprime(1787)
    assert factor_from_carry_parity(n, templates=(0,)) is None
    assert _divide(n, 1123)[1] == sum(1 << i for i in (9, 7, 5, 4, 3, 2, 0))
    assert _divide(n, 1787)[1] == sum(1 << i for i in (5, 4, 3, 2, 1))

    for t in (1, 2, 4, 8, 16, 32, 64, 128):
        b = 1 << t
        a, d = int(f.eval(b)), int(g.eval(b))
        target = a*d
        assert target == dilate(n, t) and 1 < a < d < 2*a
        assert _divide(target, a)[1] and _divide(target, d)[1]
        if t >= 2:
            # Integer bit strings below encode polynomials over F_2.
            s = (1 << t)-1
            af = (1 << (10*t)) ^ (1 << (5*t)) ^ 1
            assert a == af ^ binary_product(s, (1 << t) ^ (1 << (6*t)))
            assert _gcd(a, binary_derivative(a)) == 1
            pg = binary_product(3, 1 ^ (1 << (9*t))) ^ sum(1 << (j*t) for j in range(3, 8))
            assert d == binary_product(s, pg)
            assert _gcd(pg, binary_derivative(pg)) == 1
            assert pg.bit_length()-1 == 9*t+1
        if t <= 16:
            divisor = carry_factor(target, 0)
            assert {divisor, target//divisor} == {a, d}
            print(f"t={t}, {target.bit_length()} bits, native factor {divisor}")
    print("Eight dilations through t=128, both designated binary divisors absent in every case")


if __name__ == '__main__':
    main()
