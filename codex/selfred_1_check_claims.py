"""Finite sanity checks for selfred_phase22_report.md; these are not asymptotic proofs.

Uses only Python's standard library. Run: python3 check_claims.py
"""

from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd, isqrt, lcm, prod
from random import Random


def factor_small(m):
    result = {}
    r = 2
    while r * r <= m:
        while m % r == 0:
            result[r] = result.get(r, 0) + 1
            m //= r
        r += 1
    if m > 1:
        result[m] = result.get(m, 0) + 1
    return result


def prime(m):
    return m >= 2 and len(factor_small(m)) == 1 and factor_small(m).get(m) == 1


def radical(m):
    return prod(factor_small(m))


def parity_kernel(m):
    return prod(r for r, e in factor_small(m).items() if e % 2)


def divisors_count(m):
    return prod(e + 1 for e in factor_small(m).values())


def horner(coefficients, x):
    answer = 0
    for c in reversed(coefficients):
        answer = answer * x + c
    return answer


def check_height_free():
    rng = Random(2201)
    checked = 0
    for p in (13, 17, 23, 31, 43, 59, 79):
        samples = [q for q in range(p + 1, 2 * p) if prime(q)]
        for d in range(1, min(5, len(samples))):
            for _ in range(30):
                qs = sorted(rng.sample(samples, rng.randrange(d + 1, len(samples) + 1)))
                if checked % 3 == 0:
                    coeffs = [-1] + [0] * (d - 1) + [1]
                else:
                    coeffs = [rng.randrange(-10000, 10001) for _ in range(d)]
                    coeffs.append(rng.randrange(1, 10001))
                coeffs = [c * 12 * p**2 for c in coeffs]
                content = reduce(gcd, map(abs, coeffs))
                common = reduce(gcd, (abs(horner(coeffs, p * q)) for q in qs))
                assert common > 0  # a nonzero degree-d polynomial has at most d roots
                remainder = common // gcd(common, content)
                while remainder % p == 0:
                    remainder //= p
                block = min((qs[i:i + d + 1] for i in range(len(qs) - d)),
                            key=lambda v: v[-1] - v[0])
                diameter = block[-1] - block[0]
                delta = prod(v - u for u, v in combinations(block, 2))
                assert diameter <= (d * p) // (len(qs) - d)
                assert delta % remainder == 0
                assert remainder <= diameter ** (d * (d + 1) // 2)
                assert all(r <= diameter for r in factor_small(remainder))
                checked += 1
    return checked


def interpolate_signature(m):
    b = m.bit_length() - 1
    ys = [1] + [divisors_count(m**j) for j in range(1, b + 1)]
    coefficients = [Fraction(0) for _ in range(b + 1)]
    for j, value in enumerate(ys):
        basis = [Fraction(1)]
        denominator = 1
        for k in range(b + 1):
            if k == j:
                continue
            next_basis = [Fraction(0)] * (len(basis) + 1)
            for i, a in enumerate(basis):
                next_basis[i] -= k * a
                next_basis[i + 1] += a
            basis = next_basis
            denominator *= j - k
        for i, a in enumerate(basis):
            coefficients[i] += value * a / denominator
    while len(coefficients) > 1 and not coefficients[-1]:
        coefficients.pop()
    signature = []
    for e in range(1, b + 1):
        while len(coefficients) > 1 and horner(coefficients, -Fraction(1, e)) == 0:
            # Divide by 1 + e*T, starting with the highest coefficient.
            quotient = [Fraction(0)] * (len(coefficients) - 1)
            residual = coefficients[:]
            for i in reversed(range(1, len(residual))):
                quotient[i - 1] = residual[i] / e
                residual[i - 1] -= quotient[i - 1]
            assert residual[0] == 0
            coefficients = quotient
            signature.append(e)
    assert coefficients == [1]
    return sorted(signature)


def check_oracles():
    for m in range(1, 1001):
        expected = factor_small(m)
        x, rs = m, []
        while x > 1:
            r = radical(x)
            rs.append(r)
            x //= r
        calls = len(rs)
        assert calls <= m.bit_length() - 1
        rs.append(1)
        blocks = [rs[j - 1] // rs[j] for j in range(1, len(rs))]
        assert prod(b**j for j, b in enumerate(blocks, 1)) == m
        assert all(radical(b) == b for b in blocks)
        assert all(gcd(a, b) == 1 for a, b in combinations(blocks, 2))
        x, svalues = m, []
        while x > 1:
            s = parity_kernel(x)
            svalues.append(s)
            u = isqrt(x // s)
            assert u * u * s == x
            x = u
        assert lcm(*svalues) == radical(m)
        if m > 1:
            b = m.bit_length() - 1
            r = next(r for r in range(b + 2, 2 * (b + 2) + 1) if prime(r))
            d = divisors_count(m**(r - 1))
            exponent = 0
            while d % r == 0:
                exponent += 1
                d //= r
            mu_from_d = (-1)**exponent if d == 1 else 0
            mu = (-1)**len(expected) if all(e == 1 for e in expected.values()) else 0
            assert mu_from_d == mu
        if m <= 200:
            assert interpolate_signature(m) == sorted(expected.values())
    return 1000


def reduced_forms(discriminant):
    assert discriminant < 0 and discriminant % 4 in (0, 1)
    forms = []
    for a in range(1, isqrt((-discriminant) // 3) + 1):
        for b in range(-a, a + 1):
            if (b*b - discriminant) % (4*a):
                continue
            c = (b*b - discriminant) // (4*a)
            if a > c or gcd(gcd(a, b), c) != 1:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            forms.append((a, b, c))
    return forms


def legendre(a, p):
    a %= p
    if not a:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def twist_count(d, p):
    # d*y^2 = x^3-x, including the point at infinity.
    return p + 1 + legendre(d, p) * sum(legendre(x**3 - x, p) for x in range(p))


def check_orders_and_twists():
    p, q = 5, 7
    n = p * q
    units = [d for d in range(1, n) if gcd(d, n) == 1]
    good = 0
    for d1 in units:
        for d2 in units:
            s1 = (legendre(d1, p), legendre(d1, q))
            s2 = (legendre(d2, p), legendre(d2, q))
            independent = s1 != (1, 1) and s2 != (1, 1) and s1 != s2
            if independent:
                counts = [twist_count(d, p) * twist_count(d, q)
                          for d in (1, d1, d2, d1*d2)]
                assert sum(counts) == 4 * (p + 1) * (q + 1)
                good += 1
    assert Fraction(good, len(units)**2) == Fraction(3, 8)
    for p, q in ((5, 7), (13, 17), (17, 29), (19, 23)):
        n, s = p*q, p+q
        assert (p*p+p+1)*(q*q+q+1) == s*s+(n+1)*s+n*n-n+1
        if n % 4 == 1:
            assert ((p-1)*(q-1) % 8 == 0) == (p % 4 == 1)
    return good, len(units)**2


def fermat_check(x):
    ps = [p for p in range(x, 2*x + 1) if prime(p)]
    eligible, success = 0, 0
    by_class = {1: [0, 0], 3: [0, 0]}
    for p, q in combinations(ps, 2):
        if p*q % 4 != 1:
            continue
        eligible += 1
        by_class[p % 4][0] += 1
        n = p*q
        a = isqrt(n) + 1
        b = isqrt(a*a - n)
        found = b*b == a*a - n and (a-b)*(a+b) == n
        if found:
            success += 1
            by_class[p % 4][1] += 1
        if (q-p)**2 <= x:
            assert found
    balanced_advantage = sum(Fraction(s, 4*t) for t, s in by_class.values())
    return {"X": x, "eligible_pairs": eligible, "one_test_factors": success,
            "advantage_over_fair_guess_raw_law": str(Fraction(success, 2 * eligible)),
            "advantage_over_best_constant_fair_class_mixture": str(balanced_advantage)}


def main():
    print("Height-free theorem finite cases:", check_height_free())
    print("Oracle identities, integers checked:", check_oracles())
    for p, q, expected in ((13, 17, 16), (7, 11, 8), (17, 29, 12), (19, 23, 20)):
        forms = reduced_forms(-4*p*q)
        assert len(forms) == expected
        print("Class-number example:", p, q, "h =", len(forms), "forms =", forms)
    good, total = check_orders_and_twists()
    print("Independent twist-sign pairs:", good, "/", total, "= 3/8")
    for x in (256, 512, 1024):
        print("Fermat finite distribution:", fermat_check(x))
    print("All finite sanity checks passed. Asymptotic assertions require the report proofs.")


if __name__ == "__main__":
    main()
