"""Full-period character moment, modular recovery, and finite Fourier checks."""

from functools import cache
from itertools import combinations
from math import gcd, isqrt

from sympy import isprime, kronecker_symbol, primerange, primitive_root


@cache
def character(delta):
    return tuple(int(kronecker_symbol(delta, a)) for a in range(abs(delta)))


def cubic_constant(chi):
    f = len(chi)
    # f times the third generalized Bernoulli number of the short character.
    return sum(weight * a * (a*a-f*f) for a, weight in enumerate(chi))


def direct_moments(n, chi):
    """Input-only reference evaluation over one complete induced period."""
    f = len(chi)
    moments = [0, 0, 0, 0]
    for a in range(1, f*n):
        weight = chi[a % f]
        if weight and gcd(a, n) == 1:
            moments[0] += weight
            moments[1] += weight*a
            moments[2] += weight*a*a
            moments[3] += weight*a*a*a
    return tuple(moments)


def split_from_cubic_digit(n, chi, reading):
    """Recover the positive-character prime from M3 mod |C|N² on the balanced promise."""
    c = cubic_constant(chi)
    if not c:
        raise ValueError("expected a nonzero cubic constant")
    residue = reading % (abs(c)*n*n)
    scaled, remainder = divmod(residue, abs(c)*n)
    if remainder:
        raise ValueError("inconsistent cubic moment")
    gap_mod = ((scaled if c > 0 else -scaled) - 1) % n
    for lift in range(-2, 2):
        gap = gap_mod + lift*n
        if 2*abs(gap) >= 3*n:
            continue
        radicand = gap*gap + 4*n*n
        total = isqrt(radicand)
        if total*total != radicand or (total-gap) % 2:
            continue
        r2 = (total-gap)//2
        r = isqrt(r2)
        if 1 < r < n and r*r == r2 and n % r == 0 and chi[r % len(chi)] == 1:
            return r
    raise ValueError("the reading does not give a split prime")


def check_moments():
    cases = 0
    fixtures = {}
    for delta in (-3, -4, -7, -8, -20, -68):
        chi = character(delta)
        c = cubic_constant(chi)
        assert c != 0
        for p, q in combinations(primerange(3, 36), 2):
            n = p*q
            if q >= 2*p or chi[n % len(chi)] != -1:
                continue
            r, d = (p, q) if chi[p % len(chi)] == 1 else (q, p)
            moments = direct_moments(n, chi)
            assert moments == (0, 0, 0, n*c*(1-r*r)*(1+d*d))
            digit = moments[3] % (abs(c)*n*n)
            assert split_from_cubic_digit(n, chi, digit) == r
            if gcd(c, n) == 1:
                gap_mod = ((moments[3] % (n*n))//n * pow(c, -1, n) - 1) % n
                assert gap_mod == (d*d-r*r) % n
            if (delta, n) in ((-4, 15), (-8, 77), (-68, 437)):
                fixtures[n] = (delta, c, moments[3], digit, r)
            cases += 1
    return cases, fixtures


def roots_of_unity(period):
    """Specialize complex-root identities in a prime field for exact tests."""
    for multiplier in range(1, 10001):
        modulus = multiplier*period+1
        if isprime(modulus):
            root = pow(int(primitive_root(modulus)), multiplier, modulus)
            return modulus, [pow(root, j, modulus) for j in range(period)]
    raise ValueError("test field search exhausted")


def cot_kernel_sum(weights, stride, powers, modulus):
    # cot(theta) = i*(w+1)/(w-1) with w = exp(2i theta); i cancels in the ratio.
    total = 0
    for j, weight in enumerate(weights):
        if weight:
            w = powers[j*stride]
            total += weight*(w+1)*pow(w-1, -1, modulus)
    return total % modulus


def check_fourier(delta, p, q):
    n, chi = p*q, character(delta)
    f = len(chi)
    assert chi[n % f] == -1
    r, d = (p, q) if chi[p % f] == 1 else (q, p)
    m = f*n
    modulus, powers = roots_of_unity(m)
    weights = [chi[a % f] if gcd(a, n) == 1 else 0 for a in range(m)]
    supported = [(a, weight) for a, weight in enumerate(weights) if weight]
    gauss = sum(w*powers[a*n] for a, w in enumerate(chi)) % modulus
    assert gauss
    for k in range(m):
        observed = sum(w*powers[(k*a) % m] for a, w in supported) % modulus
        ramanujan = (p-1 if k % p == 0 else -1)*(q-1 if k % q == 0 else -1)
        assert observed == -chi[k % f]*gauss*ramanujan % modulus
        if gcd(k, n) == 1:
            assert observed == -chi[k % f]*gauss % modulus
    base = cot_kernel_sum(chi, n, powers, modulus)
    induced = cot_kernel_sum(weights, 1, powers, modulus)
    assert base and induced == ((n+r-d-1)*base) % modulus
    return m


if __name__ == "__main__":
    count, fixtures = check_moments()
    print(f"{count} complete periods yielded the positive-character prime from one cubic reading")
    print(f"Fixtures N -> (discriminant, cubic constant, cubic moment, modular reading, prime): {fixtures}")
    total = sum(check_fourier(*case) for case in ((-4, 3, 5), (-7, 11, 13), (-8, 7, 11)))
    print(f"{total} finite-field Fourier coefficients and three cotangent sums verified")
