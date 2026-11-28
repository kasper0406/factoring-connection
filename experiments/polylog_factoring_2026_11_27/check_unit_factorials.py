"""Compare Bhargava factorial valuations for l-adic units and all integers.

For degrees 0,...,D-1, a decrease in determinant by a factor K permits an
increase in radius of K**(2/(D*(D-1))) in the triangular coefficient volume.
The finite support is l<=D. This measures the extra unit condition alone.
"""

from fractions import Fraction
from math import exp, factorial, fsum, isqrt, log, prod


def primes_through(limit):
    sieve = bytearray([1])*(limit+1)
    for p in range(2, isqrt(limit)+1):
        if sieve[p]:
            sieve[p*p::p] = bytes((limit-p*p)//p+1)
    return [p for p in range(2, limit+1) if sieve[p]]


def valuation(value, p):
    if not value:
        raise ValueError("zero has no finite valuation")
    count = 0
    while value % p == 0:
        value //= p
        count += 1
    return count


def unit_powers(p, degree):
    """Yield (phi(p^k), p^k) while their floor difference can contribute."""
    a, b = p-1, p
    while a <= degree:
        yield a, b
        a, b = a*p, b*p


def floor_sum(d, step):
    """Sum floor(j/step) for j=0,...,d-1."""
    q, rem = divmod(d, step)
    return step*q*(q-1)//2 + q*rem


def determinant_valuation(d, p):
    return sum(floor_sum(d, a)-floor_sum(d, b) for a,b in unit_powers(p, d-1))


def radius_multiplier(d):
    log_det = fsum(determinant_valuation(d,p)*log(p) for p in primes_through(d))
    return exp(2*log_det/(d*(d-1)))


def check_orderings():
    count = 0
    for p in primes_through(53):
        nodes = [x for x in range(1, 130) if x % p]
        for j in range(33):
            ordinary = valuation(factorial(j), p)
            unit = sum(j//a for a,_ in unit_powers(p, j))
            # The next unit attains the minimum. Compare other units as well.
            values = [valuation(prod(x-y for y in nodes[:j]),p) for x in nodes[j:j+20]]
            assert values[0] == min(values) == unit
            assert determinant_valuation(j+1,p)-determinant_valuation(j,p) == unit-ordinary
            count += 1
    return count


def check_determinants():
    count = 0
    for d in list(range(2, 128))+[255, 512, 1023, 4096]:
        for p in primes_through(139):
            total = 0
            for a,b in unit_powers(p, d-1):
                contribution = floor_sum(d,a)-floor_sum(d,b)
                main = Fraction(d*d*(b-a), 2*a*b)
                # Exact remainder and integrated derivative bounds, respectively.
                assert abs(contribution-main) <= min(Fraction(b,8), Fraction(d,2*(p-1)))
                total += contribution
            assert total == determinant_valuation(d,p)
            main = Fraction(d*d,2*(p-1)**2)
            if p <= d:
                assert abs(total-main) <= Fraction(5*d,4*(p-1))
            else:
                assert total == 0
            count += 1
    return count


if __name__ == '__main__':
    print('Local orderings and determinant increments:', check_orderings())
    print('Exact degree/prime determinant bounds:', check_determinants())
    for dimension in (8, 32, 128, 512, 2048, 8192):
        print(f'D={dimension:4d}   extra radius multiplier {radius_multiplier(dimension):.8f}')
