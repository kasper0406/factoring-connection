"""Construct a packed ternary sphere with an explicitly factored conductor."""

from dataclasses import dataclass
from math import gcd, isqrt, prod

from sympy import factorint, isprime, jacobi_symbol

if __package__:
    from .fixed_field_packing import Parameters, packing_params
else:
    from fixed_field_packing import Parameters, packing_params


def certified_prime(n):
    """Prove primality; beyond the deterministic MR range use Lucas's criterion."""
    if not isprime(n):
        return False
    if n < 2**64:
        return True
    factors = {int(p): e for p, e in factorint(n - 1).items()}
    if prod(p**e for p, e in factors.items()) != n - 1:
        return False
    for p in factors:
        if not certified_prime(p):
            return False
        for a in range(2, n):
            if pow(a, n - 1, n) != 1:
                return False
            if gcd(pow(a, (n - 1) // p, n) - 1, n) == 1:
                break
        else:
            return False
    return True


def _norm(v):
    return v[0] ** 2 + v[1] ** 2


def _reduced_line(n, a, b):
    # Preimage of the projective line spanned by (a,b) modulo n.
    g = gcd(a, n)
    m = n // g
    slope = b * pow(a // g, -1, m) % m if m > 1 else 0
    u, v = (g, slope), (0, m)
    while True:
        if _norm(v) < _norm(u):
            u, v = v, u
        square = _norm(u)
        dot = u[0] * v[0] + u[1] * v[1]
        k = (2 * dot + square) // (2 * square)
        v = (v[0] - k * u[0], v[1] - k * u[1])
        if _norm(v) >= square:
            return u, v


def line_octant(n, a, b, height):
    """All primitive (i,j), 0<j<i, of norm in (height²/2,height²] on this line."""
    if n < 2 or height < 1 or gcd(a, b, n) != 1:
        raise ValueError("expected a primitive projective direction and a positive height")
    u, v = _reduced_line(n, a, b)
    square = _norm(u)
    dot = u[0] * v[0] + u[1] * v[1]
    bound = height * height
    last = isqrt(square * bound) // n
    for j in range(-last, last + 1):
        if j == 0:
            # All other multiples of u are imprimitive.
            indices = (-1, 1)
        else:
            radius = isqrt(square * bound - (j * n) ** 2)
            lo = -((j * dot + radius) // square)
            hi = (radius - j * dot) // square
            indices = range(lo, hi + 1)
        for i in indices:
            x, y = i * u[0] + j * v[0], i * u[1] + j * v[1]
            if 0 < y < x and bound < 2 * (x*x + y*y) <= 2 * bound and gcd(x, y) == 1:
                yield x, y


@dataclass(frozen=True)
class FactoredPacking:
    parameters: Parameters
    gaussian_primes: tuple
    primes: tuple
    t: int
    height: int
    lines: int

    @property
    def norm_factorization(self):
        """Return the full factorization of M, including its short base norm."""
        base = self.parameters.M // self.parameters.T**2
        factors = {}
        divisor = 2
        while divisor * divisor <= base:
            while base % divisor == 0:
                factors[divisor] = factors.get(divisor, 0) + 1
                base //= divisor
            divisor += 1 if divisor == 2 else 2
        if base > 1:
            factors[base] = 1
        for prime in self.primes:
            factors[prime] = 2
        return factors


def factored_packing(n, x, rng, *, prime_test=certified_prime):
    """Input-only construction for a mixed odd semiprime. prime_test must be exact."""
    a = x*x + 1
    if n < 3 or n % 2 == 0 or x < 0 or jacobi_symbol(-a, n) != -1:
        raise ValueError("expected an odd input and a nonnegative mixed parameter")
    bins = 8 * (x + 1)
    height = (1 + isqrt(n * bins)) * n.bit_length()
    trials = 0
    while True:
        for _ in range(16 * height.bit_length() ** 2):
            u, v = rng.randrange(n), rng.randrange(n)
            if gcd(u, v, n) != 1:
                continue
            trials += 1
            seen = {}
            for point in line_octant(n, u, v, height):
                p = _norm(point)
                if gcd(p, 2 * a * n) != 1 or not prime_test(p):
                    continue
                key = bins * point[1] // point[0], point[0] % 2
                if key not in seen:
                    seen[key] = point
                    continue
                other = seen[key]
                real = other[0] * point[0] + other[1] * point[1]
                imag = abs(other[0] * point[1] - other[1] * point[0])
                t, remainder = divmod(imag, n)
                assert remainder == 0 and t >= 2 and t % 2 == 0
                params = packing_params(n, x, real, t)
                primes = (_norm(other), p)
                assert primes[0] != primes[1] and prod(primes) == params.T
                return FactoredPacking(params, (other, point), primes, t, height, trials)
        # Cover the finitely many small-input exceptions to the asymptotic prime count.
        height *= 2
