"""Factor from smooth primitive Shanks lifts, using the witnesses throughout."""

from math import gcd, isqrt
from random import SystemRandom

from sympy import primerange

from check_shanks_split_seeds import split_seed
from shanks_conductor_relations import binary_norm, kernel_basis, shortest_rational_lift
from shanks_witness_fusion import ONE, ROOTS, coefficients_mod, conjugate, gamma, mul


def primary_generator(prime):
    """Choose a norm-prime Eisenstein generator congruent to 1 modulo 3."""
    for b in range(1, isqrt(4 * prime // 3) + 1):
        disc = 4 * prime - 3 * b*b
        root = isqrt(disc)
        if root*root == disc and (b + root) % 2 == 0:
            base = ((b + root)//2, b)
            for rotation in ROOTS:
                a, c = mul(base, rotation)
                for sign in (1, -1):
                    if (sign*a) % 3 == 1 and c % 3 == 0:
                        return sign*a, sign*c
    raise ValueError('expected a prime congruent to 1 modulo 3')


def divide_if_integral(value, factor, norm):
    a, b = mul(value, conjugate(factor))
    if a % norm or b % norm:
        return None
    return a // norm, b // norm


def smooth_row(pair, base):
    norm = binary_norm(pair)
    if norm % 3 == 0:
        return None
    beta = pair[0], -3*pair[1]
    exponents = []
    for prime, generator in base:
        power = 0
        if norm % prime == 0:
            orientation = generator
            sign = 1
            if divide_if_integral(beta, orientation, prime) is None:
                orientation = conjugate(orientation)
                sign = -1
            while norm % prime == 0:
                beta = divide_if_integral(beta, orientation, prime)
                assert beta is not None
                norm //= prime
                power += sign
        exponents.append(power)
    if norm != 1:
        return None
    assert beta in (ONE, (-1, 0))
    return tuple(exponents)


def modular_power(value, exp, n):
    if exp < 0:
        value, exp = conjugate(value), -exp
    answer = ONE
    while exp:
        if exp & 1:
            answer = tuple(a % n for a in mul(answer, value))
        value = tuple(a % n for a in mul(value, value))
        exp //= 2
    return answer


def mixed_factor(n, rows, witnesses, class_generators, dep):
    """Return a factor when the dependency's two CRT cube-root labels differ."""
    z = ONE
    for u, e in zip(witnesses, dep, strict=True):
        z = tuple(a % n for a in mul(z, modular_power(u, e, n)))
    for i, generator in enumerate(class_generators):
        total = sum(e * row[i] for row, e in zip(rows, dep, strict=True))
        assert total % 3 == 0
        z = tuple(a % n for a in mul(z, modular_power(generator, -(total//3), n)))
    assert modular_power(z, 3, n) == ONE
    for root in ROOTS:
        factor = gcd(n, z[0] - root[0], z[1] - root[1])
        if 1 < factor < n:
            return factor
    return None


def factor_smooth_batch(n, bound, *, extra_rows=2, rng=None):
    """On a balanced semiprime >3, collect batches until a verified split occurs."""
    if n <= 9 or n % 2 == 0 or bound < 7 or extra_rows < 1:
        raise ValueError('expected an odd balanced semiprime >9, bound≥7 and extra_rows≥1')
    rng = rng if rng is not None else SystemRandom()
    split_primes = [int(r) for r in primerange(7, bound + 1) if r % 3 == 1]
    for prime in (3, *split_primes):
        if n % prime == 0:
            return prime, {'source': 'factor base', 'seeds': 0}
    base = [(r, primary_generator(r)) for r in split_primes]
    class_generators = [coefficients_mod(gamma(g), n) for _, g in base]
    stats = {'seeds': 0, 'smooth_rows': 0, 'batches': 0}
    while True:
        rows, witnesses = [], []
        while len(rows) < len(base) + extra_rows:
            x = rng.randrange(n)
            stats['seeds'] += 1
            divisor, t = split_seed(n, x)
            if t is not None:
                divisor, pair = shortest_rational_lift(n, t)
                if pair is not None:
                    row = smooth_row(pair, base)
                    if row is not None:
                        rows.append(row)
                        witnesses.append(coefficients_mod(gamma((x, -1)), n))
                        stats['smooth_rows'] += 1
            if divisor is not None and 1 < divisor < n:
                return divisor, dict(stats, source='seed gcd')
        deps = kernel_basis(list(zip(*rows)), len(rows))
        stats['batches'] += 1
        stats['nullity'] = len(deps)
        for dep in deps:
            divisor = mixed_factor(n, rows, witnesses, class_generators, dep)
            if divisor is not None:
                return divisor, dict(stats, source='mixed dependency')
