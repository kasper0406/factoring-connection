"""Measure promises on unplanted pairs. All scores here use the known test factors."""

from math import gcd, isqrt, prod

from sympy import factorint, primerange

if __package__:
    from .check_root_scale_order_hints import prime_pairs
else:
    from check_root_scale_order_hints import prime_pairs


def local_powers(value, primes):
    powers = []
    for ell in primes:
        prime_power = 1
        while value % ell == 0:
            prime_power *= ell
            value //= ell
        powers.append(prime_power)
    return powers


def coverage(p, q, pool):
    n = p * q
    fixed = pool[:4]
    minus = [prod(local_powers(r - 1, fixed)) for r in (p, q)]
    plus = [prod(local_powers(r + 1, fixed)) for r in (p, q)]
    both = [a * b for a, b in zip(minus, plus)]
    shared = gcd(*both)
    assert shared == prod(local_powers(gcd(p*p - 1, q*q - 1), fixed))
    assert max(both) <= max(r*r - 1 for r in (p, q))
    best_four = max(prod(sorted(local_powers(r*r - 1, pool), reverse=True)[:4])
                    for r in (p, q))

    def qualifies(d):
        return d**16 > n**5

    fixed_mixed = qualifies(max(both))
    single_sign = qualifies(max(minus + plus))
    common = qualifies(shared)
    best = qualifies(best_four)
    assert max(common, single_sign) <= fixed_mixed <= best
    return common, single_sign, fixed_mixed, fixed_mixed and not (common or single_sign), best


def polynomial(n):
    return (n*n - 1) * (n*n + 1) * (n*n + n + 1)


def supported_part(value, source):
    """Keep every power in value whose prime divides source, regardless of multiplicity."""
    rest = value
    while (common := gcd(rest, source)) > 1:
        rest //= common
    return value // rest


def check_polynomial_witnesses(radius):
    threshold = isqrt(radius)
    assert threshold**2 == radius
    primes = list(map(int, primerange(radius, 2*radius)))
    factors = {d: tuple(map(int, factorint(d))) for d in range(threshold, radius+1)}
    # Count integer pairs by independent congruences on the factor and its partner.
    incidences = {d: 0 for d in factors}
    hits = large = 0
    for p in primes:
        prime_divisors = tuple(map(int, factorint(p*p-1)))
        for q in primes:
            if p == q:
                continue
            source = polynomial(p*q)
            powers = local_powers(p*p-1, prime_divisors)
            supported = [(ell, v) for ell, v in zip(prime_divisors, powers) if source % ell == 0]
            target = prod(v for _, v in supported)
            assert target == supported_part(p*p-1, source)
            if target < threshold:
                continue
            hits += 1
            if any(ell > threshold for ell, _ in supported):
                large += 1
                continue
            d = 1
            for ell, v in supported:
                while v > 1 and d < threshold:
                    v //= ell
                    d *= ell
            assert threshold <= d <= radius and (p*p-1) % d == 0
            assert source % prod(factors[d]) == 0
            incidences[d] += 1
    for d, count in incidences.items():
        radical = prod(factors[d])
        roots_at_p = sum(a*a % d == 1 for a in range(d))
        roots_at_q = prod(sum(polynomial(a) % ell == 0 for a in range(ell))
                          for ell in factors[d])
        assert roots_at_q <= 6**len(factors[d])
        cap = roots_at_p * roots_at_q * ((radius+d-1)//d) * ((radius+radical-1)//radical)
        assert count <= cap
    assert hits == large + sum(incidences.values())
    return hits, large, hits-large


def main():
    pool = tuple(map(int, primerange(2, 312)))  # first 64 primes
    assert len(pool) == 64 and pool[:4] == (2, 3, 5, 7)
    print("bits pairs common single either_sign extra_promise best_four_of_64 polynomial_support")
    for bits in (16, 32, 48, 64):
        count = 512
        pairs = list(prime_pairs(bits, count, seed_offset=8_134_611))
        cases = [coverage(p, q, pool) for p, q in pairs]
        polynomial_hits = sum(max(supported_part(r*r-1, polynomial(p*q))
                                  for r in (p, q))**16 > (p*q)**5 for p, q in pairs)
        print(bits, count, *(sum(column) for column in zip(*cases)), polynomial_hits)
    for radius in (16, 64, 256):
        print("polynomial witnesses", radius, *check_polynomial_witnesses(radius))


if __name__ == "__main__":
    main()
