"""Deterministic mixed-sign, dyadic and exhaustive checks for the public reader."""

from fractions import Fraction
from itertools import combinations
from math import gcd, prod

from sympy import isprime, primerange

if __package__:
    from .smooth_order_cover import _profile_cover, recover_smooth_order, unit_roots
else:
    from smooth_order_cover import _profile_cover, recover_smooth_order, unit_roots


def supported_part(value, primes):
    part = 1
    for prime in primes:
        while value % prime == 0:
            part *= prime
            value //= prime
    return part


def certify_large_fixtures():
    # Complete Lucas certificates with base 3. All auxiliary primes are below
    # 2^64, the range in which SymPy's primality test is deterministic.
    fixtures = (
        (158263163925781249999, {2: 1, 3: 19, 13: 1, 61: 1, 3253: 1, 26393: 1}),
        (237394745888671876507, {2: 1, 3: 1, 17: 1, 31: 1, 13841: 1, 5424275864593: 1}),
    )
    for prime, factors in fixtures:
        assert all(ell < 2**64 and isprime(ell) for ell in factors)
        assert prod(ell**v for ell, v in factors.items()) == prime-1
        assert pow(3, prime-1, prime) == 1
        assert all(gcd(pow(3, (prime-1)//ell, prime)-1, prime) == 1 for ell in factors)


def check_mixed_signs():
    certify_large_fixtures()
    for r, q in ((701, 1087), (9121249, 13681883),
                 (158263163925781249999, 237394745888671876507)):
        n = r*q
        if q < 2**64:
            assert isprime(r) and isprime(q)
        assert r < q < 2*r
        minus, plus = (supported_part(v, (3, 5)) for v in (r-1, r+1))
        assert max(minus, plus)**4 < n and (minus*plus)**8 > n**3
        assert gcd(n*n-1, r*r-1) == gcd(q*q-1, r*r-1) == 24
        assert 24**4 < n and supported_part(q*q-1, (3, 5))**4 < n
        work = {}
        assert recover_smooth_order(n, (3, 5), stats=work) == r
        d = work['winning_modulus']
        assert n <= d**4 and (r*r-1) % d == 0 and work['winning_residue'] == r % d
        assert gcd(d, r-1) > 1 and gcd(d, r+1) > 1
        assert work['tiles'] <= 10*work['residues']
        if n.bit_length() >= 47:
            assert work['winning_method'] == 'lattice' and work['winning_radius'] > 2
        print(f"{n.bit_length()} bits, factor {r}, sides {minus} and {plus}, common 24, "
              f"quantum {work['quantum']}, profiles {work['tried']}/{work['profiles']}, "
              f"winning D {d}, {work['lattices']} lattice calls, "
              f"winning radius {work['winning_radius']}")


def check_roots():
    for v in range(1, 10):
        for factors in ({2: v}, {2: v, 3: 2}):
            d = prod(p**e for p, e in factors.items())
            assert set(unit_roots(factors)) == {x for x in range(d) if (x*x-1) % d == 0}
    work = {}
    assert recover_smooth_order(257*379, (2,), stats=work) == 257
    assert (work['winning_modulus'], work['winning_residue']) == (512, 257)
    assert recover_smooth_order(101*103, (7,)) is None
    assert recover_smooth_order(101*103, ()) is None
    assert recover_smooth_order(101*103, (3, 101)) == 101


def check_small_pairs():
    pairs = [(p, q) for p, q in combinations(map(int, primerange(67, 192)), 2) if q <= 2*p]
    cases = promised = larger_only = 0
    for p, q in pairs:
        n = p*q
        for bases in ((2,), (3, 5), (2, 3, 5)):
            cases += 1
            eligible = [r for r in (p, q) if supported_part(r*r-1, bases)**8 > n**3]
            result = recover_smooth_order(n, bases)
            assert result is None or result in (p, q)
            if eligible:
                promised += 1
                larger_only += eligible == [q]
                assert result in (p, q)
                _, _, profiles = _profile_cover(n, bases, Fraction(1, 8))
                assert len({d for d, _ in profiles}) == len(profiles)
                for r in eligible:
                    assert any((r*r-1) % d == 0 and r % d in unit_roots(factors)
                               for d, factors in profiles)
    assert larger_only > 0 and cases > promised
    print(f"{len(pairs)} balanced pairs, {cases} support cases, {promised} promised, "
          f"{larger_only} promised only at the larger prime")


if __name__ == '__main__':
    check_roots()
    check_mixed_signs()
    check_small_pairs()
