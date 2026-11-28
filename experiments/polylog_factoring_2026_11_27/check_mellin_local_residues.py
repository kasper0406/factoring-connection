"""Compute and check the fixed local residue of the cubic order on semiprimes."""

from math import gcd, isqrt


if __package__:
    from .mellin_local_residues import forced_power_residue, local_h_residue
else:
    from mellin_local_residues import forced_power_residue, local_h_residue


def _phi_three(x):
    return x * x + x + 1


def _check_local_completions():
    checked = 0
    for base, modulus in ((2, 256), (3, 27), (5, 5), (7, 7), (11, 11), (13, 13), (17, 17), (19, 19)):
        units = [x for x in range(1, modulus) if gcd(x, modulus) == 1]
        for product in units:
            # The largest divisor of this modulus on which the local answer
            # is constant is the gcd of all differences from the x=1 answer.
            variation = modulus
            baseline = 3 * _phi_three(product)
            for x in units:
                y = product * pow(x, -1, modulus) % modulus
                variation = gcd(variation, _phi_three(x) * _phi_three(y) - baseline)
            # Lift this local product to an eligible odd residue coprime to 3.
            n = product
            while n <= 1 or gcd(n, 6) != 1:
                n += modulus
            _, supplied = local_h_residue(n)
            assert variation == gcd(supplied, modulus), (base, product, variation)
            checked += 1
    return checked


def _check_prime_pairs():
    primes = [p for p in range(5, 512, 2) if all(p % d for d in range(3, isqrt(p) + 1, 2))]
    checked = 0
    moduli = set()
    for i, p in enumerate(primes):
        for q in primes[i + 1:]:
            if q > 2 * p:
                break
            n = p * q
            residue, modulus = local_h_residue(n)
            # Direct four-divisor convolution defining H, independent of the
            # local congruence and the condensed prime-product formula.
            exact = n * n + q * q * (p + 1) + p * p * (q + 1) + (p + 1) * (q + 1)
            assert exact % modulus == residue, (n, exact, residue, modulus)
            moduli.add(modulus)
            checked += 1
    assert moduli == {8, 32, 64, 128, 72, 288, 576, 1152}
    return checked


def _prime_power_value(x, j):
    # Sum (x^j / d)^2 psi(d) over the formal divisors d=1,x,...,x^j.
    return x**(2*j) + (x+1)*sum(x**(2*j-k-1) for k in range(1, j+1))


def _check_powers():
    checked = 0
    for p, q in ((41, 79), (53, 61), (1009, 1231)):
        n = p*q
        assert forced_power_residue(n, 1) == local_h_residue(n)
        for j in (1, 2, 3, 4, 6, 10, 16, 31, 64):
            residue, modulus = forced_power_residue(n, j)
            assert _prime_power_value(p, j) * _prime_power_value(q, j) % modulus == residue
            for test in (256, 27, 125, 49, 121, 13, 17, 19, 23, 31):
                common = test
                baseline = _prime_power_value(1, j) * _prime_power_value(n % test, j)
                for x in range(1, test):
                    if gcd(x, test) == 1:
                        y = n * pow(x, -1, test) % test
                        common = gcd(common, _prime_power_value(x, j) * _prime_power_value(y, j) - baseline)
                assert common == gcd(test, modulus), (p, q, j, test)
                checked += 1
    return checked


if __name__ == '__main__':
    print(f'{_check_prime_pairs()} balanced prime pairs; '
          f'{_check_local_completions()} exhaustive local product classes')
    print(f"{_check_powers()} exhaustive powered-modulus comparisons")
