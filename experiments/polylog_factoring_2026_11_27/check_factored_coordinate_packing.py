"""Verify exact coordinate sampling, both supplied factorizations, and a native recovery."""

from fractions import Fraction
from math import prod
from random import Random
from sympy import factorint, isprime, jacobi_symbol

if __package__:
    from .factored_coordinate_packing import factored_coordinate_packing, uniform_factored
    from .check_fixed_field_packing import sphere_progression
    from .fixed_field_packing import recover_factor
else:
    from factored_coordinate_packing import factored_coordinate_packing, uniform_factored
    from check_fixed_field_packing import sphere_progression
    from fixed_field_packing import recover_factor


def check_sampling():
    rng = Random(98271)
    for h in (1, 2, 6, 16, 63):
        seen = set()
        for _ in range(3000):
            value, factors = uniform_factored(h, rng)
            assert 1 <= value <= h and factors == factorint(value)
            seen.add(value)
        assert seen == set(range(1, h+1))

    # Independently sum one trial's exact accepted probabilities, conditioning
    # on the next strictly smaller state to remove the chain's self-loops.
    for h in range(1, 19):
        def recurse(top, product, mass):
            if top == 1:
                accepted[product] += mass * Fraction(product, h)
                return
            for repeats in range(1, h.bit_length()+1) if isprime(top) else (1,):
                value = product * (top**repeats if isprime(top) else 1)
                if value > h:
                    break
                transition = Fraction(1, top**repeats) if isprime(top) else Fraction(1, top-1)
                for smaller in range(1, top):
                    recurse(smaller, value, mass * transition)
        accepted = [Fraction() for _ in range(h+1)]
        for first in range(1, h+1):
            recurse(first, 1, Fraction(1, h))
        assert min(accepted[1:]) > 0 and len(set(accepted[1:])) == 1


def multiply(factors):
    return prod(p**e for p, e in factors.items())


def check_packings():
    for n in (15, 35, 65, 77, 143, 221, 323, 437, 899, 10403, 1022117):
        x = next(i for i in range(40) if jacobi_symbol(-(i*i+1), n) == -1)
        result = factored_coordinate_packing(n, x, Random(3100+n))
        params = result.parameters
        assert isprime(params.T) and n < params.B <= result.height
        assert params.T <= result.height**2 and result.t >= 2
        assert multiply(result.norm_factorization) == params.M
        assert multiply(result.scale_factorization) == params.Q
        assert all(isprime(p) for p in result.norm_factorization | result.scale_factorization)
        print(f"N={n}, x={x}, B={params.B}, t={result.t}, prime T={params.T}, "
              f"radius bits={params.M.bit_length()}, coordinate draws={result.draws}")
        if n == 77:
            assert (params.B, result.t, params.T) == (85, 16, 1525049)
            first_vector = (1, 1, 0)
            modulus = params.Q * 13
            for residue, expected in (
                ((53401, 37081, 16320), [(-14, -63, -42)]),
                ((58841, 31641, 65280), []),
            ):
                assert (sum(c*c for c in residue) - params.M) % modulus == 0
                label = []
                for u, c in zip(first_vector, residue):
                    z, remainder = divmod(params.P*u-c, params.Q)
                    assert remainder == 0
                    label.append(z % 13)
                assert (sum(z*z for z in label)-n*n) % 13 == 0
                assert (sum(u*z for u, z in zip(first_vector, label))+n*x) % 13 == 0
                actual = [tuple((params.P*u-c)//params.Q for u, c in zip(first_vector, w))
                          for w in sphere_progression(params.M, residue, modulus)]
                assert actual == expected
            counts = [0, 0]

            def decide(norm, modulus, residue):
                assert norm == params.M
                ell, rem = divmod(modulus, params.Q)
                assert not rem and isprime(ell) and multiply(result.query_factorization(ell)) == modulus
                answer = next(sphere_progression(norm, residue, modulus), None) is not None
                counts[answer] += 1
                return answer
            factor = recover_factor(n, x, decide, gaussian=(params.B, result.t))
            assert factor == 7 and all(counts)
            print(f"  complete N=77 recovery, {counts[1]} occupied and {counts[0]} empty queries")


if __name__ == "__main__":
    check_sampling()
    check_packings()
