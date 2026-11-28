"""Compare both native parity paths with independent scalar sphere enumeration."""

from collections import Counter
from math import isqrt
from random import Random

from sympy import jacobi_symbol

if __package__:
    from .native_ternary_parity import native_ternary_parity, packed_parity_shortcut
    from .packed_ternary_decision import brute_sphere_count, brute_sphere_parity, recover_packed_factor
else:
    from native_ternary_parity import native_ternary_parity, packed_parity_shortcut
    from packed_ternary_decision import brute_sphere_count, brute_sphere_parity, recover_packed_factor


def scalar_pair_tests(norm, modulus, residue):
    """Count the integer square tests in the two-coordinate scalar reference."""
    rx, ry, _ = residue
    limit = isqrt(norm)
    tests = 0
    for w in range(-limit + (rx+limit) % modulus, limit+1, modulus):
        ymax = isqrt(norm-w*w)
        first = -ymax + (ry+ymax) % modulus
        tests += max(0, (ymax-first)//modulus + 1)
    return tests


def compare(n):
    x = next(t for t in range(n) if jacobi_symbol(-(t*t+1), n) == -1)
    tally = Counter()

    def grade(norm, modulus, residue):
        stats = {}
        parity = native_ternary_parity(norm, modulus, residue, stats=stats)
        count = brute_sphere_count(norm, modulus, residue)
        shortcut = packed_parity_shortcut(norm, modulus, residue)
        method = stats["method"]
        assert method in ("mirror", "conic")
        assert parity == bool(count % 2)
        assert shortcut == (parity if method == "mirror" else None)
        tally[method] += 1
        tally[method + "_odd"] += parity
        tally["conic_tests"] += stats["coordinate_tests"]
        tally["scalar_tests"] += scalar_pair_tests(norm, modulus, residue)
        return parity

    result = recover_packed_factor(n, x, grade)
    assert 1 < result.factor < n and n % result.factor == 0
    return (n, x, result.factor, tally["mirror"], tally["mirror_odd"], tally["conic"],
            tally["conic_odd"], tally["conic_tests"], tally["scalar_tests"])


def check_fallback():
    rng = Random(603319)
    for _ in range(256):
        n, modulus = rng.randrange(100), rng.randrange(1, 15)
        c = tuple(rng.randrange(-modulus, 2*modulus) for _ in range(3))
        assert native_ternary_parity(n, modulus, c) == brute_sphere_parity(n, modulus, c)
    # A cube reflection cancels two distinct non-coordinate partners. It can
    # return false on a nonempty class, just as the scalar parity must.
    n, x, scale = 15, 1, 7**4
    norm = (x+scale*n)**2 + 1
    u = (0, 1, 1)
    assert brute_sphere_count(norm, scale, u) == 6
    assert packed_parity_shortcut(norm, scale, u) is False


if __name__ == "__main__":
    check_fallback()
    print("N x divisor mirror mirror-odd generic generic-odd native-tests scalar-tests")
    for n in (15, 35, 65, 77, 91, 183):
        print(*compare(n))
