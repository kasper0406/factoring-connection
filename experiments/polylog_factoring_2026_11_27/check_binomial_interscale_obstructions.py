"""Verify shift quotients, their orbit markers, and literal modular examples."""

from collections import Counter
from fractions import Fraction
from math import comb, gcd, prod


def normalized_binomial(t):
    return Fraction(comb(2 * t, t), 4**t)


def check_shift_ratios():
    profiles = {}
    comparisons = 0
    for m in range(1, 41):
        numerator = [Fraction(2 * j + 1, 2 * m) for j in range(m)]
        denominator = [Fraction(j + 1, m) for j in range(m)]
        profile = Counter(value % 1 for value in numerator)
        profile.subtract(value % 1 for value in denominator)
        marker = Fraction(1, 2 * m)
        assert profile[marker] == 1
        assert all(earlier[marker] == 0 for earlier in profiles.values())
        profiles[m] = profile
        for n in range(9):
            ratio = prod(n + a for a in numerator) / prod(n + a for a in denominator)
            assert ratio == normalized_binomial(m * (n + 1)) / normalized_binomial(m * n)
            comparisons += 1
    return comparisons


def check_modular_families():
    fixtures = [(16, 31, 47), (16, 71, 103), (32, 127, 191), (32, 47, 79)]
    comparisons = 0
    for finest, p, q in fixtures:
        n = p * q
        assert finest // 2 < p < q < 2 * p
        assert n % finest == 1 and p % finest == q % finest
        high = p % finest == finest - 1
        assert high or p % finest == finest // 2 - 1
        d = 4
        while d <= finest:
            for h in range(1, d // 4 + 1, 2):
                t = h * n // d
                residue = comb(2 * t, t) % n
                if d < finest or high:
                    assert residue == 0, (n, finest, h, d)
                else:
                    assert gcd(n, residue) == 1, (n, finest, h, d)
                comparisons += 1
            d *= 2
    return comparisons


if __name__ == "__main__":
    print(f"Verified {check_shift_ratios()} shift ratios and 40 orbit markers.")
    print(f"Checked {check_modular_families()} literal binomials on the modular families.")
