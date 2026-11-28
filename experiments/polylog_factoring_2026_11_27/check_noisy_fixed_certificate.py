"""Exact weighted carry bounds for the two fully observed fixed multipliers."""

from fractions import Fraction
from itertools import product

if __package__:
    from .check_two_fixed_and_streams import advance, carry_matrix
else:
    from check_two_fixed_and_streams import advance, carry_matrix


def all_transitions(h):
    rows = {}
    for state in product(range(h), repeat=4):
        rows[state] = [
            tuple((s + h*x)//2 for s, x in zip(state, bits))
            for bits in product(range(2), repeat=4)
        ]
    return rows


def certificate(h, upper, error_denominator):
    exact = carry_matrix(h)
    total = all_transitions(h)
    v = dict.fromkeys(total, 1)
    # 32 times the tilted matrix gives weight 32 to a match and 1 to a miss.
    for _ in range(12):
        matches = advance(exact, v)
        all_choices = advance(total, v)
        following = {s: 31*matches[s] + all_choices[s] for s in v}
        previous, v = v, following
    ratio = max(Fraction(v[s], 32*previous[s]) for s in v)
    assert ratio < upper
    assert upper**error_denominator * 32**2 < 8**error_denominator
    print(h, ratio, "<", upper, "at error rate", Fraction(1, error_denominator))


if __name__ == "__main__":
    certificate(3, Fraction(37, 5), 90)
    certificate(5, Fraction(71, 10), 60)
