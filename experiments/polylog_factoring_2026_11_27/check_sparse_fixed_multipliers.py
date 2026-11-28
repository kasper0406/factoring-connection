"""Exact certificates for periodically selected bits of fixed-multiplier AND streams."""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import log2

if __package__:
    from .check_two_fixed_and_streams import advance, carry_matrix
else:
    from check_two_fixed_and_streams import advance, carry_matrix


def unscaled_matrix(h):
    """Count four input digits when only the unscaled AND outputs must match."""
    rows = {}
    for state in product(range(h), repeat=4):
        rows[state] = [
            tuple((carry+h*bit)//2 for carry, bit in zip(state, bits))
            for bits in product(range(2), repeat=4)
            if (bits[0] & bits[1]) == (bits[2] & bits[3])
        ]
    return rows


def check_period(h, pattern, expected):
    matrices = (unscaled_matrix(h), carry_matrix(h))
    def apply_period(vector):
        # Digits arrive low to high; the leftmost matrix acts on the vector last.
        for selected in reversed(pattern):
            vector = advance(matrices[selected], vector)
        return vector
    vector = apply_period(dict.fromkeys(matrices[0], 1))
    two_periods = apply_period(vector)
    ratio = max(Fraction(two_periods[s], vector[s]) for s in vector)
    assert min(vector.values()) > 0 and ratio == expected < 8**len(pattern)

    # Independent complete one-period comparison: group actual pairs of words
    # by the two masked observations and multiply the bucket multiplicities.
    modulus = 1 << len(pattern)
    mask = sum(bit << j for j, bit in enumerate(pattern))
    counts = {}
    for cx, cy in product(range(h), repeat=2):
        counts[cx, cy] = Counter(
            (x & y, (h*x+cx) & (h*y+cy) & mask)
            for x, y in product(range(modulus), repeat=2)
        )
    for (c1, c2, c3, c4), value in vector.items():
        others = counts[c3, c4]
        assert value == sum(n*others[key] for key, n in counts[c1, c2].items())
    return len(vector), ratio, float(ratio)**(1/len(pattern))


def binary_entropy(x):
    return -x*log2(x)-(1-x)*log2(1-x)


def conditional_entropy(h):
    # Independent stationary carries and next input bits have uniform joint law.
    counts = Counter()
    for c, d in product(range(h), repeat=2):
        for x, y in product(range(2), repeat=2):
            counts[x & y, (c+h*x) & (d+h*y) & 1] += 1
    total = 4*h*h
    joint_one = Fraction(counts[1, 1], total)
    even_carry = Fraction(h+1, 2*h)
    assert joint_one == even_carry**2/4
    assert counts[1, 0]+counts[1, 1] == counts[0, 1]+counts[1, 1] == total//4
    z = float(even_carry**2)
    conditional = (binary_entropy(z)+3*binary_entropy((1-z)/3))/4
    return joint_one, conditional, (1-binary_entropy(.25))/conditional


def context_entropy(h, length):
    """Exact counts for a scaled AND bit given a length-bit unscaled AND suffix."""
    states = Counter((c, d, 0) for c, d in product(range(h), repeat=2))
    outputs = Counter()
    for j in range(length):
        after = Counter()
        for (c, d, history), count in states.items():
            for x, y in product(range(2), repeat=2):
                pattern = history + ((x & y) << j)
                left, right = c+h*x, d+h*y
                if j+1 == length:
                    outputs[pattern, left & right & 1] += count
                else:
                    after[left//2, right//2, pattern] += count
        states = after
    total = h*h*(4**length)
    assert sum(outputs.values()) == total
    entropy = 0
    for pattern in range(1 << length):
        no, yes = outputs[pattern, 0], outputs[pattern, 1]
        assert no+yes == h*h*3**(length-pattern.bit_count())
        if no and yes:
            entropy += (no+yes)/total*binary_entropy(yes/(no+yes))
    return entropy, (1-binary_entropy(.25))/entropy


if __name__ == "__main__":
    for h, pattern, expected in ((3, (0, 1, 1), Fraction(64591, 127)),
                                 (5, (0, 1, 0, 1, 1), Fraction(286956018, 9253))):
        print("h, pattern, carry rows, max ratio, rate:", h, pattern, *check_period(h, pattern, expected))
        print("Joint one frequency, conditional entropy, necessary density:", *conditional_entropy(h))
        print("Eight-bit context entropy and necessary density:", *context_entropy(h, 8))
