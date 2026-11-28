"""Check the carry-matrix certificate and the two-stream fiber argument."""

from collections import Counter, defaultdict
from itertools import product


def carry_matrix(h):
    """Integer transition matrix, represented by its outgoing edges with multiplicity."""
    rows = {}
    for state in product(range(h), repeat=4):
        successors = []
        for bits in product(range(2), repeat=4):
            if bits[0] & bits[1] != bits[2] & bits[3]:
                continue
            values = tuple(carry + h * bit for carry, bit in zip(state, bits))
            if values[0] & values[1] & 1 != values[2] & values[3] & 1:
                continue
            successors.append(tuple(value // 2 for value in values))
        rows[state] = successors
    return rows


def advance(rows, vector):
    return {state: sum(vector[successor] for successor in row) for state, row in rows.items()}


def matrix_certificate(h):
    rows = carry_matrix(h)
    vector = dict.fromkeys(rows, 1)
    maxima = [1]
    for _ in range(7):
        vector = advance(rows, vector)
        maxima.append(max(vector.values()))
    next_vector = advance(rows, vector)
    gap = max(5 * next_vector[state] - 36 * vector[state] for state in rows)
    if h == 3:
        assert maxima[3] == 508 < 8**3
        assert (min(vector.values()), max(vector.values()), gap) == (584368, 1319652, -16424)
    else:
        assert h == 5
        assert maxima[5] == 28744 < 8**5
        assert (min(vector.values()), max(vector.values()), gap) == (410192, 1390424, -724264)
    return maxima


def observation(x, y, h, modulus, skip=0):
    mask = (modulus - 1) ^ ((1 << skip) - 1)
    return x * y % modulus, (x & y) & mask, (h * x & h * y) & mask


def direct_first_differences(h, d, b):
    """Count ordered pairs of compatible leaves whose first coordinates first differ at d."""
    modulus = 1 << (d + b)
    groups = defaultdict(lambda: [0, 0])
    for x in range(1, modulus, 2):
        for y in range(1, modulus, 2):
            key = observation(x, y, h, modulus), x % (1 << d)
            groups[key][(x >> d) & 1] += 1
    return sum(2 * left * right for left, right in groups.values())


def block_first_differences(h, d, b):
    """Count the same pairs via low slopes and explicit high words."""
    modulus = 1 << b
    low_modulus = 1 << d
    low_counts = defaultdict(Counter)
    for x in range(1, low_modulus, 2):
        for y in range(1, low_modulus, 2):
            slope = -y * pow(x, -1, modulus) % modulus
            low_counts[slope][h * x // low_modulus, h * y // low_modulus] += 1
    count = 0
    for r, s, other_r, other_s in product(range(modulus), repeat=4):
        if not ((r ^ other_r) & (s ^ other_s) & 1) or r & s != other_r & other_s:
            continue
        slope = (other_s - s) * pow(other_r - r, -1, modulus) % modulus
        for (carry_x, carry_y), multiplicity in low_counts[slope].items():
            first = ((h * r + carry_x) & (h * s + carry_y)) % modulus
            other = ((h * other_r + carry_x) & (h * other_s + carry_y)) % modulus
            if first == other:
                count += multiplicity
    return count


def check_first_differences(h, maxima):
    total = 0
    for d in range(1, 5):
        for b in range(1, d + 1):
            count = direct_first_differences(h, d, b)
            assert count == block_first_differences(h, d, b), (h, d, b)
            pairs = 1 << (2 * (d + b - 1))
            assert count * (8**b) <= pairs * 2 * h**2 * maxima[b]
            total += count
    return total


def check_good_leaves(h):
    checked = 0
    census = {}
    for k in range(4, 10):
        modulus = 1 << k
        for skip in (0, 2):
            fibers = defaultdict(list)
            for x in range(1, modulus, 2):
                for y in range(1, modulus, 2):
                    fibers[observation(x, y, h, modulus, skip)].append(x)
            if k == 9:
                census[skip] = max(map(len, fibers.values()))
            for t in range(max(1, skip), k // 2 + 1):
                bad_pairs = large_fiber_pairs = 0
                for leaves in fibers.values():
                    # If a low-t class contains two longer prefixes, each leaf in
                    # it has another leaf whose first difference is in [t,k-t].
                    low_classes = defaultdict(set)
                    for x in leaves:
                        low_classes[x % (1 << t)].add(x % (1 << (k - t + 1)))
                    good = sum(len(low_classes[x % (1 << t)]) == 1 for x in leaves)
                    assert good <= 4**t
                    bad_pairs += len(leaves) - good
                    if len(leaves) > 2 * 4**t:
                        large_fiber_pairs += len(leaves)
                    checked += 1
                assert large_fiber_pairs <= 2 * bad_pairs
    return checked, census


def check_planted_menu():
    tested = 0
    for s in range(2, 9):
        half_word = 1 << s
        modulus = half_word**2
        for a in range(1, half_word, 2):
            x, y = 1 + half_word * a, (1 - half_word * a) % modulus
            assert x * y % modulus == 1
            for h in range(1, half_word, 2):
                assert (h * x & h * y) % modulus == h + half_word
                tested += 1
    return tested


if __name__ == "__main__":
    for multiplier in (3, 5):
        count_maxima = matrix_certificate(multiplier)
        differences = check_first_differences(multiplier, count_maxima)
        comparisons, maxima_at_nine = check_good_leaves(multiplier)
        print(f"h={multiplier}: exact carry certificate, {differences} first-difference pairs, "
              f"{comparisons} fiber/scale checks.")
        print(f"At 9 bits, maximum ordered fiber was {maxima_at_nine[0]} "
              f"with complete streams and {maxima_at_nine[2]} omitting their two lowest bits.")
    print(f"Verified {check_planted_menu()} planted worst-case multiplier/parameter combinations.")
