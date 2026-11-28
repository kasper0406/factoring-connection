"""Check carry-syndrome-bounds.md, grade zero-template decoding, and replay public calls."""

from argparse import ArgumentParser
from collections import Counter, defaultdict
from itertools import combinations

from sympy import primerange

if __package__:
    from .carry_parity import binary_product, factor_from_carry_parity
    from .check_root_scale_order_hints import prime_pairs
else:
    from carry_parity import binary_product, factor_from_carry_parity
    from check_root_scale_order_hints import prime_pairs


BITS = (12, 16, 20, 24, 32, 48, 64)
SAMPLE = 1024
SCREEN_OFFSET = 7_000_087
PUBLIC_OFFSET = 9_000_089


def remainder(value, modulus):
    """Remainder of a coefficient bit word, by polynomial long division."""
    width = modulus.bit_length()
    while value.bit_length() >= width:
        value ^= modulus << (value.bit_length()-width)
    return value


def syndrome_radius(n, smaller_prime):
    """Minimum mask weight up to two for which the smaller factor divides the word.

    This offline test is equivalent to the uncapped public decoder on balanced
    distinct semiprimes. Its inputs include the actual smaller factor.
    """
    target = remainder(n, smaller_prime)
    if not target:
        return 0
    prior = set()
    can_pair = False
    monomial = 1
    for _ in range(1, n.bit_length()):
        monomial <<= 1
        if monomial.bit_length() == smaller_prime.bit_length():
            monomial ^= smaller_prime
        if monomial == target:
            return 1
        if target ^ monomial in prior:
            can_pair = True
        prior.add(monomial)
    return 2 if can_pair else None


def check_small():
    """Exhaust independent corrected-polynomial division against the syndrome test."""
    count = 0
    splits = [0]*3
    for p, q in combinations(map(int, primerange(64, 128)), 2):
        n = p*q
        direct = None
        for weight in range(3):
            for subset in combinations(range(1, n.bit_length()), weight):
                mask = sum(1 << j for j in subset)
                if remainder(n ^ mask, p) == 0:
                    direct = weight
                    break
            if direct is not None:
                break
        assert direct == syndrome_radius(n, p)
        if direct is not None:
            for radius in range(direct, 3):
                splits[radius] += 1
        count += 1
    assert splits == [1, 16, 60]
    return count


def check_collision_bound():
    instances = collisions = 0
    for bits in range(4, 9):
        radius = 1 << (bits-1)
        for k in range(1, min(3, bits-2)+1):
            width = 1 << (bits+k)
            middle_bits = bits-k-2
            multipliers = range(1, 1 << (k+1))
            mean_weights = Counter()
            for p in range(radius+1, 2*radius, 2):
                blocks = defaultdict(lambda: defaultdict(list))
                for q in range(radius, 2*radius):
                    blocks[p*q//width][remainder(p*q, p)].append(q)
                # Each nonzero multiplier names exactly one possible XOR word.
                words = {binary_product(p, v): v for v in multipliers}
                gaps = defaultdict(Counter)
                for fibers in blocks.values():
                    for partners in fibers.values():
                        for q, other in combinations(partners, 2):
                            word = (p*q) ^ (p*other)
                            gap = other-q
                            assert word in words and 0 < gap < 1 << (k+1)
                            assert word - p*gap == 2*((p*q) & word)
                            gaps[word][gap] += 1
                actual = sum(sum(counts.values()) for counts in gaps.values())
                collisions += actual
                upper = 0
                for word, v in words.items():
                    low_weight = (word % radius).bit_count()
                    gap_weight = (word % (1 << k)).bit_count()
                    assert len(gaps[word]) <= 1 << gap_weight
                    assert all(count <= radius >> low_weight for count in gaps[word].values())
                    upper += radius >> (low_weight-gap_weight)
                    middle_weight = ((word % radius) >> (k+1)).bit_count()
                    mean_weights[v] += 1 << (middle_bits-middle_weight)
                assert actual <= upper
                for menu_size in (1, 3):
                    selected = sum(sum(sorted(map(len, fibers.values()), reverse=True)[:menu_size])
                                   for fibers in blocks.values())
                    assert selected <= menu_size*len(blocks)+actual
                instances += 1
            # The average of 2^(-middle_weight) is exactly (3/4)^middle_bits.
            assert all(mean_weights[v] * (1 << middle_bits) == (radius//2)*3**middle_bits
                       for v in multipliers)
    return instances, collisions


def screen(bits):
    pairs = list(prime_pairs(bits, SAMPLE, seed_offset=SCREEN_OFFSET))
    predicted = [0]*3
    covered = [0]*3
    positives = []
    example = None
    for p, q in pairs:
        n = p*q
        radius = syndrome_radius(n, p)
        true_weight = (n ^ binary_product(p, q)).bit_count()
        assert true_weight > 2 or (radius is not None and radius <= true_weight)
        for h in range(3):
            predicted[h] += radius is not None and radius <= h
            covered[h] += true_weight <= h
        if radius is not None:
            positives.append((p, q))
            if example is None and true_weight > radius:
                example = (n, p, q, radius, true_weight)
    return pairs, predicted, covered, positives, example


def replay(pairs):
    """The decoder receives N alone, with its default divisor cap and public radius."""
    successes = [0]*3
    skips = [0]*3
    for p, q in pairs:
        n = p*q
        for radius in range(3):
            stats = {}
            factor = factor_from_carry_parity(n, errors=radius, stats=stats)
            assert factor is None or factor == p
            successes[radius] += factor is not None
            skips[radius] += stats['skipped']
            # Compare only after the actual public call returns.
            grade = syndrome_radius(n, p)
            if not stats['skipped']:
                assert (factor is not None) == (grade is not None and grade <= radius)
    return successes, skips


def run(bits, native):
    for size in bits:
        pairs, predicted, covered, positives, example = screen(size)
        print(f'{size} bits, offline {len(pairs)}: uncapped splits {predicted}, '
              f'true words {covered}, first accidental (N,p,q,radius,true weight) {example}',
              flush=True)
        if native:
            actual, skipped = replay(positives)
            print(f'{size} bits, actual calls on all {len(positives)} screen positives: '
                  f'splits {actual}, skipped masks {skipped}', flush=True)
            unfiltered = 4 if size <= 24 else 1
            independent = list(prime_pairs(size, unfiltered, seed_offset=PUBLIC_OFFSET,
                                          exclude=pairs))
            assert set(independent).isdisjoint(pairs)
            actual, skipped = replay(independent)
            print(f'{size} bits, actual calls on {unfiltered} disjoint unfiltered inputs: '
                  f'splits {actual}, skipped masks {skipped}', flush=True)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--bits', type=int, nargs='+', choices=BITS, default=BITS)
    parser.add_argument('--replay', action='store_true', help='also run the full public decoder')
    args = parser.parse_args()
    print(f'Checked all {check_small()} previously measured small inputs by direct mask enumeration.')
    cases, collisions = check_collision_bound()
    print(f'Checked {cases} fixed-factor/block-size cases and {collisions} collision pairs.')
    run(args.bits, args.replay)
