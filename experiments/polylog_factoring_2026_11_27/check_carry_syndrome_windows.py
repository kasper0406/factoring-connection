"""Exhaust the grouping, bit constraints and moment step for interior-window syndromes."""

from collections import Counter, defaultdict
from fractions import Fraction

if __package__:
    from .carry_parity import _divide, binary_product
else:
    from carry_parity import _divide, binary_product


def check(bits, k):
    radius = 1 << (bits-1)
    width = bits+k
    word_mask = (1 << width)-1
    patterns = 1 << (k+1)
    inputs = range(radius+1,2*radius,2)
    partners = range(radius,2*radius)
    windows = 2*bits-width+1
    energy = [0]*windows
    marked_counts = [0]*windows
    cardinality = witnesses = 0
    worst_single = [0]*windows
    probabilities = Counter()
    for p in inputs:
        marked = set()
        for v in range(1,patterns):
            multiple = binary_product(p,v)
            mask = multiple & (radius-1)
            low_ones = (multiple & ((1 << k)-1)).bit_count()
            weight = mask.bit_count()
            mass = Fraction(1, 1 << (weight-low_ones))
            probabilities[v] += mass
            per_multiple = set()
            possible_h = set()
            for h in range(1,patterns):
                delta = multiple-p*h
                if delta < 0 or delta % 2:
                    continue
                fixed = delta//2
                if fixed & ~multiple:
                    continue
                possible_h.add(h)
                per_multiple.update(z for z in range(radius) if z & mask == fixed & mask)
            assert len(possible_h) <= 1 << low_ones
            assert len(per_multiple) <= radius*mass
            marked.update(per_multiple)
        cardinality += len(marked)
        for shift in range(windows):
            buckets = defaultdict(lambda: defaultdict(list))
            frequencies = Counter()
            for q in partners:
                n = p*q
                outside = n & ~(word_mask << shift)
                residue = _divide(n,p)[1]
                buckets[outside][residue].append(q)
                frequencies[(n >> shift) & (radius-1)] += 1
            energy[shift] += sum(t*t for t in frequencies.values())
            eligible_q = sum(frequencies[z] for z in marked)
            marked_counts[shift] += eligible_q
            best = 0
            for group in buckets.values():
                best += max(map(len, group.values()))
                for qs in group.values():
                    for q, other in zip(qs,qs[1:]):
                        diff = p*q ^ (p*other)
                        assert diff % (1 << shift) == 0
                        multiple = diff >> shift
                        v, rem = _divide(multiple,p)
                        h, low = divmod(other-q,1 << shift)
                        assert rem == low == 0 and 0 < v < patterns and 0 < h < patterns
                        assert 2*((p*q >> shift) & multiple) == multiple-p*h
                        assert ((p*q >> shift) & (radius-1)) in marked
                        witnesses += 1
            assert best <= len(buckets)+eligible_q
            worst_single[shift] += best

    for mass in probabilities.values():
        assert mass <= len(inputs)*Fraction(3,4)**(bits-2-k)
    assert all(count**2 <= cardinality*moment for count,moment in zip(marked_counts,energy))
    return windows, witnesses, min(worst_single), max(worst_single)


if __name__ == '__main__':
    for bits, k in ((5,1),(6,1),(7,2),(8,2)):
        print("bits, k, positions, checked collisions, min/max best-table hits:",
              bits, k, *check(bits,k))
