"""Exact finite checks for causal and lookahead predictions of binary carry parity."""

from collections import Counter
from fractions import Fraction
from math import comb


def carry_parities(p, q, bits):
    previous = 0
    output = []
    for i in range(bits-1):
        convolution = sum(((p >> j) & 1)*((q >> (i-j)) & 1) for j in range(i+1))
        total = convolution + previous
        assert total % 2 == ((p*q >> i) & 1)
        previous = total//2
        output.append(previous % 2)
    return output


def oracle_causal(p, q, bits):
    """Predict each outgoing parity with the past factor bits and current N bit."""
    incoming = 0
    predictions = []
    for i in range(1, bits-1):
        inside = sum(((p >> j) & 1)*((q >> (i-j)) & 1) for j in range(1, i))
        offset = incoming + inside
        # If the end bits differ their sum is one. When equal, guess 00.
        forced_parity = ((p*q >> i) - offset) % 2
        predictions.append(((offset + forced_parity)//2) % 2)
        incoming = (offset + ((p >> i) & 1) + ((q >> i) & 1))//2
    return predictions


def best_low_word_counts(word, positions):
    """Optimum leaf count for each Hamming budget, with u and the true past supplied."""
    def visit(i, prefix):
        if i > positions:
            return [1] * (positions+1)
        children = []
        for digit in (0, 1):
            p = prefix + (digit << i)
            q = word * pow(p, -1, 1 << (i+1)) % (1 << (i+1))
            # True parity equals the causal guess XOR the indicator of 11.
            label = digit * ((q >> i) & 1)
            children.append((label, visit(i+1, p)))
        counts = []
        for budget in range(positions+1):
            candidates = []
            for prediction in (0, 1):
                leaves = 0
                for label, child in children:
                    remaining = budget - int(label != prediction)
                    if remaining >= 0:
                        leaves += child[remaining]
                candidates.append(leaves)
            counts.append(max(candidates))
        return counts
    return visit(1, 1)


def binomial_half_cdf(trials, maximum):
    return Fraction(sum(comb(trials, i) for i in range(min(trials, maximum)+1)), 2**trials)


def check_predictions(bits):
    m = bits-2
    lo, hi = 1 << (bits-1), 1 << bits
    histogram = Counter()
    equal_histogram = Counter()
    oracle_errors = [0]*m
    lookahead_errors = [0]*m
    jump_histograms = [Counter(), Counter()]
    for p in range(lo+1, hi, 2):
        for q in range(lo+1, hi, 2):
            carries = carry_parities(p, q, bits)
            parities = carries[1:]
            predicted = oracle_causal(p, q, bits)
            errors = [int(a != b) for a, b in zip(predicted, parities, strict=True)]
            for j, error in enumerate(errors, 1):
                assert error == (((p & q) >> j) & 1)
                oracle_errors[j-1] += error
                # On an equal-bit occasion the predictor can flip its guess
                # according to any future low-word bit, keeping its mean error.
                equal_bits = ((p^q) >> j) & 1 == 0
                extra = (p*q >> ((j+1) % (m+1))) & 1
                lookahead_errors[j-1] += error ^ (extra if equal_bits else 0)
            histogram[sum(errors)] += 1
            equal = sum(((p^q) >> j) & 1 == 0 for j in range(1, bits-1))
            equal_histogram[equal] += 1
            for order in (1, 2):
                # Modulo two, (1-X)^order = 1+X^order for these two orders.
                odd_jumps = sum(carries[j] ^ (carries[j-order] if j >= order else 0)
                                for j in range(1, bits-1))
                jump_histograms[order-1][odd_jumps] += 1
    assert histogram == {i: comb(m, i)*3**(m-i) for i in range(m+1)}
    assert equal_histogram == {i: comb(m, i)*2**m for i in range(m+1)}
    assert oracle_errors == [4**(m-1)]*m
    assert lookahead_errors == oracle_errors
    for jumps in jump_histograms:
        for maximum in range(m+1):
            assert sum(jumps[j] for j in range(maximum+1)) <= sum(
                histogram[j] for j in range(maximum+1))

    best = [0]*(m+1)
    for u in range(1, 2**(m+1), 2):
        counts = best_low_word_counts(u, m)
        for budget, count in enumerate(counts):
            best[budget] += count
    for budget, count in enumerate(best):
        for k in range(1, m+1):
            bound = binomial_half_cdf(m, k-1) + binomial_half_cdf(k, budget)
            assert Fraction(count, 4**m) <= bound
        assert sum(histogram[j] for j in range(budget+1)) <= count <= 4**m
    return histogram[0], best[0], 4**m


def check_high_boxes():
    rows = 0
    max_ratio = Fraction(0)
    # These are finite occupancy measurements of the curved strip, not
    # a choice of the implicit constant in the asymptotic estimate.
    for radius in (16, 64, 256, 1024):
        for modulus in (2, 4, 8):
            if modulus**4 > radius:
                continue
            for width in (radius, 2*radius):
                occupancy = Counter((p*q//width, p % modulus, q % modulus)
                                    for p in range(radius, 2*radius)
                                    for q in range(radius, 2*radius))
                assert sum(occupancy.values()) == radius**2
                max_ratio = max(max_ratio, Fraction(max(occupancy.values())*modulus**2, radius))
                rows += 1
    return rows, max_ratio


if __name__ == '__main__':
    print('factor bits  causal perfect  low-word optimum  ordered odd pairs')
    for b in range(3, 10):
        causal, lookahead, total = check_predictions(b)
        print(f'{b:11} {causal:15} {lookahead:17} {total:18}')
    rows, ratio = check_high_boxes()
    print(f'{rows} high-word occupancy tables. Maximum cell size / (R/M²) = {ratio}.')
