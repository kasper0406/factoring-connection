"""Use shifted count bits modulo powers of two to locate component residues."""

import sys
from collections import Counter
from itertools import combinations
from math import prod
from pathlib import Path

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.polylog_factoring_2026_11_27.large_mod_residue_evaluator import (
    find_component_in_residue,
)


def valuation(value):
    """Return the 2-adic valuation of a nonzero integer."""
    return (value & -value).bit_length() - 1


def two_component_zero_bit(modulus, shift, exponent):
    """Evaluate a zero bit natively, assuming two odd prime-power components."""
    if modulus <= 1 or not modulus % 2 or shift <= 0 or exponent < 1:
        raise ValueError('Expected an odd modulus, positive shift and positive exponent')
    work = Counter()
    if not shift % 2:
        return False, work
    difference = modulus - shift**2
    depth = (exponent + 1) // 2
    if difference:
        depth = max(depth, exponent - valuation(difference))
    work['progression_depth'] = depth
    result = find_component_in_residue(modulus, 1 << depth, -shift, work)
    return result is not None, work


def two_component_branch_bit(modulus, residue, depth, zero_bit):
    """Use one shifted zero bit to test membership in a residue class."""
    shift = -residue % (1 << depth)
    difference = modulus - shift**2
    overlap = min(depth, valuation(difference)) if difference else depth
    return zero_bit(shift, depth + overlap)


def branch_bit(residue, depth, max_components, zero_bit):
    """Test an odd residue class by looking for variation among all short lifts."""
    lift_bits = max_components.bit_length()
    step = 1 << depth
    span = 1 << lift_bits
    first = -residue % step
    # The empty side of any mixed pair has valuation strictly below this cap.
    cap = max_components * (depth + lift_bits)
    low, high = 0, cap
    while low < high:
        middle = (low + high + 1) // 2
        if zero_bit(first, middle):
            low = middle
        else:
            high = middle - 1
    reference = low
    for lift in range(1, span):
        shift = first + lift * step
        if reference and not zero_bit(shift, reference):
            return True
        if reference < cap and zero_bit(shift, reference + 1):
            return True
    return False


def split_with_zero_bits(modulus, zero_bit, *, two_components=False):
    """Recover one component using shifted zero bits followed by a native lattice."""
    if modulus <= 1 or not modulus % 2:
        raise ValueError('Expected an odd modulus greater than one')
    bits = modulus.bit_length()
    stop = (bits + 3) // 4
    residue = 1
    for depth in range(2, stop + 1):
        if two_components:
            occupied = two_component_branch_bit(modulus, residue, depth, zero_bit)
        else:
            occupied = branch_bit(residue, depth, bits, zero_bit)
        if not occupied:
            residue += 1 << (depth - 1)
    work = Counter()
    component = find_component_in_residue(modulus, 1 << stop, residue, work)
    if component is None:
        raise ValueError('Zero-bit replies did not locate a complete component')
    return component, work


def check():
    cases = 0
    fixtures = list(combinations((3, 5, 7, 11, 13, 17, 25, 27), 2))
    fixtures = [pair for pair in fixtures if set(pair) not in ({5, 25}, {3, 27})]
    fixtures += [(9, 25), (101, 1009), (101**2, 1009**3), (1000003, 1000033)]
    for components in fixtures:
        modulus = prod(components)
        shifts = range(1, 25) if modulus < 10**4 else (1, 3, 5, 11, 25)
        for shift in shifts:
            left, right = (valuation(q + shift) for q in components)
            for exponent in range(1, max(4, left + right + 2)):
                actual, _ = two_component_zero_bit(modulus, shift, exponent)
                assert actual == (left + right >= exponent)
                cases += 1

    trie_cases = 0
    fixtures += [(3, 5, 7, 11, 13), (9, 25, 49, 11), (27,),
                 (9, 25, 49, 121, 169, 289, 361, 529)]
    for components in fixtures:
        modulus = prod(components)
        calls = Counter()

        def zero_bit(shift, exponent):
            assert shift > 0 and exponent > 0
            calls['queries'] += 1
            return sum(valuation(q + shift) for q in components) >= exponent

        for depth in range(1, 7):
            for residue in range(1, 1 << depth, 2):
                expected = any(q % (1 << depth) == residue for q in components)
                assert branch_bit(residue, depth, len(components), zero_bit) == expected
                if len(components) == 2:
                    assert two_component_branch_bit(modulus, residue, depth, zero_bit) == expected
                trie_cases += 1
        modes = (False, True) if len(components) == 2 else (False,)
        for two_components in modes:
            calls.clear()
            component, _ = split_with_zero_bits(modulus, zero_bit, two_components=two_components)
            assert component in components
            if two_components:
                assert calls['queries'] == (modulus.bit_length() + 3) // 4 - 1
    print(f'{cases} native thresholds and {trie_cases} residue branches matched')


if __name__ == '__main__':
    check()
