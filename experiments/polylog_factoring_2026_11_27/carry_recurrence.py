"""Recover a balanced divisor by reconstructing two linearly recurrent carry segments."""

from fractions import Fraction
from math import isqrt, lcm

if __package__:
    from .carry_polynomial import _read_factor
else:
    from carry_polynomial import _read_factor


def _fit(seed):
    """Berlekamp–Massey over Q. Return integral recurrence coefficients and divisor."""
    relation = [Fraction(1)]
    backup = [Fraction(1)]
    order, since_update, last_error = 0, 1, Fraction(1)
    for j, value in enumerate(seed):
        error = value + sum(relation[i]*seed[j-i] for i in range(1, order+1))
        if error == 0:
            since_update += 1
            continue
        previous = relation.copy()
        target_length = len(backup) + since_update
        relation.extend([Fraction(0)]*max(0, target_length-len(relation)))
        for i, coefficient in enumerate(backup):
            relation[i+since_update] -= (error/last_error)*coefficient
        if 2*order <= j:
            order = j+1-order
            backup, last_error, since_update = previous, error, 1
        else:
            since_update += 1
    denominator = lcm(*(c.denominator for c in relation))
    return tuple(int(-c*denominator) for c in relation[1:order+1]), denominator


def _extend(seed, order, length):
    """Keep all available seed terms and extend for as long as admissible."""
    values = list(seed)
    weights, denominator = _fit(seed)
    if len(weights) <= order:
        while len(values) < length:
            numerator = sum(c*values[-i] for i, c in enumerate(weights, 1))
            value, remainder = divmod(numerator, denominator)
            if remainder or value < 0 or value > length:
                break
            values.append(value)
    return tuple(values)


def _raw_column(a, b, j):
    """Convolution coefficient at j for two low-aligned bit strings."""
    return sum(((a >> i) & 1) * ((b >> (j-i)) & 1) for i in range(j+1))


def _low_seeds(n, count):
    modulus = 1 << count
    seeds = set()
    for a in range(1, modulus, 2):
        b = n * pow(a, -1, modulus) % modulus
        if a > b:  # The other ordering generates the same low carries.
            continue
        carry = 0
        row = []
        for j in range(count):
            carry = (carry + _raw_column(a, b, j))//2
            row.append(carry)
        seeds.add(tuple(row))
    return seeds


def _reverse_bits(value, count):
    return int(f'{value:0{count}b}'[::-1], 2)


def _high_seeds(n, count):
    """Enumerate possible high bits and reconstruct carries downwards from zero."""
    size = n.bit_length()
    smallest = isqrt((n-1)//2)+1
    largest = isqrt(n)
    seeds = set()
    for a_bits in range(smallest.bit_length(), largest.bit_length()+1):
        shift = max(0, a_bits-count)
        lo = max(smallest, 1 << (a_bits-1))
        hi = min(largest, (1 << a_bits)-1)
        for prefix in range(lo >> shift, (hi >> shift)+1):
            a_lo = max(lo, prefix << shift)
            a_hi = min(hi, ((prefix+1) << shift)-1)
            b_lo, b_hi = (n+a_hi-1)//a_hi, n//a_lo
            a_reversed = _reverse_bits(prefix, a_bits-shift)
            for b_bits in range(b_lo.bit_length(), b_hi.bit_length()+1):
                b_shift = max(0, b_bits-count)
                left = max(b_lo, 1 << (b_bits-1)) >> b_shift
                right = min(b_hi, (1 << b_bits)-1) >> b_shift
                for b_prefix in range(left, right+1):
                    b_reversed = _reverse_bits(b_prefix, b_bits-b_shift)
                    degree = a_bits+b_bits-2
                    if degree not in (size-2, size-1):
                        continue
                    carry = 0
                    row = []
                    for j in range(size-1, size-count-1, -1):
                        offset = degree-j
                        raw = _raw_column(a_reversed, b_reversed, offset) if offset >= 0 else 0
                        carry = ((n >> j) & 1) + 2*carry-raw
                        if not 0 <= carry <= size:
                            break
                        row.append(carry)
                    if len(row) == count:
                        seeds.add(tuple(row))
    return seeds


def _coefficients(bits, internal):
    size = len(bits)
    coefficients = []
    incoming = 0
    for j, outgoing in enumerate((*internal, 0)):
        raw = bits[j]+2*outgoing-incoming
        if not 0 <= raw <= min(j+1, size-j):
            return None
        coefficients.append(raw)
        incoming = outgoing
    while coefficients[-1] == 0:
        coefficients.pop()
    return coefficients if coefficients[-1] == 1 else None


def recurrent_carry_factor(n, order, *, stats=None):
    """Try every cut of the internal carries. A verified proper divisor or None.

    The guarantee concerns odd N=p*q with p<=q<=2p. Each segment of the true
    carry word must have order at most 'order' over Q in its outside-in direction.
    """
    if not isinstance(n, int) or n < 4 or not isinstance(order, int) or order < 1:
        raise ValueError("expected N >= 4 and positive integral recurrence order")
    if n % 2 == 0:
        return 2
    length = n.bit_length()-1
    count = min(2*order, length)
    if stats is None:
        stats = {}
    forward = [_extend(row, order, length) for row in sorted(_low_seeds(n, count))]
    reverse = [_extend(row, order, length) for row in sorted(_high_seeds(n, count))]
    stats.update(low_seeds=len(forward), high_seeds=len(reverse), candidates=0,
                 factorizations=0, winning_cut=None, faithful_lift=None)
    digits = [(n >> j) & 1 for j in range(length+1)]
    seen = set()
    # A centered cut is often sufficient; the guarantee covers every cut.
    for split in sorted(range(length+1), key=lambda h: abs(2*h-length)):
        left = {row[:split] for row in forward if len(row) >= split}
        right = {row[:length-split][::-1] for row in reverse if len(row) >= length-split}
        for prefix in sorted(left):
            for suffix in sorted(right):
                internal = prefix + suffix
                if internal in seen:
                    continue
                seen.add(internal)
                stats['candidates'] += 1
                coefficients = _coefficients(digits, internal)
                if coefficients is None:
                    continue
                stats['factorizations'] += 1
                factor = _read_factor(n, coefficients)
                if factor is not None:
                    stats['winning_cut'] = split
                    other = n//factor
                    degree = factor.bit_length()+other.bit_length()-2
                    stats['faithful_lift'] = (len(coefficients) == degree+1 and all(
                        v == _raw_column(factor, other, j) for j, v in enumerate(coefficients)))
                    return factor
    return None
