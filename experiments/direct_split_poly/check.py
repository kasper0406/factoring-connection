"""Factor-aware diagnostics for two cutoff-product transfer representations."""

import json
from math import gcd, isqrt

from flint import fmpz, nmod_mat


def prefix_products(modulus, limit):
    values = [1]
    for k in range(1, limit + 1):
        values.append(values[-1] * k % modulus)
    return values


def packed_comparator(k, border_bits, modulus):
    less, equal = 0, 1
    for digit, border in zip(f"{k:0{len(border_bits)}b}", border_bits):
        if digit == "0":
            less = (less + border * equal) % modulus
            equal = (1 - border) * equal % modulus
        else:
            equal = border * equal % modulus
    return less


def cyclic_middle_diagnostic(n, cap, max_width):
    # Materialize exact coefficients for the diagnostic. Their simple
    # recurrence divides ordinary integers exactly, never residues modulo N.
    coefficients = [1]
    for k in range(1, 2 * cap + 1):
        coefficients.append(coefficients[-1] * (2 * cap - k + 1) // k)
    middle_gcd = gcd(coefficients[cap], n)
    nonunits = {}
    for width in range(1, max_width + 1):
        position = cap % width
        folded_gcd = gcd(sum(coefficients[position::width]), n)
        if folded_gcd != 1:
            nonunits[width] = folded_gcd
    return middle_gcd, nonunits


def diagnostic(p, q):
    assert fmpz(p).is_prime() and fmpz(q).is_prime() and p < q < 2 * p
    n = p * q
    cap = isqrt(n)
    cutoff_gcd = gcd(prefix_products(n, cap)[-1], n)
    assert cutoff_gcd == p
    values = prefix_products(p, cap)

    size = min(p, cap - p + 2)
    start = p - size
    labels = [1 + (k * k + 7 * k + 5) % (p - 1) if k < p else 0 for k in range(cap + 1)]
    unary_ranks = []
    for sequence in (values, labels):
        block = [sequence[start + i : start + i + size] for i in range(size)]
        unary_ranks.append(nmod_mat(block, p).rank())
    assert unary_ranks == [size, size]

    m = 1 << (isqrt(min(p, cap + 1)).bit_length() - 1)
    block = [[values[m * j + i] for j in range(m)] for i in range(m)]
    matrix = nmod_mat(block, p)
    determinant = 1
    for j in range(m):
        determinant = determinant * values[m * j] % p
        for a in range(j):
            determinant = determinant * m * (j - a) % p
    assert determinant and int(matrix.det()) == determinant and matrix.rank() == m

    # The local threshold bits require the factors. This is a control, never
    # an operational compiler. Its only inverse is modulo the known prime q.
    bits = q.bit_length()
    inverse = pow(p, -1, q)
    packed = []
    for pos in reversed(range(bits)):
        bp, bq = (p >> pos) & 1, (q >> pos) & 1
        packed.append(bp + p * ((bq - bp) * inverse % q))
    for k in range(1 << bits):
        value = packed_comparator(k, packed, n)
        assert (value % p, value % q) == (int(k < p), int(k < q))
    bit_divisors = [gcd(b, n) for b in packed]
    proper = [d for d in bit_divisors if 1 < d < n]
    assert proper and all(d in (p, q) for d in proper)
    assert gcd(packed_comparator(cap, packed, n), n) == p
    fold_limit = min(cap, 128)
    middle_gcd, nonunit_folds = cyclic_middle_diagnostic(n, cap, fold_limit)
    assert middle_gcd == q
    if (p, q) in ((503, 937), (1009, 1877)):
        assert not nonunit_folds
    return {
        "N": n,
        "bits": n.bit_length(),
        "p": p,
        "q": q,
        "cap": cap,
        "cutoff_gcd": cutoff_gcd,
        "unary_width_bound": size,
        "unary_ranks": unary_ranks,
        "binary_exact_rank": m,
        "binary_determinant_mod_p": determinant,
        "first_packed_bit_divisor": proper[0],
        "exact_middle_gcd": middle_gcd,
        "max_fold_width": fold_limit,
        "nonunit_folds": nonunit_folds,
    }


if __name__ == "__main__":
    pairs = [(101, 103), (101, 179), (251, 467), (503, 937), (1009, 1877)]
    for p, q in pairs:
        print(json.dumps(diagnostic(p, q), separators=(",", ":")))
