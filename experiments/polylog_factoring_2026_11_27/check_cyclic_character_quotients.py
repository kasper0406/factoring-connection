"""Test character and low-order quotient shortcuts for the cyclic query family."""

from cyclic_class_queries import encode


def mixed_roots(p, q):
    n = p * q
    inverse = pow(p, -1, q)
    roots = {
        (left + p * ((right - left) * inverse % q)) % n
        for left in (1, p - 1)
        for right in (1, q - 1)
    }
    return sorted(roots - {1, n - 1})


def window_parity(roots, start, length, n):
    return sum((root - start - 1) % n < length for root in roots) % 2


def quotient_collision(p, q):
    n = p * q
    half = (n - 1) // 2
    roots = mixed_roots(p, q)
    values = []
    for start in range(half // 2 + 1):
        _, exponent, _ = encode(n, start)
        assert exponent % 2 == 1
        values.append((exponent, window_parity(roots, start, half, n), start))

    deepest = None
    for bits in range(1, 2 * n.bit_length() + 8):
        seen = {}
        collision = None
        for exponent, parity, start in values:
            residue = exponent % (1 << bits)
            if residue in seen and seen[residue][1] != parity:
                collision = (seen[residue], (exponent, parity, start))
                break
            seen[residue] = (exponent, parity, start)
        if collision is not None:
            deepest = (bits, collision)
    return n, deepest


if __name__ == "__main__":
    controls = [(101, 109), (127, 227), (503, 937)]
    results = [quotient_collision(p, q) for p, q in controls]
    expected = [(11009, 13), (28829, 15), (471311, 19)]
    assert [(n, result[0]) for n, result in results] == expected
    for n, (bits, ((left_e, left_v, left_a), (right_e, right_v, right_a))) in results:
        assert left_e % (1 << bits) == right_e % (1 << bits)
        assert left_v != right_v
        print(
            f"N={n}: quotient modulo 2^{bits} collides at starts "
            f"{left_a},{right_a} with parities {left_v},{right_v}"
        )
