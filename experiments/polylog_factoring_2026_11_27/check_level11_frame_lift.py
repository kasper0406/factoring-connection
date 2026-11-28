"""Check both source and target classes of the conditional intermediate-frame decoder."""

from collections import Counter
from random import Random

from sympy import nextprime

if __package__:
    from .level11_frame_lift import factor_from_frame, frame_lift_exponent
    from .level11_marking import framed_line_classifier
else:
    from level11_frame_lift import factor_from_frame, frame_lift_exponent
    from level11_marking import framed_line_classifier


def run():
    class_pairs = Counter()
    exhaustive = 0
    for p, q in ((3, 5), (5, 7), (7, 13), (13, 19), (19, 23)):
        n = p*q
        k = frame_lift_exponent(n)
        for r in (p, q):
            for base in ("principal", "nonprincipal"):
                # Only the fixture generator is told which prime produces
                # the intermediate frame. The decoder receives n and the frame.
                classify = framed_line_classifier(r, k, Random(r), base=base)
                for u, v in [(1, s) for s in range(r)] + [(0, 1)]:
                    frame = classify(u, v)
                    assert factor_from_frame(n, k, frame, base=base) == p
                    class_pairs[base, frame[0]] += 1
                    exhaustive += 1
    assert len(class_pairs) == 4

    large = largest_bits = 0
    for bits in (31, 63, 127):
        p = int(nextprime(1 << bits))
        q = int(nextprime(p + 100))
        n = p*q
        largest_bits = max(largest_bits, n.bit_length())
        k = frame_lift_exponent(n)
        assert (1 << (4*(k-1))) < n <= (1 << (4*k))
        for r in (p, q):
            for base in ("principal", "nonprincipal"):
                classify = framed_line_classifier(r, k, Random(r), base=base)
                frame = classify(1, 37)
                assert factor_from_frame(n, k, frame, base=base) == p
                large += 1

    # A public whole-norm classifier answers a different question. On
    # this input its determinant sector excludes both intermediate primes.
    n = 19*23
    k = frame_lift_exponent(n)
    assert all(n % (1 << k) != r % (1 << k) for r in (19, 23))
    for base in ("principal", "nonprincipal"):
        endpoint = framed_line_classifier(n, k, Random(53), base=base)
        for u, v in ((0, 1), (1, 0), (1, 7)):
            assert factor_from_frame(n, k, endpoint(u, v), base=base) is None

    print(f"{exhaustive} exhaustive intermediate frames across both classes at both ends")
    print(f"{large} larger intermediate frames through {largest_bits} bits")
    print("6 whole-modulus frames rejected by the lift")


if __name__ == "__main__":
    run()
