"""Validate changes of carries and slopes. --quick skips the 85-bit profile search."""

from argparse import ArgumentParser
from itertools import combinations_with_replacement
from math import isqrt

from sympy import isprime, primerange

if __package__:
    from .carry_polynomial import _uncarry, carry_factor
else:
    from carry_polynomial import _uncarry, carry_factor


def schoolbook(a, b):
    n = a*b
    raw = [0]*n.bit_length()
    for i in range(a.bit_length()):
        if a >> i & 1:
            for j in range(b.bit_length()):
                raw[i+j] += b >> j & 1
    carries = []
    prev = 0
    for i, value in enumerate(raw):
        prev, digit = divmod(value+prev, 2)
        assert digit == (n >> i & 1)
        carries.append(prev)
    assert prev == 0
    assert sum(carries) == a.bit_count()*b.bit_count()-n.bit_count()
    return raw, carries


def differences(carries, order):
    values = carries[:]
    for _ in range(order):
        values = [new-old for new, old in zip(values+[0], [0]+values)]
    return tuple((i, v) for i, v in enumerate(values) if v)


def check_small():
    facts = checks = 0
    for p, q in combinations_with_replacement(map(int, primerange(3, 44)), 2):
        n = p*q
        raw, carries = schoolbook(p, q)
        while not raw[-1]:
            raw.pop()
        bits = [n >> i & 1 for i in range(n.bit_length())]
        for order in (1, 2):
            jumps = differences(carries, order)
            assert _uncarry(bits, jumps, order) == raw
            checks += 1
            if len(jumps) <= 4 and all(abs(v) <= 2 for _, v in jumps):
                assert carry_factor(n, len(jumps), order=order, jump_bound=2) in (p, q)
                facts += 1
    assert carry_factor(131, 3) is None
    assert carry_factor(22, 0) == 2
    print(f"{checks} exact schoolbook reconstructions and {facts} promised input-only splits")


def check_families():
    # A has two adjacent high bits and separated lower bits including bit 0.
    # B contributes a single low interval and the translate by 2^height.
    for span in (4, 9, 16, 32, 64):
        positions = [j*(span+2) for j in range(span)]
        height = (span+1)*(span+2)
        a = (3 << (height-1)) + sum(1 << i for i in positions)
        b = (1 << height) + (1 << span)-1
        raw, carries = schoolbook(a, b)
        assert differences(carries, 1) == ((height, 1), (height+span, -1))
        assert sum(carries) == sum(max(0, value-1) for value in raw) == span
        assert (a.bit_count(), b.bit_count()) == (span+2, span+1)
    for a, b in ((4, 4), (8, 13), (31, 32), (127, 129)):
        _, carries = schoolbook(2**a-1, 2**b-1)
        assert differences(carries, 2) == ((1, 1), (a, -1), (b+1, -1), (a+b, 1))
        assert max(carries) == a-1 and sum(carries) == b*(a-1)


def check_fixture(p, q, order, expected):
    assert p < 2**64 and q < 2**64 and isprime(p) and isprime(q)
    n = p*q
    raw, carries = schoolbook(p, q)
    assert differences(carries, order) == expected
    assert carry_factor(n, 0) is None
    work = {}
    result = carry_factor(n, len(expected), order=order, jump_bound=1, stats=work)
    assert result in (p, q) and work['winning_jumps'] == expected
    fermat_offset = (p+q)//2-isqrt(n)-1
    print(f"{n.bit_length()} bits, factor {result}, difference order {order}, "
          f"total carries {sum(carries)}, peak {max(carries)}, "
          f"raw excess {sum(max(0, v-1) for v in raw)}, "
          f"{work['profiles']} profiles, {work['factorizations']} polynomial factorizations, "
          f"Fermat offset {fermat_offset}")


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true')
    args = parser.parse_args()
    check_small()
    check_families()
    check_fixture(1048571, 1048573, 2, ((3, 1), (20, -1), (23, -1), (40, 1)))
    if not args.quick:
        check_fixture(4398046773247, 6597070815233, 1, ((42, 1), (60, -1)))
