"""Independent recurrence checks plus public recovery of balanced factors."""

from itertools import combinations, product

from sympy import Matrix, isprime, primerange

if __package__:
    from .carry_recurrence import _extend, _fit, _high_seeds, _low_seeds, recurrent_carry_factor
    from .carry_polynomial import _read_factor
else:
    from carry_recurrence import _extend, _fit, _high_seeds, _low_seeds, recurrent_carry_factor
    from carry_polynomial import _read_factor


def schoolbook(p, q):
    raw = [0]*(p*q).bit_length()
    for j in range(p.bit_length()):
        for i in range(q.bit_length()):
            raw[i+j] += ((p >> j) & 1) * ((q >> i) & 1)
    carries = []
    incoming = 0
    for j, value in enumerate(raw):
        incoming, digit = divmod(value+incoming, 2)
        assert digit == ((p*q >> j) & 1)
        carries.append(incoming)
    assert carries.pop() == 0
    return carries


def direct_order(row):
    if not any(row):
        return 0
    for k in range(1, len(row)):
        lhs = Matrix([row[i-k:i] for i in range(k,len(row))])
        if lhs.rank() == lhs.row_join(Matrix(row[k:])).rank():
            return k
    return len(row)


def check_fits():
    for row in product(range(2), repeat=6):
        assert direct_order(row) == len(_fit(row)[0])
    assert _extend((8,4), 1, 12) == (8,4,2,1)
    assert _extend((0,0), 1, 5) == (0,)*5
    # The final two terms of the known seed can be identical without
    # making the extrapolation ambiguous.
    periodic = (2,2,0,2,2,0)
    assert _extend(periodic, 3, 12) == periodic+periodic
    return 67


def check_small():
    pairs = cases = promised = 0
    for p, q in combinations(map(int, primerange(11, 64)), 2):
        if q > 2*p:
            continue
        pairs += 1
        n = p*q
        row = schoolbook(p,q)
        for k in (1,2,3):
            size = min(2*k,len(row))
            assert tuple(row[:size]) in _low_seeds(n,size)
            assert tuple(row[-size:][::-1]) in _high_seeds(n,size)
            eligible = any(max(len(_fit(row[:h])[0]), len(_fit(row[h:][::-1])[0])) <= k
                           for h in range(len(row)+1))
            result = recurrent_carry_factor(n,k)
            cases += 1
            if eligible:
                assert result in (p,q)
                promised += 1
            assert result is None or result in (p,q)
    return pairs, cases, promised


def check_fixtures():
    # Selected after scanning short periodic digit patterns.
    p, q = 78643, 87383
    assert isprime(p) and isprime(q)
    n = p*q
    row = schoolbook(p,q)
    assert len(_fit(row[:15])[0]) <= 5 and len(_fit(row[15:][::-1])[0]) <= 5
    assert _read_factor(n, [(n >> j) & 1 for j in range(n.bit_length())]) is None
    stats = {}
    assert recurrent_carry_factor(n,5,stats=stats) in (p,q)
    assert stats['faithful_lift']
    print("Faithful split:", n, p, q, "carry sum", sum(row), stats)

    p, q = 1789569707, 2147483647
    assert isprime(p) and isprime(q)
    n = p*q
    row = schoolbook(p,q)
    assert _fit(row[:30]) == _fit(row[30:][::-1]) == ((1,1,-1),1)
    assert recurrent_carry_factor(n,3,stats=stats) in (p,q)
    assert not stats['faithful_lift']
    print("Split already available with zero carries:", n, "true carry sum", sum(row), stats)

    # For every tested odd m, the number of changes of slope grows with m
    # although the order-three certificate stays fixed.
    for m in range(7, 128, 2):
        a, b = (5*(1 << (m-1))+1)//3, (1 << m)-1
        row = schoolbook(a,b)
        assert row[:m-1] == [j//2 for j in range(1,m)]
        right = tuple(row[m-1:][::-1])
        assert right == tuple(1+j//2 for j in range(1,m+1))
        assert right == _extend(right[:6],3,m)


if __name__ == '__main__':
    print("Recurrence controls:", check_fits())
    print("Small pairs, public trials, promised cases:", *check_small())
    check_fixtures()
