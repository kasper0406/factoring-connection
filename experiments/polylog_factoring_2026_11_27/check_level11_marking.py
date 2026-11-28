"""Exhaustive checks of the five-state level-11 marking and its prime marginals."""

from collections import Counter
from fractions import Fraction
from itertools import product
from random import Random

if __package__:
    from . import level11_marking as marking
    from .check_level11_marginals import eta_coefficients, image_line, norm_points
else:
    import level11_marking as marking
    from check_level11_marginals import eta_coefficients, image_line, norm_points


def trace_44(r):
    # LMFDB curve 44.a1: y² = x³ + x² - 77x - 289.
    total = 0
    for x in range(r):
        rhs = (x**3 + x*x - 77*x - 289) % r
        if rhs:
            total += -1 if pow(rhs, (r-1)//2, r) == 1 else 1
    return total


def direct_witnesses(n, factors, rows):
    """Enumerate every integral witness independently of the ideal shell solver."""
    by_cell = {}
    for principal, shell in ((True, n), (False, 3*n)):
        for point in norm_points(shell):
            if not principal:
                # For (3,n)=1, I*bar(J) = I ∩ bar(J): the two
                # congruences have disjoint prime supports and equal index.
                x, y, z, w = marking._conjugate(point)
                if (x+y+z) % 3 or (z+w-x) % 3:
                    continue
            cell = tuple(image_line(point, rows, r) for r in factors)
            by_cell.setdefault(cell, []).append((point, principal))
    return by_cell


def direct_marks(n, factors, rows):
    by_cell = direct_witnesses(n, factors, rows)
    result = {}
    for cell in product(*(range(r+1) for r in factors)):
        marks = [marking._witness_mark(point, principal) for point, principal in by_cell[cell]]
        assert len(set(marks)) == 1 and len(marks) == (4 if marks[0] < 3 else 6)
        result[cell] = marks[0]
    return result


def expected(r, a11, source, target):
    # Projectors onto constants, the single A-orientation sign, and the
    # orthogonal three-dimensional old space give the three coefficients.
    sign = (0, 0, 0, 1, -1)
    return (Fraction(r+1, 5) + a11*(int(source == target) - Fraction(1, 5))
            + Fraction(trace_44(r)-a11, 2)*sign[source]*sign[target])


def crt_line(p, q, i, j):
    left = (1, i) if i < p else (0, 1)
    right = (1, j) if j < q else (0, 1)
    return tuple(u + p*((v-u)*pow(p, -1, q) % q) for u, v in zip(left, right))


def pair_check(p, q, coefficients):
    n = p*q
    inverse = marking._inverse_splitting(n, Random(n))
    forward = inverse.adjugate()*pow(int(inverse.det()), -1, n)
    rows = [[int(x) % n for x in row] for row in forward.tolist()]
    table = direct_marks(n, (p, q), rows)
    example = {}
    for axis, r, s in ((0, p, q), (1, q, p)):
        labels = direct_marks(r, (r,), rows)
        for (i,), local in labels.items():
            counts = Counter(mark for cell, mark in table.items() if cell[axis] == i)
            assert [counts[k] for k in range(5)] == [
                expected(s, coefficients[s], local, k) for k in range(5)]
            if r == 29 and local == 3:
                example = dict(sorted(counts.items()))

    # Run the public, factor-free classifier on the same coherent splitting.
    classifier = marking.marked_line_classifier(n, Random(n))
    for (i, j), mark in table.items():
        assert classifier(*crt_line(p, q, i, j)) == mark
    return len(table), example


def nonprincipal_pair(p, q, coefficients):
    n = p*q
    inverse = marking._inverse_splitting(n, Random(n))
    classify = marking.marked_line_classifier(n, Random(n), base="nonprincipal")
    table = {}
    for i, j in product(range(p+1), range(q+1)):
        table[i, j] = classify(*crt_line(p, q, i, j))
    for axis, r, s in ((0, p, q), (1, q, p)):
        for i in range(r+1):
            line = (1, i) if i < r else (0, 1)
            local = marking._marked_line(r, inverse, line, "nonprincipal")
            counts = Counter(mark for cell, mark in table.items() if cell[axis] == i)
            assert all(counts[k] == expected(s, coefficients[s], local, k) for k in range(5))
    counts = Counter(table.values())
    assert counts[3] - counts[4] == trace_44(p)*trace_44(q)
    return len(table), counts[3]-counts[4]


def direct_frame(point, principal, k):
    modulus = 1 << k
    s = 0
    while (s*s+s+4) % modulus:
        s = (s - (s*s+s+4)*pow(2*s+1, -1, modulus)) % modulus
    # Multiplication by omega, and j=(s+omega) times conjugation, on (1,omega).
    x, y, z, w = point
    raw = (x+s*z-3*w, -3*y+(s+3)*z+3*s*w,
           y+z+(s+1)*w, x+y-s*z+3*w)
    scale = 1 if principal else pow(3, -1, modulus)
    return tuple(a*scale % modulus for a in raw)


def check_frames(p, q, k, native_all):
    n = p*q
    inverse = marking._inverse_splitting(n, Random(n))
    forward = inverse.adjugate()*pow(int(inverse.det()), -1, n)
    rows = [[int(x) % n for x in row] for row in forward.tolist()]
    cells = direct_witnesses(n, (p,q), rows)
    assert len(cells) == (p+1)*(q+1)
    classify = marking.framed_line_classifier(n, k, Random(n))
    states = Counter()
    row_support = [set() for i in range(p+1)]
    column_support = [set() for j in range(q+1)]
    for (i,j), witnesses in cells.items():
        principal = witnesses[0][1]
        orbit = {direct_frame(point, flag, k) for point, flag in witnesses}
        assert len(orbit) == len(witnesses) == (4 if principal else 6)
        state = ("principal" if principal else "nonprincipal", min(orbit))
        if native_all or not states[state]:
            assert classify(*crt_line(p, q, i, j)) == state
        a,b,c,d = state[1]
        assert ((a*d-b*c)*(1 if principal else 3) - n) % (1 << k) == 0
        states[state] += 1
        row_support[i].add(state)
        column_support[j].add(state)
    # Quantify distance to the observed full law, independent of any mixing bound.
    supports = row_support + column_support
    largest = max(sum(states[t] for t in support) for support in supports)
    sector = 5*(1 << (3*k))//16
    full_tv = (sum(abs(Fraction(count, len(cells)) - Fraction(1,sector))
                   for count in states.values()) + Fraction(sector-len(states),sector))/2
    return len(states), sector, full_tv, 1-Fraction(largest,len(cells))


if __name__ == "__main__":
    coefficients = eta_coefficients(29)
    total = 0
    for p, q in ((5, 7), (7, 13), (13, 19), (19, 29)):
        cells, example = pair_check(p, q, coefficients)
        total += cells
        if (p,q) == (19,29):
            print("19*29: 29-line marked A+, counts", [example.get(i, 0) for i in range(5)])
    # The public construction also applies to the prime of the reference ideal.
    by_three = 0
    for n in (1, 3, 9, 15, 27):
        classifier = marking.marked_line_classifier(n, Random(n))
        for v in range(n):
            mark = classifier(1, v)
            assert 0 <= mark < 5
            by_three += 1
    print(f"{total} CRT cells checked against native markings, {by_three} extra native lines")
    for p, q in ((3, 5), (5, 7), (19, 29)):
        cells, signed = nonprincipal_pair(p, q, coefficients)
        print(f"{p}*{q}: A+ start, {cells} native cells, signed sum {signed}")
    for p, q, k in ((5,7,2), (5,7,7), (19,29,3), (73,79,3)):
        occupied, sector, distance, complement = check_frames(p, q, k, (p,q) != (73,79))
        print(f"{p}*{q}, modulus {1<<k}: {occupied}/{sector} states, "
              f"full TV {distance}, minimum full mass outside any marginal support {complement}")
    # Large k without listing the matrix group; also check the source norm factor.
    for base, factor in (("principal",1), ("nonprincipal",3)):
        n = 1000003*1000033
        classify = marking.framed_line_classifier(n, 128, Random(17), base=base)
        for u,v in ((1,7),(0,1)):
            target, (a,b,c,d) = classify(u,v)
            assert ((a*d-b*c)*(1 if target == "principal" else 3) - factor*n) % (1<<128) == 0
