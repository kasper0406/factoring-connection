"""Exhaust cubic characters and compare exact marked splitting probabilities to extracted gcds."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import gcd, lcm

from sympy import primerange

if __package__:
    from .eisenstein_torus import ONE, OMEGA, cubic_jacobi, norm, torus_power, torus_product
    from .marked_torus_ladder import torus_cubic_split
else:
    from eisenstein_torus import ONE, OMEGA, cubic_jacobi, norm, torus_power, torus_product
    from marked_torus_ladder import torus_cubic_split


def order_parts(r):
    epsilon = 1 if r % 3 == 1 else -1
    part = r-epsilon
    height = 0
    while part % 3 == 0:
        height += 1
        part //= 3
    return epsilon, height, part


def tag_law(a, b, d, mark):
    if mark:
        return {(a, d*pow(b,-1,3)*mark % 3): Fraction(1)}
    size = 3**(a-1)
    law = {(0,0): Fraction(1,size)}
    for j in range(1,a):
        for sign in (1,2):
            law[j,sign] = Fraction(3**(j-1),size)
    assert sum(law.values()) == 1
    return law


def predicted(p, q, mark, bound):
    ep, ap, bp = order_parts(p)
    eq, aq, bq = order_parts(q)
    g = gcd(bp,bq)
    m = lcm(*(t for t in range(1,bound+1) if t % 3))
    ratio_p, ratio_q = bp//g, bq//g
    hit_p = Fraction(gcd(ratio_p,m),ratio_p)
    hit_q = Fraction(gcd(ratio_q,m),ratio_q)
    d = p*q-ep*eq
    while d % 3 == 0:
        d //= 3
    d *= m
    fail_both = 0
    for u in range(3):
        v = (-mark-ep*u)*eq % 3
        left = tag_law(ap,bp,d,u)
        right = tag_law(aq,bq,d,v)
        fail_both += sum(x*right.get(key,0) for key,x in left.items())/3
    return hit_p + hit_q - hit_p*hit_q*(1+fail_both)


def predicted_euler(p,q):
    ep, eq = (1 if r % 3 == 1 else -1 for r in (p,q))
    mp,mq = (p-ep)//3, (q-eq)//3
    g = gcd(mp,mq)
    hp,hq = Fraction(g,mp), Fraction(g,mq)
    e = (p*q-ep*eq)//3
    same = Fraction(1,3) if (e//g) % 3 else Fraction(1)
    return hp+hq-hp*hq*(1+same)


def all_elements(r):
    return [(a,b) for a,b in product(range(r),repeat=2) if norm((a,b)) % r == 1]


def crt_element(p, q, left, right):
    inverse = pow(p, -1, q)
    return tuple(x + p*((y-x)*inverse % q) for x, y in zip(left, right))


def marked_fibers(p, q):
    by_mark = {mark: [] for mark in range(3)}
    for u, v in product(all_elements(p), all_elements(q)):
        z = crt_element(p, q, u, v)
        by_mark[cubic_jacobi(p*q, z)].append(z)
    return by_mark


def enters_cubic_roots(r, value, exponent):
    roots = (ONE, OMEGA, (r-1, r-1))
    z = torus_power(r, value, exponent)
    for _ in range(r.bit_length()+1):
        if z in roots:
            return True
        z = torus_power(r, z, 3)
    return False


def check_word_law(p, q, bound, words):
    """Enumerate local word distributions inside every public mark vector."""
    n = p*q
    fibers = marked_fibers(p, q)
    multiplier = lcm(*(t for t in range(1, bound+1) if t % 3))
    ep, _, bp = order_parts(p)
    eq, _, bq = order_parts(q)
    exponent = n-ep*eq
    while exponent % 3 == 0:
        exponent //= 3
    exponent *= multiplier
    cases = 0
    for r, sign, b in ((p,ep,bp), (q,eq,bq)):
        remaining = b//gcd(bp,bq)
        deficit = remaining//gcd(remaining,multiplier)
        local_order = r-sign
        entries = {z: enters_cubic_roots(r,z,exponent) for z in all_elements(r)}
        for word in words:
            expected = Fraction(gcd(deficit,*word),deficit)
            # Keep each multiplicity from enumerating full public-element fibers.
            terms = [{m: Counter(torus_power(r,z,v % local_order) for z in fiber)
                      for m, fiber in fibers.items()} for v in word]
            for marks in product(range(3), repeat=len(word)):
                distribution = Counter({ONE: 1})
                for term, mark in zip(terms, marks):
                    combined = Counter()
                    for x, count_x in distribution.items():
                        for y, count_y in term[mark].items():
                            combined[torus_product(r,x,y)] += count_x*count_y
                    distribution = combined
                hits = sum(count for z,count in distribution.items() if entries[z])
                assert Fraction(hits,distribution.total()) == expected, (n,r,bound,word,marks)
                cases += 1
    return cases


def check_143():
    n = 143
    fibers = marked_fibers(11,13)
    assert [len(fibers[m]) for m in range(3)] == [48,48,48]
    for mark, fiber in fibers.items():
        for z in fiber:
            assert torus_power(n,z,48) == ONE
            assert torus_cubic_split(n,z,full_ladder=False) is None
            divisor = torus_cubic_split(n,z)
            if mark:
                assert divisor in (11,13)
            else:
                assert divisor is None
    z = (1,66)
    assert cubic_jacobi(n,z) == 1 and torus_cubic_split(n,z) == 11
    endpoint = torus_power(n, z, 16)
    square = torus_power(n, endpoint, 2)
    projector = tuple(
        (one + x + y) * pow(3, -1, n) % n
        for one, x, y in zip(ONE, endpoint, square)
    )
    assert projector == (78, 0)
    root = (2 * projector[0] - 1) % n
    assert root == 12 and root**2 % n == 1 and gcd(root - 1, n) == 11


def exhaustive():
    primes = list(primerange(5,44))
    points = {r: all_elements(r) for r in primes}
    count = cases = 0
    for p,q in combinations(primes,2):
        n = p*q
        hits = Counter()
        marks = Counter()
        for left,right in product(points[p],points[q]):
            z = crt_element(p, q, left, right)
            mark = cubic_jacobi(n,z)
            marks[mark] += 1
            for kind,bound in [('euler',1),('ladder',1),('ladder',5)]:
                f = torus_cubic_split(n,z,full_ladder=kind=='ladder',bound=bound)
                if f is not None:
                    assert f in (p,q)
                    hits[kind,bound,mark] += 1
            count += 1
        assert len(set(marks.values())) == 1
        for mark in range(3):
            assert Fraction(hits['euler',1,mark],marks[mark]) == predicted_euler(p,q)
            for bound in (1,5):
                assert Fraction(hits['ladder',bound,mark],marks[mark]) == predicted(p,q,mark,bound)
        cases += 1
    print(cases, "semiprimes and", count, "torus elements verified")
    for p,q in ((7,13),(11,23),(13,29),(17,19),(19,37),(23,31),(29,43)):
        print((p,q), "Euler", predicted_euler(p,q), "ladder by mark",
              tuple(str(predicted(p,q,j,1)) for j in range(3)))


if __name__ == '__main__':
    exhaustive()
    check_143()
    # d_23=4 at B=1 and d_23=2 at B=2. For N=11*29 both deficits exceed 1.
    words = ((1,-2), (2,-4), (0,4), (5,-10), (0,0), (2,0,-3))
    cases = sum(check_word_law(p,q,bound,words) for p,q,bound in
                [(5,23,1), (5,23,2), (11,29,1)])
    print(cases, "local word probabilities verified across all public mark vectors")
