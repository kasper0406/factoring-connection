"""Factor from homogeneous quadratic character constraints without enumerating all CRT classes."""

from bisect import bisect_left
from itertools import groupby
from math import gcd, isqrt, prod

if __package__:
    from .cfrac_character_rows import character_word
    from .sqrt_prime_marks import _DecimalRing
else:
    from cfrac_character_rows import character_word
    from sqrt_prime_marks import _DecimalRing
from src.trackC.polyback import multipoint_eval, product_tree_linear


def _unique(values):
    return [u for u, _ in groupby(sorted(values))]


def _group_factor(n, offsets, points, work):
    """Detect proper divisors among all offsets minus points, ignoring full equalities."""
    roots = _unique(u % n for u in offsets)
    values = _unique(v % n for v in points)
    work["polynomial_coefficients"] += len(roots)+1
    work["evaluation_points"] += len(values)
    ring = _DecimalRing(n)
    poly = product_tree_linear(ring, [(-u) % n for u in roots])[-1][0]
    matches, others = [], []
    for v in values:
        i = bisect_left(roots,v)
        (matches if i < len(roots) and roots[i] == v else others).append(v)
    work["derivative_points"] += len(matches)
    derivative = ring.make([i*c for i,c in enumerate(poly) if i]) if matches else []
    for f, positions in ((poly,others),(derivative,matches)):
        for v,value in zip(positions, multipoint_eval(ring, f, positions)):
            d = gcd(value,n)
            if 1 < d < n:
                return d
            if d == n:
                work["saturated"] += 1
                for u in roots:
                    d = gcd(u-v,n)
                    if 1 < d < n:
                        return d
                raise ArithmeticError("saturated product without a proper difference")
    return None


def _simple_basis(rows):
    pivots = {}
    for mask in rows:
        while mask:
            j = mask.bit_length()-1
            if j not in pivots:
                pivots[j] = mask
                break
            mask ^= pivots[j]
    return pivots


def _project(rows, side):
    """Equations supported entirely on side, including linear combinations of input rows."""
    outer, inside = {}, []
    for row in rows:
        while row & ~side:
            j = (row & ~side).bit_length()-1
            if j not in outer:
                outer[j] = row
                break
            row ^= outer[j]
        else:
            inside.append(row)
    return _simple_basis(inside)


def _half_classes(blocks, rows, modulus):
    """Admissible CRT contributions indexed by the shared partial parity of the rows."""
    side = sum(mask for _,mask,_ in blocks)
    local = _project(rows, side)
    free = [i for i in range(side.bit_length()) if side & (1 << i) and i not in local]
    base_size = prod(p for p,_,_ in blocks)
    idempotent = (modulus//base_size)*pow(modulus//base_size, -1, base_size)
    result = []
    for number in range(1 << len(free)):
        word = sum(((number >> j)&1) << i for j,i in enumerate(free))
        for pivot, mask in sorted(local.items()):
            word |= ((word & mask).bit_count() % 2) << pivot
        syndrome = sum(((word & row).bit_count()%2) << i for i,row in enumerate(rows))
        representatives, current = [0], 1
        for prime, mask, choices in blocks:
            allowed = choices[word & mask]
            inverse = pow(current,-1,prime)
            representatives = [u+current*((v-u)*inverse % prime)
                               for u in representatives for v in allowed]
            current *= prime
        result.extend((syndrome, u*idempotent % modulus) for u in representatives)
    result.sort()
    return [(key,[u for _,u in group]) for key,group in groupby(result,key=lambda item:item[0])]


def _local_blocks(primes):
    bins8 = {character_word(r,()):[r] for r in (1,3,5,7)}
    blocks = [(8,3,bins8)]
    for i,p in enumerate(primes,2):
        bit = 1 << i
        squares = {k*k % p for k in range(1,(p+1)//2)}
        nonsquares = [k for k in range(1,p) if k not in squares]
        blocks.append((p,bit,{0:sorted(squares),bit:nonsquares}))
    return blocks


def character_wheel_factor(n, constraints, *, stats=None):
    """Return a proper divisor when the character rows hold at the smaller balanced prime.

    constraints is the CharacterRows result of collect_character_rows; its conductor is
    factored and its equations are homogeneous. A missed class may yield an incidental gcd.
    """
    if n < 15 or n % 2 == 0:
        raise ValueError("expected an odd balanced semiprime")
    if constraints.divisor is not None:
        if 1 < constraints.divisor < n and n % constraints.divisor == 0:
            return constraints.divisor
        raise ValueError("unverified supplied divisor")
    m = constraints.conductor
    if gcd(m,n) != 1:
        for p in constraints.primes:
            if n % p == 0:
                return p
    rows = [mask for mask,_ in constraints.rows]
    sides, sizes = [[],[]], [1,1]
    for block in sorted(_local_blocks(constraints.primes), reverse=True):
        which = int(sizes[1] < sizes[0])
        sides[which].append(block)
        sizes[which] *= sum(map(len,block[2].values()))
    left,right = (_half_classes(blocks, rows, m) for blocks in sides)
    work = {} if stats is None else stats
    work.update(projected_left=sum(len(items) for _,items in left),
                projected_right=sum(len(items) for _,items in right),
                syndromes=len(left), polynomials_tested=0, polynomial_coefficients=0,
                evaluation_points=0, derivative_points=0, saturated=0)
    lower, upper = isqrt((n-1)//2)+1, isqrt(n)
    # 0 <= u+v <= 2*m-2. Every p=u+v+t*m in [lower,upper] is in this lift range.
    start, stop = -((2*m-2-lower)//m), upper//m
    lifts = stop-start+1
    work["lifts"] = lifts
    for (key,A), (other,B) in zip(left,right,strict=True):
        if key != other:
            raise ArithmeticError("inconsistent projections of the linear equations")
        stride = min(lifts,max(1,isqrt(len(B)*lifts//len(A))))
        offsets = (a+i*m for a in A for i in range(stride))
        points = (-b-t*m for b in B for t in range(start,stop+1,stride))
        work["polynomials_tested"] += 1
        divisor = _group_factor(n,offsets,points,work)
        if divisor is not None:
            return divisor
    return None
