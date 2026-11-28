"""Check the triangular ideals and replay unlabelled lists on planted factors."""

from itertools import product
from math import gcd, isqrt, prod
from random import Random

from sympy import Matrix, Poly, integer_nthroot, nextprime, primerange

if __package__:
    from .list_order_window import T, _basis, polynomial_list_divisors, residue_list_divisors
else:
    from list_order_window import T, _basis, polynomial_list_divisors, residue_list_divisors


def check_basis():
    equations = ((77,Poly(T,T)),(9,Poly(T**2+3*T-1,T)),(13,Poly(T**3+5,T)))
    d, radius = 8, 3
    count = 0
    for center, weights in product((-2,0,4), ((3,2,1),(2,0,2),(1,3,0))):
        rows = _basis(equations,weights,center,radius,d)
        determinant = radius**(d*(d-1)//2) * prod(
            mod**(f.degree()*t*(t+1)//2) for (mod,f),t in zip(equations,weights))
        assert Matrix(rows).det() == determinant
        for row in rows:
            polynomial = Poly.from_list([v//radius**j for j,v in enumerate(row)][::-1],T)
            for (mod,f),t in zip(equations,weights):
                quotient = polynomial
                for power in range(t,0,-1):
                    quotient, rem = quotient.div(f.shift(center),auto=False)
                    assert all(int(c) % mod**power == 0 for c in rem.all_coeffs())
                    count += 1
    return count


def binary_fixture(k):
    rng = Random(35361+k)
    moduli = list(primerange(257,700))[:k]
    combined = prod(moduli)
    scale = int(integer_nthroot(combined**3,5)[0])
    p = int(nextprime(scale+rng.randrange(scale//16)))
    q = int(nextprime(4*scale//3+rng.randrange(scale//20)))
    n, hints = p*q, []
    for modulus in moduli:
        true = p % modulus
        decoy = rng.randrange(modulus-1)
        decoy += decoy >= true
        residues = [true,decoy]
        rng.shuffle(residues)
        hints.append((modulus,residues))
    stats = {}
    assert residue_list_divisors(n,hints,max_dimension=12,stats=stats) == (p,)
    # Even after centering the full balanced box, its radius exceeds the
    # list-only degree-two determinant threshold sqrt(product of moduli).
    lower, upper = isqrt((n-1)//2)+1, isqrt(2*n)
    assert ((upper-lower+1)//2)**2 > combined
    assert all(gcd(n,modulus) == 1 for modulus in moduli)
    return n.bit_length(), combined.bit_length(), stats['dimension'], p


def check_local_lists():
    # With a composite local modulus, polynomial roots can mix listed residues.
    assert polynomial_list_divisors(77,[(15,(2,-3,1))],center=7,radius=0) == (7,)
    assert residue_list_divisors(77,[(15,(1,2))],center=7,radius=0) == ()
    # Two different local degrees and all balanced divisors of each toy n.
    count = 0
    for p,q in ((17,23),(29,31),(43,59)):
        n = p*q
        for r in (p,q):
            hints = [(1009,(r,1)), (1013,(r,11,13))]
            stats = {}
            got = residue_list_divisors(n,hints,center=(p+q)//2,radius=q-p,stats=stats)
            expected = tuple(v for v in range(2,n) if n % v == 0 and
                             all(v % m in values for m,values in hints))
            assert got == expected == (r,)
            count += 1
    for hints in ([(6,(1,1)),(9,(1,1))], [(7,(0,1))], [(5,(1,2))]):
        try:
            polynomial_list_divisors(77,hints,max_dimension=3)
        except ValueError:
            pass
        else:
            raise AssertionError('expected invalid moduli or leading coefficient')
    return count


def corrupted_fixture(count, bad):
    rng = Random(7214+count)
    moduli = list(primerange(257,700))[:count]
    combined = prod(moduli)
    base = int(integer_nthroot(combined,4)[0])
    p = int(nextprime(base+rng.randrange(base//16)))
    q = int(nextprime(4*base//3+rng.randrange(base//20)))
    n = p*q
    wrong = set(rng.sample(range(count),bad))
    # This power of two follows solely from count and the number of allowed
    # errors. Its factorization reveals none of the matching odd moduli.
    agreement = 2**(8*(count-bad))
    hints = []
    for i, mod in enumerate(moduli):
        true = p % mod
        decoys = rng.sample([v for v in range(mod) if v != true],2)
        hints.append((mod, (decoys[0], decoys[1] if i in wrong else true)))
    # Shuffle entries independently after planting so their position is no hint.
    hints = [(mod, values if rng.randrange(2) else values[::-1]) for mod,values in hints]
    assert gcd(agreement,combined) == 1
    assert {i for i,(mod,values) in enumerate(hints) if p % mod not in values} == wrong
    stats = {}
    assert residue_list_divisors(n,hints,min_agreement=agreement,max_dimension=12,stats=stats) == (p,)
    assert len(set(stats['multiplicities'][1:])) == 1
    assert stats['multiplicities'][1] > 0
    assert residue_list_divisors(n,hints,max_dimension=12) == ()
    return n, p, agreement, stats['dimension'], stats['multiplicities'][:2]


def check_enumerated_corruption():
    p, q = 29, 41
    n = p*q
    moduli = (101,103,107,109)
    choices = ((0,1),(p,1),(0,q),(p,q))
    agreement = min(moduli)**2
    counts = [0,0,0]
    for lists in product(choices,repeat=len(moduli)):
        hints = tuple(zip(moduli,lists))
        expected = tuple(r for r in range(p,q+1) if n % r == 0 and
                         prod(m for m,values in hints if r % m in values) >= agreement)
        stats = {}
        got = residue_list_divisors(n,hints,center=(p+q)//2,radius=(q-p)//2,
                                    max_dimension=8,min_agreement=agreement,stats=stats)
        assert got == expected and stats['dimension'] >= 2
        counts[len(got)] += 1
    assert all(counts)
    # The absence of a public certificate is a different outcome from ().
    try:
        residue_list_divisors(n,hints,min_agreement=agreement,max_dimension=2)
    except ValueError as error:
        assert 'certificate' in str(error)
    else:
        raise AssertionError('expected the deliberately insufficient cap to fail')
    return counts


def check_partial_gcd():
    # At r=7 the first polynomial contributes gcd(15,6)=3, the second
    # contributes its whole modulus 13. Literal lists count only the 13.
    equations = [(15,(-1,1)),(13,(-7,1))]
    for threshold, expected in ((13,(7,)),(39,(7,)),(40,())):
        assert polynomial_list_divisors(77,equations,center=7,radius=0,
                                        min_agreement=threshold) == expected
    lists = [(15,(1,)),(13,(7,))]
    assert residue_list_divisors(77,lists,center=7,radius=0,min_agreement=13) == (7,)
    assert residue_list_divisors(77,lists,center=7,radius=0,min_agreement=14) == ()
    for invalid in (0,196):
        try:
            polynomial_list_divisors(77,equations,min_agreement=invalid)
        except ValueError:
            pass
        else:
            raise AssertionError('expected an out-of-range agreement bound to fail')


if __name__ == '__main__':
    print('Local ideal divisibility checks:', check_basis())
    print('Exhaustively matched mixed-degree cases:', check_local_lists())
    for k in (4,8,16,32):
        bits, hint_bits, dimension, factor = binary_fixture(k)
        print(f'{k} binary lists, {2**k:,} CRT combinations, N={bits} bits, '
              f'combined modulus={hint_bits} bits, dimension={dimension}, factor={factor}')
    check_partial_gcd()
    print('Exhaustive corrupted-list outcomes with zero, one or two matches:', check_enumerated_corruption())
    for count, bad in ((16,5),(24,8)):
        n, factor, agreement, dimension, weights = corrupted_fixture(count,bad)
        print(f'{bad}/{count} corrupted binary lists, N={n}, H=2^{agreement.bit_length()-1}, '
              f'dimension={dimension}, multiplicities={weights}, factor={factor}')
