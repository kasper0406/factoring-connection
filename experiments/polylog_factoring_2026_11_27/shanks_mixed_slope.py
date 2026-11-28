"""Mixed Shanks labels at a short rational seed, with arithmetic controls."""

from math import gcd, isqrt
from sympy import integer_nthroot, isprime, sqrt_mod

from check_shanks_split_seeds import split_seed
from shanks_conductor_relations import shortest_rational_lift
from shanks_witness_fusion import fuse


def branches(a, b, x, y):
    first = a*y-b*x
    second = a*x+b*(x+y)
    return first, second, first+second


def centered_class(residue, modulus, bound):
    start = -bound+(residue+bound) % modulus
    return range(start,bound+1,modulus)


def quadratic_reader(n, a, b, height, multiplier_bound):
    """Find linear factors of a bounded homogeneous quadratic evaluated at (a,b)."""
    if gcd(a,b) != 1 or min(height,multiplier_bound) < 1:
        raise ValueError("expected a primitive seed and positive bounds")
    if abs(a) < abs(b):
        a,b = b,a
    if abs(a) < 2:
        raise ValueError("expected seed height at least two")
    modulus = abs(a)
    inv = pow(b,-1,modulus)
    for magnitude in range(1,multiplier_bound**2 + 1):
        for m in (magnitude,-magnitude):
            for right in centered_class(m*n*inv**2,modulus,4*height**2):
                rest = (m*n-right*b*b)//a
                for middle in centered_class(rest*inv,modulus,8*height**2):
                    left = (rest-middle*b)//a
                    if not left or abs(left) > 4*height**2:
                        continue
                    discr = middle**2-4*left*right
                    if discr < 0:
                        continue
                    root = isqrt(discr)
                    if root*root != discr:
                        continue
                    for sign in (root,-root):
                        content = gcd(2*left,middle+sign)
                        linear = (2*left*a+(middle+sign)*b)//content
                        found = gcd(n,linear)
                        if 1 < found < n:
                            return found
    return None


def family_factor(n):
    """Recognize the fixed two-parameter family by integer roots alone."""
    a = int(integer_nthroot(8*n,3)[0])
    if a % 8 != 1 or a <= 1:
        return None
    product = (a-1)//8
    difference = -((a**3-8*n)//(8*a))
    disc = difference**2+4*product**2
    total = isqrt(disc)
    if total**2 != disc or total < abs(difference):
        return None
    x = isqrt((total-difference)//2)
    y = isqrt((total+difference)//2)
    p = a*y-x
    return p if 1 < p < n and n % p == 0 else None


def shanks_factor(n):
    a = int(integer_nthroot(8*n,3)[0])
    divisor,t = split_seed(n,a)
    if divisor or t is None:
        return divisor
    divisor,pair = shortest_rational_lift(n,t)
    if divisor:
        return divisor
    return fuse(n,[pair],[a],[1],3)[0]


def check():
    fixtures = [(224,419), (7704,12661), (761104,1169707), (246,371)]
    for x,y in fixtures:
        a = 8*x*y+1
        p,q,_ = branches(a,1,x,y)
        n = p*q
        assert isprime(p) and isprime(q) and max(p,q) < 2**64
        assert a == integer_nthroot(8*n,3)[0]
        for found in (family_factor(n),shanks_factor(n),quadratic_reader(n,a,1,y,1)):
            assert found in (p,q)
        print(f"{n.bit_length()} bits: all three readers split {n}")

    checked = 0
    for x in range(2,65,2):
        for y in range((3*x+1)//2,2*x+1):
            if y%2 == 0 or gcd(x,y) != 1:
                continue
            a = 8*x*y+1
            p,q,_ = branches(a,1,x,y)
            n = p*q
            assert p == family_factor(n) and n % shanks_factor(n) == 0
            checked += 1
    print("Two-variable family, exhaustive parameter controls:",checked)

    # The rational label alone supplies the native factor-residue list.
    from hh_residue_lists import factor_from_residue_list
    x,y,n = 246,371,48652964904198163
    d = x*x+x*y+y*y
    assert isprime(d) and d**4 > n
    rhs = n*y*pow(x,-1,d) % d
    roots = sqrt_mod(rhs,d,all_roots=True)
    candidates = set(roots) | {n*pow(r,-1,d)%d for r in roots}
    divisor = factor_from_residue_list(n,d,candidates,modulus_primes=[d])
    assert divisor in (270877613,179612351)
    print(f"Label alone: {len(candidates)} residue candidates modulo {d}, returned {divisor}")

    # Non-integral seeds and all choices of colliding branch pair.
    cases = 0
    for a in range(80,90):
        for b in range(3,9):
            if gcd(a,b) != 1:
                continue
            for x,y in [(2,5),(3,4),(4,7)]:
                vals = branches(a,b,x,y)
                for i,j in [(0,1),(1,2),(0,2)]:
                    p,q = vals[i],vals[j]
                    if gcd(p,q) == 1:
                        n = p*q
                        f = quadratic_reader(n,a,b,max(x,y),1)
                        assert f and 1 < f < n and n % f == 0
                        cases += 1
    print("Homogeneous rational reader, coprime branch-pair controls:",cases)


if __name__ == "__main__":
    check()
