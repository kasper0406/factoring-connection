"""Exact power-coset window, order-sweep and collision controls on selected prime pairs."""

from collections import Counter
from fractions import Fraction
from math import gcd, lcm

from sympy import divisors, isprime, nextprime, primerange, primitive_root

if __package__:
    from .affine_prime_window import prime_window_divisor, window_radius
else:
    from affine_prime_window import prime_window_divisor, window_radius


def centered(a, p):
    return (a + p//2) % p - p//2


def power_image(p, e):
    fibers = Counter(pow(y,e,p) for y in range(1,p))
    assert set(fibers.values()) == {gcd(e,p-1)}
    return set(fibers)


def bounded_order_factor(n, z):
    """Public comparator, with at most bit_length(n) gcds."""
    value = 1
    for t in range(1,n.bit_length()+1):
        value = value*z % n
        g = gcd(value-1,n)
        if 1 < g < n:
            return t, g
    return None


def local_orders(h, p):
    candidates = divisors(len(h))
    return {x: next(t for t in candidates if pow(x,t,p) == 1) for x in h}


def proper(f, g):
    return f*(1-g) + g*(1-f)


def measure(p, q, e, a=1):
    n = p*q
    radius = window_radius(n)
    hp, hq = power_image(p,e), power_image(q,e)
    assert gcd(a,n) == 1 and p < q <= 2*p
    ap = {centered(a*x,p) for x in hp if abs(centered(a*x,p)) <= radius}
    aq = {centered(a*x,q) for x in hq if abs(centered(a*x,q)) <= radius}
    rates = (Fraction(len(ap),len(hp)), Fraction(len(aq),len(hq)))
    window = proper(*rates)
    direct = proper(Fraction(1,len(hp)),Fraction(1,len(hq)))
    op, oq = local_orders(hp,p), local_orders(hq,q)
    wins = Counter()
    for x in hp:
        for y in hq:
            hit = (centered(a*x,p) in ap) != (centered(a*y,q) in aq)
            # The first exponent that can separate two different orders
            # is the smaller order; equal orders never separate.
            sweep = op[x] != oq[y] and min(op[x],oq[y]) <= n.bit_length()
            wins[hit,sweep] += 1
    total = len(hp)*len(hq)
    assert window == Fraction(wins[True,False]+wins[True,True],total)
    return (hp,hq), (ap,aq), (window,direct,
            Fraction(wins[False,True]+wins[True,True],total), Fraction(wins[True,False],total))


def check_crt():
    p,q,e = 257,263,2520
    n, radius = p*q, window_radius(p*q)
    assert radius == 2
    image_pairs = Counter()
    for y in range(1,n):
        if gcd(y,n) == 1:
            z = pow(y,e,n)
            image_pairs[z%p,z%q] += 1
    assert set(image_pairs.values()) == {gcd(e,p-1)*gcd(e,q-1)}
    assert len(image_pairs) == (p-1)*(q-1)//(gcd(e,p-1)*gcd(e,q-1))
    for a in (1,3,5):
        images, _, rates = measure(p,q,e,a)
        if a == 5:
            assert all(a not in image for image in images)
        hit_count = 0
        for (x,y), count in image_pairs.items():
            z = x + p * (((y-x)*pow(p,-1,q)) % q)
            expected = (p if abs(centered(a*x,p)) <= radius else 1)
            expected *= q if abs(centered(a*y,q)) <= radius else 1
            assert prime_window_divisor(n,a*z) == expected
            hit_count += count * (1 < expected < n)
        assert Fraction(hit_count,(p-1)*(q-1)) == rates[0]
    return len(image_pairs)


def check_product_bound():
    cosets = nonempty = 0
    for p in (*primerange(239,500), 587, 593, 8191):
        q = next(iter(primerange(p+1,2*p+1)))
        radius = window_radius(p*q)
        max_product = radius**2
        assert 2*max_product < p
        tau = [0]*(max_product+1)
        for j in range(1,max_product+1):
            for multiple in range(j,max_product+1,j):
                tau[multiple] += 1
        # The h cosets of the unique index-h subgroup are labeled by the
        # discrete log modulo h. These logs serve only the finite census.
        gen = int(primitive_root(p))
        logs = {}
        x = 1
        for j in range(p-1):
            logs[x] = j
            x = gen*x % p
        for index in divisors(p-1):
            buckets = {}
            for t in range(-radius,radius+1):
                if t:
                    buckets.setdefault(logs[t%p]%index, []).append(t)
            cosets += index
            nonempty += len(buckets)
            size = (p-1)//index
            for points in buckets.values():
                products = Counter(x*y for x in points for y in points)
                assert len({z % p for z in products}) == len(products) <= size
                assert all(multiplicity <= 2*tau[abs(z)] for z,multiplicity in products.items())
                assert len(points)**2 <= 2*size*max(tau)
    return cosets, nonempty


def run():
    e = lcm(*range(1,11))
    assert e == 2520
    pairs = check_crt()
    cosets, nonempty = check_product_bound()
    fixtures = (
        (257,263,2, {-2,-1,1,2}, {1,2},
         (Fraction(143,1048), Fraction(161,4192), Fraction(2111,4192), Fraction(1,262))),
        (587,593,32, {-2,1,3}, {-1,1},
         (Fraction(398,10841), Fraction(365,21682), Fraction(657,21682), Fraction(72,10841))),
        (8191,8209,18, {1,2,4,8}, {-1,1},
         (Fraction(233,741), Fraction(125,1482), Fraction(1481,1482), Fraction(0))),
    )
    for p,q,y,ap,aq,expected in fixtures:
        n = p*q
        images, windows, rates = measure(p,q,e)
        assert windows == (ap,aq) and rates == expected
        z = pow(y,e,n)
        stats = {}
        assert gcd(z-1,n) == 1 and prime_window_divisor(n,z,stats=stats) == p
        if window_radius(n) >= 3:
            assert stats["method"] == "lattice"
        if p == 587:
            assert (z%p,z%q) == (p-2,152)
            assert bounded_order_factor(n,z) is None
        if p == 8191:
            assert images[0] == {pow(2,t,p) for t in range(13)}
            assert (len(images[0]),len(images[1])) == (13,114)
            assert bounded_order_factor(n,z) == (13,p)
        print(f"{p}*{q}, X={window_radius(n)}: window {rates[0]}, gcd {rates[1]}, "
              f"order sweep {rates[2]}, window-only vs sweep {rates[3]}")

    # Mersenne primes let the geometric orbit be checked at larger sizes
    # without listing the field. This is a family selected using the prime.
    for ell in (13,17,19,31,61):
        p = 2**ell-1
        assert isprime(p)
        q = int(nextprime(p))
        radius = window_radius(p*q)
        orbit = {pow(2,t,p) for t in range(ell)}
        hits = {centered(x,p) for x in orbit if abs(centered(x,p)) <= radius}
        assert len(orbit) == ell and hits == {2**t for t in range(radius.bit_length())}
        assert ell <= (p*q).bit_length()
    print(f"Mersenne orbit at exponent {ell}: {len(hits)}/{len(orbit)} hits")

    # Both local power images collapse. A public shift can still split, but
    # the powered center is then the same constant on every unit input.
    p,q = 257,263
    n = p*q
    universal = lcm(*range(1,301))
    assert power_image(p,universal) == power_image(q,universal) == {1}
    for a, expected in ((1,n), (13,1), (256,p)):
        assert prime_window_divisor(n,a) == expected
        for y in (2,3,19,41):
            assert pow(y,universal,n) == 1
            assert prime_window_divisor(n,a*pow(y,universal,n)) == expected
    print(f"{pairs} exact CRT image pairs in three cosets; product bound on "
          f"{cosets} cosets ({nonempty} nonempty window intersections)")


if __name__ == "__main__":
    run()
