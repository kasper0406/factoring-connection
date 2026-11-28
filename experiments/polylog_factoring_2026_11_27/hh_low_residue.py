"""Deterministic balanced factoring with a specified low residue of the smaller prime."""

from bisect import bisect_left
from math import gcd, isqrt

if __package__:
    from .affine_prime_window import affine_divisors, window_radius
    from .sqrt_prime_marks import _DecimalRing
else:
    from affine_prime_window import affine_divisors, window_radius
    from sqrt_prime_marks import _DecimalRing
from src.cycloharvey.calibration import ceil_sqrt_ratio
from src.thmA.arith import gauss_reduce, iroot, primes_upto
from src.trackC.polyback import multipoint_eval, product_tree_linear


def _bounded_order(n, base, bound):
    """Exact multiplicative order if it is at most bound, else None."""
    width = isqrt(bound)+1
    babies = []
    value = 1
    for offset in range(width):
        if offset and value == 1:
            return offset
        babies.append((value, offset))
        value = value*base % n
    # At step h a match gives h*width-offset. These ranges are consecutive.
    babies.sort()
    jump = value
    for h in range(1, (bound+width-1)//width+1):
        i = bisect_left(babies, (value, -1))
        if i < width and babies[i][0] == value:
            order = h*width-babies[i][1]
            return order if order <= bound else None
        value = value*jump % n
    return None


def _prime_powers(value):
    powers = {}
    d = 2
    while d*d <= value:
        if value % d == 0:
            power = 1
            while value % d == 0:
                power *= d
                value //= d
            powers[d] = power
        d += 1 if d == 2 else 2
    if value > 1:
        powers[value] = value
    return powers


def hh_order_base(n, bound):
    """On an odd composite return (factor,None) or (None,base) with order > bound."""
    if bound < 1 or n < 9 or n % 2 == 0:
        raise ValueError("expected a positive order bound and an odd composite")
    if bound < n.bit_length()-1:
        return None, 2
    limit = 4
    while True:
        primes = primes_upto(limit)
        if len(primes)*(len(primes)+1)//2 > bound:
            break
        limit *= 2
    # The generated element is the product of components of coprime prime-power orders.
    components = {}
    combined_order, combined_base = 1, 1
    for base in primes:
        if n % base == 0:
            return base, None
        if pow(base, combined_order, n) == 1:
            continue
        order = _bounded_order(n, base, bound)
        if order is None:
            return None, base
        for prime, power in _prime_powers(order).items():
            divisor = gcd(pow(base, order//prime, n)-1, n)
            if 1 < divisor < n:
                return divisor, None
            old_power, old_component = components.get(prime, (1, 1))
            if power > old_power:
                component = pow(base, order//power, n)
                combined_order = combined_order//old_power*power
                combined_base = combined_base*pow(old_component, -1, n)*component % n
                components[prime] = power, component
        if combined_order > bound:
            return None, combined_base
    # All p|n satisfy p=1 mod combined_order and p<=limit^2, by distinct prime products.
    for candidate in range(1+combined_order, min(isqrt(n), limit*limit)+1, combined_order):
        if n % candidate == 0:
            return candidate, None
    raise ValueError("order construction requires composite input")


def hh_parameters(n, m, buckets):
    """Certified half-width in units of m² for the Harvey–Hittmeir search."""
    return ceil_sqrt_ratio(16*n, (m*buckets)**3)


def hh_giants(n, m, residue, buckets, lam):
    """Yield (a,b,j,exponent) per geometric bucket for the SMALLER prime."""
    g = gcd(residue-1, m)
    ratio = n*pow(residue*residue, -1, m) % m
    m2 = m*m
    reciprocal = n*pow(residue, -1, m2) % m2
    lower, upper = ceil_sqrt_ratio(n, 2), isqrt(n)
    while lower <= upper:
        square = lower*lower
        # Images of (1,ratio),(0,m) in HH Proposition 3.3's scaled lattice.
        (v,w), _ = gauss_reduce((n, buckets*(square*ratio-n)), (0,buckets*square*m))
        a = v//n
        b = (w+buckets*v)//(buckets*square)
        tau = (a*reciprocal+b*residue) % m2
        tau0 = (a*n+b*square)//(lower*m2)
        j = (tau0-lam)*m2 + tau
        yield a, b, j, (a*n+b-j)//g
        lower = (lower*(buckets+1)+buckets-1)//buckets


def _recover_sum(n, a, b, total):
    if a == 0 or b == 0:
        coefficient = a or b
        d, rem = divmod(total, coefficient)
        return d if rem == 0 and 1 < d < n and n % d == 0 else None
    disc = total*total-4*a*b*n
    if disc < 0:
        return None
    root = isqrt(disc)
    if root*root == disc:
        for signed in (root, -root):
            d, rem = divmod(total+signed, 2*b)
            if rem == 0 and 1 < d < n and n % d == 0:
                return d
    return None


def _two_part_split(n, base, annihilator):
    odd = annihilator//(annihilator & -annihilator)
    value = pow(base, odd, n)
    while value != 1:
        d = gcd(value-1, n)
        if 1 < d < n:
            return d
        value = value*value % n
    return None


def hh_collision(n, m, base, g, lam, giants, *, stats=None):
    """Share one normalized baby list among all supplied giants with the same g."""
    work = {} if stats is None else stats
    work.update(status="empty", baby_count=0, giant_count=0)
    d = gcd(base, n)
    if d != 1:
        work["status"] = "base_gcd" if d < n else "short_order"
        return d if d < n else None
    stride = m*m//g
    jump = pow(base, stride, n)
    value = 1
    babies = []
    for i in range(2*lam+1):
        babies.append((value, i))
        value = value*jump % n
        d = gcd(value-1, n)
        if 1 < d < n:
            work.update(status="baby_gcd", baby_count=i+1)
            return d
        if d == n:
            divisor = _two_part_split(n, base, stride*(i+1))
            work.update(status="two_part" if divisor else "short_order", baby_count=i+1)
            return divisor
    work["baby_count"] = len(babies)
    babies.sort()
    remaining = []
    for a, b, j, exponent in giants:
        giant = pow(base, exponent, n)
        work["giant_count"] += 1
        index = bisect_left(babies, (giant, -1))
        if index < len(babies) and babies[index][0] == giant:
            total = j+m*m*babies[index][1]
            divisor = _recover_sum(n, a, b, total)
            if divisor:
                work["status"] = "integer_sum"
                return divisor
        else:
            remaining.append(giant)
    if not remaining:
        return None
    ring = _DecimalRing(n)
    polynomial = product_tree_linear(ring, [(-v) % n for v in remaining])[-1][0]
    points = [v for v,_ in babies]
    for value, evaluation in zip(points, multipoint_eval(ring, polynomial, points)):
        divisor = gcd(evaluation, n)
        if 1 < divisor < n:
            work["status"] = "collision"
            return divisor
        if divisor == n:
            # All full-modulus equalities were removed, so some proper gcd is present.
            for giant in remaining:
                divisor = gcd(giant-value, n)
                if 1 < divisor < n:
                    work["status"] = "saturated"
                    return divisor
            raise ArithmeticError("saturated product without a proper divisor")
    return None


def _affine_class(n, m, residue, stats):
    lower, upper = ceil_sqrt_ratio(n, 2), isqrt(n)
    lo, hi = (lower-residue+m-1)//m, (upper-residue)//m
    radius = max(2, window_radius(n))
    shift = residue*pow(m, -1, n) % n
    stats.update(method="affine", intervals=0)
    for start in range(lo, hi+1, 2*radius+1):
        end = min(start+2*radius, hi)
        width = (end-start+1)//2
        center = start+width
        stats["intervals"] += 1
        for divisor in affine_divisors(n, shift+center, width, lower):
            if 1 < divisor < n:
                return divisor
    return None


def hh_low_factor(n, m, residue, *, stats=None):
    """Split n=pq, p<q<2p, given p mod the dyadic m. An incorrect class can miss."""
    if n < 15 or n % 2 == 0 or m < 2 or m & (m-1) or residue % 2 == 0:
        raise ValueError("expected an odd balanced semiprime, odd residue and dyadic modulus >=2")
    residue %= m
    work = {} if stats is None else stats
    if m**8 > n:
        return _affine_class(n, m, residue, work)
    # The balanced interval occupies about log(sqrt(2))*buckets geometric bins.
    # A factor four balances those bins against 2*lambda+1 baby steps.
    buckets = 4*(iroot(n//m**3, 5)+1)
    lam = hh_parameters(n, m, buckets)
    threshold = m*(2*lam+1)
    work.update(method="HH", buckets=buckets, lam=lam, order_bound=threshold)
    divisor, base = hh_order_base(n, threshold)
    if divisor:
        work["status"] = "order_setup"
        return divisor
    return hh_collision(n, m, base, gcd(residue-1, m), lam,
                        hh_giants(n, m, residue, buckets, lam), stats=work)
