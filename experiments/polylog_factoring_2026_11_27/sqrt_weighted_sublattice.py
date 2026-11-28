"""Sample H <= (Z/n)^d with weight sqrt(|H|), for balanced squarefree semiprime n."""

from fractions import Fraction
from random import SystemRandom

if __package__:
    from .affine_prime_window import prime_window_divisor, window_radius
    from .sublattice_sampler import lattice_from_types, type_counts
else:
    from affine_prime_window import prime_window_divisor, window_radius
    from sublattice_sampler import lattice_from_types, type_counts


WINDOWS = 96


def _sign(value):
    return (value > 0) - (value < 0)


def _quadratic_sign(a, b, p):
    """Sign of a+b*sqrt(p), by squaring only when its two terms have opposite signs."""
    x, y = _sign(a), _sign(b)
    if not x or not y or x == y:
        return x or y
    return x*_sign(a*a - p*b*b)


def _field_sign(coefficients, p, q):
    """Exact sign in the basis (1, sqrt(p), sqrt(q), sqrt(p*q))."""
    a, b, c, d = coefficients
    x, y = _quadratic_sign(a, b, p), _quadratic_sign(c, d, p)
    if not x or not y or x == y:
        return x or y
    # Compare the absolute values of X=a+b√p and √q Y=√q(c+d√p).
    return x*_quadratic_sign(a*a + p*b*b - q*(c*c + p*d*d), 2*(a*b-q*c*d), p)


def _hit_probability(p, q):
    n = p*q
    width = 2*window_radius(n)+1
    one_hit = Fraction(width*(p+q-2*width), n)
    return 1-(1-one_hit)**WINDOWS


def _rank_buckets(p, q, dimension, hit):
    """Cumulative algebraic masses for the conditional draw once a factor is known."""
    middle = (dimension+1)//2
    local = []
    for prime in (p, q):
        gaussian = type_counts(prime, 1, dimension)[0][0]
        local.append([count*prime**(r//2) for r, count in enumerate(gaussian)])
    total, cumulative, buckets = [0]*4, [0]*4, []
    for r, at_p in enumerate(local[0]):
        for s, at_q in enumerate(local[1]):
            coefficient = (r % 2) + 2*(s % 2)
            weight = at_p*at_q
            total[coefficient] += hit.numerator*weight
            if (r, s) != (middle, middle):
                cumulative[coefficient] += hit.denominator*weight
                buckets.append(((r, s), tuple(cumulative)))
    if _field_sign([a-b for a, b in zip(total, cumulative)], p, q) < 0:
        raise ArithmeticError("proper-window probability is too small for the rank law")
    buckets.append(((middle, middle), tuple(total)))
    return buckets


def _sample_ranks(p, q, dimension, hit, rng):
    buckets = _rank_buckets(p, q, dimension, hit)
    total = buckets[-1][1]
    ticket, scale = 0, 1
    while True:
        # The same bit stream is refined if its interval straddles a boundary.
        ticket = 256*ticket + rng.randrange(256)
        scale *= 256
        lo, hi = 0, len(buckets)-1
        while lo < hi:
            mid = (lo+hi)//2
            boundary = buckets[mid][1]
            difference = [scale*x-ticket*y for x, y in zip(boundary, total)]
            if _field_sign(difference, p, q) > 0:
                hi = mid
            else:
                lo = mid+1
        ranks, boundary = buckets[lo]
        difference = [scale*x-(ticket+1)*y for x, y in zip(boundary, total)]
        if _field_sign(difference, p, q) >= 0:
            return ranks


def sample_sqrt_weighted(n, dimension, rng=None, *, stats=None):
    """Return the column HNF lifting a sqrt(|H|)-weighted subspace H.

    The input promises n=p*q with distinct primes p<q<=2p. The dimension is
    positive and odd. The algorithm takes n and dimension alone.
    """
    if n < 6 or dimension < 1 or dimension % 2 == 0:
        raise ValueError("expected a balanced squarefree semiprime and positive odd dimension")
    if rng is None:
        rng = SystemRandom()
    if stats is None:
        stats = {}
    stats.update(windows=0, lattice_windows=0, max_lattice_dimension=0, gcd_offsets=0)
    factor = next((v for v in range(2, 257) if n % v == 0), None)
    small = factor is not None
    if not small:
        for _ in range(WINDOWS):
            work = {}
            mark = prime_window_divisor(n, rng.randrange(n), stats=work)
            stats['windows'] += 1
            stats['lattice_windows'] += int(work['method'] == 'lattice')
            stats['max_lattice_dimension'] = max(stats['max_lattice_dimension'], work['dimension'])
            stats['gcd_offsets'] += work['offsets_tested']
            if 1 < mark < n:
                factor = mark
                break
    middle = (dimension+1)//2
    exceptions = {}
    if factor is not None:
        p, q = sorted((factor, n//factor))
        if p == q or q > 2*p:
            raise ValueError("input violates the balanced squarefree semiprime promise")
        hit = Fraction(1) if small else _hit_probability(p, q)
        ranks = _sample_ranks(p, q, dimension, hit, rng)
        exceptions = {prime: (r,) for prime, r in zip((p, q), ranks) if r != middle}
    stats['factor_found'] = factor is not None
    return lattice_from_types(n, dimension, middle, exceptions, rng)
