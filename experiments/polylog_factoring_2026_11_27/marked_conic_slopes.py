"""Exact fine-residue decisions for mixed ternary Grams on balanced semiprimes."""

from math import gcd, isqrt, prod
from pathlib import Path
import sys

from sympy import jacobi_symbol

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.thmA.arith import gauss_reduce

if __package__:
    from .packed_quaternion_collapse import factored_partners, unpack_candidates
    from .pooled_marked_slopes import pooled_slope_factor
else:
    from packed_quaternion_collapse import factored_partners, unpack_candidates
    from pooled_marked_slopes import pooled_slope_factor


def _dot(a, b):
    return sum(i * j for i, j in zip(a, b))


def _transport(u, x, modulus):
    """An integer symmetric G with G²=D² and G*u=D*(x,1,0)."""
    canonical = (x, 1, 0)
    for sign in (1, -1):
        if tuple(sign * t for t in u) == canonical:
            return 1, lambda v: tuple(sign * t for t in v)
    a = x*x + 1
    choices = [(a - sign * _dot(u, canonical), sign) for sign in (1, -1)]
    denominator, sign = min(choices, key=lambda item: (gcd(modulus, item[0]), item[0]))
    normal = tuple(i - sign*j for i, j in zip(u, canonical))

    def transform(v):
        pairing = _dot(normal, v)
        return tuple(sign * (denominator*i - pairing*j) for i, j in zip(v, normal))

    return denominator, transform


def _slope_line(triple, denominator, factors):
    """Recover the marked projective slope at all available prime-power digits."""
    k_square, cross, j_square = triple
    left = right = 0
    modulus = 1
    for prime, exponent in factors.items():
        power = prime**exponent
        common = gcd(power, k_square, cross, j_square)
        if denominator % common:
            return None
        local = power // common
        if local == 1:
            continue
        kk, kj, jj = (v // common for v in triple)
        if jj % prime:
            u, v = kj * pow(jj, -1, local) % local, 1
        elif kk % prime:
            u, v = 1, kj * pow(kk, -1, local) % local
        else:
            return None
        # Glue a unimodular projective generator, allowing different local charts.
        inverse = pow(modulus, -1, local)
        left += modulus * ((u - left) * inverse % local)
        right += modulus * ((v - right) * inverse % local)
        modulus *= local
    return left, right, modulus


def _slopes(line, scale, bound, work):
    """Search a disk after scaling the second coordinate of the slope lattice."""
    left, right, modulus = line
    g = gcd(left, modulus)
    m = modulus // g
    residue = right * pow(left // g, -1, m) % m if m > 1 else 0
    u, v = gauss_reduce((g, scale*residue), (0, scale*m))
    determinant = scale * modulus
    square = _dot(u, u)
    dot = _dot(u, v)
    last = isqrt(square * bound) // determinant
    for j in range(-last, last + 1):
        if j == 0:
            # A primitive unscaled vector on this line can only be ±u.
            coefficients = (-1, 1) if square <= bound else ()
        else:
            radius = isqrt(square * bound - (j*determinant)**2)
            coefficients = range(-((radius+j*dot)//square), (radius-j*dot)//square+1)
        for i in coefficients:
            work["disk_slots"] += 1
            k, scaled_l = i*u[0] + j*v[0], i*u[1] + j*v[1]
            ell = scaled_l // scale
            if (ell > 0 or (ell == 0 and k > 0)) and gcd(k, ell) == 1:
                yield k, ell


def marked_slope_partners(n, x, u, z, fine_factors, *, method="pooled", stats=None):
    """Return every v of Gram (x²+1,−nx,n²) with v=z mod the supplied fine modulus.

    Promise n=pq, distinct odd primes with q≤2p after ordering. fine_factors
    supplies certified prime powers coprime to 2n(x²+1). Empty means modulus 1.
    "enumerate" lists the possible slopes. "pooled" batches their norm gcds.
    """
    a = x*x + 1
    if (n <= 1 or n % 2 == 0 or x < 0 or len(u) != 3 or len(z) != 3
            or _dot(u, u) != a or jacobi_symbol(-a, n) != -1):
        raise ValueError("invalid mixed Gram or marked ternary vector")
    if method not in ("enumerate", "pooled"):
        raise ValueError("expected pooled or enumerate")
    if any(p < 2 or e < 1 for p, e in fine_factors.items()):
        raise ValueError("expected positive exponents of known primes")
    s = prod(p**e for p, e in fine_factors.items())
    if gcd(s, 2*a*n) != 1:
        raise ValueError("the fine modulus must be coprime to twice the input and first norm")
    work = {} if stats is None else stats
    work.pop("pool", None)
    work.update(disk_slots=0, denominators_tested=0, slope_modulus=0, transport_denominator=0,
                proper_factor=None)
    if (_dot(u, z) + n*x) % s or (_dot(z, z) - n*n) % s:
        return ()
    den, transform = _transport(u, x, s)
    work["transport_denominator"] = den
    lifted_n = den*n
    z0, _, z2 = transform(z)
    t = z0 + lifted_n
    line = _slope_line((2*lifted_n - a*t, z2, t), den, fine_factors)
    if line is None:
        return ()
    work["slope_modulus"] = line[2]
    if method == "pooled":
        divisor = gcd(den, n)
        if divisor in (1, n):
            pool_stats = {}
            divisor = pooled_slope_factor(n, a, den, line[2], line[:2], stats=pool_stats)
            work["pool"] = pool_stats
        if divisor is not None:
            work["proper_factor"] = divisor
            return factored_partners(n, x, u, z, s, divisor)
        # If an integral slope has norm coprime to n, its transformed partner
        # is a multiple of n. Invertibility of den mod n makes it axial here.
        axes = [tuple(sign*n if j == i else 0 for j in range(3))
                for i in range(3) for sign in (-1, 1) if sign*u[i] == -x]
        return tuple(sorted(v for v in axes if all((c-d) % s == 0 for c, d in zip(v, z))))
    # Since the inert prime cannot divide k²+A*j² for a primitive slope,
    # k²+A*j² divides 2*den*r, and r≤sqrt(2n).
    bound = 2*den*isqrt(2*n)
    found = []
    for k, j in _slopes(line, max(1, x), bound, work):
        value = k*k + a*j*j
        if value > bound:
            continue
        work["denominators_tested"] += 1
        if 2*lifted_n % value:
            continue
        scalar = 2*lifted_n // value
        row = scalar*j*j
        numerator = transform((-lifted_n + row, -x*row, scalar*k*j))
        if any(c % (den*den) for c in numerator):
            continue
        partner = tuple(c // (den*den) for c in numerator)
        if all((c - d) % s == 0 for c, d in zip(partner, z)):
            found.append(partner)
    return tuple(sorted(found))


def prime_conductor_points(norm, modulus, residue, conductor_prime, *,
                           method="pooled", rng=None, stats=None):
    """Answer a promised packed scalar query with a supplied single prime conductor."""
    points = set()
    for mark in unpack_candidates(norm, modulus, residue, (conductor_prime,), rng):
        n, x, _, ell, params = mark.packing
        for v in marked_slope_partners(n, x, mark.u, mark.z, {ell: 1}, method=method, stats=stats):
            points.add(tuple(params.P*i - params.Q*j for i, j in zip(mark.u, v)))
    return tuple(sorted(points))
