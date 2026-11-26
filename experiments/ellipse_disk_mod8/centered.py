"""Centered projective subdivision and an unmasked, form-changing disk decoder."""

import json
import sys
from argparse import ArgumentParser
from collections import Counter
from math import gcd, isqrt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from genusleak.forms import reduce_form


def norm(f, v):
    a, b, c = f
    x, y = v
    return a * x * x + b * x * y + c * y * y


def restrict(f, u, v):
    a, b, c = f
    x, y = u
    s, t = v
    middle = 2 * a * x * s + b * (x * t + y * s) + 2 * c * y * t
    return norm(f, u), middle, norm(f, v)


def reduced_minima(f, u, v):
    """Two successive minima and a shortest vector in the original coordinates."""
    a, b, c = restrict(f, u, v)
    while True:
        shift = (a - b) // (2 * a)
        if shift:
            v = (v[0] + shift * u[0], v[1] + shift * u[1])
            c += shift * (b + shift * a)
            b += 2 * shift * a
        if a <= c:
            return a, c, u
        a, c = c, a
        u, v = v, u


def basis(chart, m, r):
    if chart == "y":
        return (m, 0), (r, 1)
    return (0, m), (1, r)


ROOT_CHARTS = (("y", 2, 0), ("y", 2, 1), ("x", 2, 0))


def radial_count(n, f, cutoff, work):
    """Signed nonzero count modulo 16. The semiprime promise gives squarefreeness."""
    if cutoff <= 0:
        return 0
    height = n * cutoff
    total = 0
    pending = [("y", 1, 0)]
    while pending:
        chart, m, r = pending.pop()
        work["nodes"] += 1
        work["largest_stride"] = max(work["largest_stride"], m)
        a, c, v = reduced_minima(f, *basis(chart, m, r))
        if a > height:
            work["empty"] += 1
        elif c > height:
            work["single_direction"] += 1
            if gcd(*v) == 1:
                step = n // gcd(n, a)
                total += 2 * (isqrt(height // a) // step)
                work["primitive_rays"] += 1
        elif m == 1:
            pending.extend(ROOT_CHARTS)
        else:
            pending.append((chart, 2 * m, r))
            pending.append((chart, 2 * m, r + m))
    return total % 16


def principal(delta):
    b = delta % 2
    return (1, b, (delta + b) // 4)


def primitive_request(f, cutoff):
    """All chart contents are powers of two, hence invertible modulo odd N."""
    a, b, c = f
    g = gcd(gcd(a, b), c)
    if g & (g - 1):
        raise ArithmeticError("unexpected odd content")
    a, b, c = reduce_form((a // g, b // g, c // g))
    return (a, abs(b), c), cutoff // g


def guided_factor(n, ell, c, delta, oracle):
    """Only unmasked signed-count residues from oracle(n, form, cutoff) are used."""
    f = principal(delta)
    cap = isqrt(delta // 3)

    def mixed(B):
        return (oracle(n, f, B) - 4 * isqrt(B)) % 16

    lo, hi = 0, 1
    while not mixed(hi):
        if hi == cap:
            raise ArithmeticError("no mixed orbit within public cap")
        lo, hi = hi, min(2 * hi, cap)
    while hi - lo > 1:
        mid = (hi + lo) // 2
        if mixed(mid):
            hi = mid
        else:
            lo = mid
    first = hi
    b = delta % 2
    # Representatives of both public antipodal pairs of norm N.
    public = (((c - b * ell) // 2, ell), ((-c - b * ell) // 2, ell))
    chart, m, r = "y", 1, 0
    levels = 0
    while True:
        a, second, v = reduced_minima(f, *basis(chart, m, r))
        if second > n * first:
            if a != n * first or gcd(*v) != 1:
                raise ArithmeticError("selected ray missed the primitive first shell")
            x, y = v
            divisor = gcd((c * y - ell * (2 * x + b * y)) // 2, n)
            if divisor <= 1 or divisor >= n or n % divisor:
                raise ArithmeticError("selected point gave no proper factor")
            return {
                "factors": sorted((divisor, n // divisor)),
                "first": first,
                "levels": levels,
                "final_stride": m,
            }
        children = ROOT_CHARTS if m == 1 else ((chart, 2 * m, r), (chart, 2 * m, r + m))
        for new_chart, stride, residue in children:
            # The lattice condition for k*point gives a single divisibility test on k.
            pub = 0
            for x, y in public:
                offset = x - residue * y if new_chart == "y" else y - residue * x
                step = stride // gcd(stride, offset)
                pub += 2 * (isqrt(first) // step)
            subform = restrict(f, *basis(new_chart, stride, residue))
            request, cutoff = primitive_request(subform, first)
            value = (oracle(n, request, cutoff) - pub) % 16
            if value:
                if value not in (2, 4, 6, 8):
                    raise ArithmeticError("invalid signed first-shell multiplicity")
                chart, m, r = new_chart, stride, residue
                break
        else:
            raise ArithmeticError("no occupied primitive child")
        levels += 1


def native_oracle():
    work = Counter()
    cache = {}

    def query(n, f, cutoff):
        work["oracle_calls"] += 1
        key = (n, f, cutoff)
        if key not in cache:
            work["disk_builds"] += 1
            cache[key] = radial_count(n, f, cutoff, work)
        return cache[key]

    return query, work


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("N", type=int)
    parser.add_argument("--cutoff", type=int)
    args = parser.parse_args()
    # Reuse the repo's N-only public branch (its module imports python-flint).
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ellipse_aggregate"))
    from factor import branch_from_n

    br, shortcut = branch_from_n(args.N)
    if shortcut:
        print(json.dumps({"factors": sorted((shortcut, args.N // shortcut))}))
    else:
        query, stats = native_oracle()
        if args.cutoff is None:
            answer = guided_factor(br.n, br.ell, br.c, br.delta, query)
        elif 0 <= args.cutoff <= br.cap:
            answer = {"signed_mod16": query(br.n, principal(br.delta), args.cutoff)}
        else:
            parser.error("cutoff outside the public search interval")
        print(json.dumps({**answer, "work": dict(stats)}, indent=2))
