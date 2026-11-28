"""Parity for a single packed sphere class, using only its three scalar inputs."""

from math import gcd, isqrt

if __package__:
    from .packed_ternary_decision import brute_sphere_parity
else:
    from packed_ternary_decision import brute_sphere_parity


def _unpack(norm, modulus, residue):
    if norm < 0 or modulus <= 0 or len(residue) != 3:
        raise ValueError("expected nonnegative norm, positive modulus and three residues")
    if norm < 2 or modulus % 2 == 0:
        return None
    h = isqrt(norm - 1)
    if h*h != norm - 1:
        return None
    ell = modulus
    while ell % 7 == 0:
        ell //= 7
    scale = modulus // ell
    n, x = divmod(h, scale)
    a = x*x + 1
    if n < 1 or n % 2 == 0 or scale < 8*(x+1)*(n+1) or gcd(a, n) != 1:
        return None
    u = tuple((c + scale//2) % scale - scale//2 for c in residue)
    if sum(c*c for c in u) != a:
        return None
    z = tuple(((c-w) // scale) % ell for c, w in zip(u, residue))
    return n, x, u, z, ell


def _mirror(modulus, residue):
    classes = [min(c % modulus, -c % modulus) for c in residue]
    return 0 in classes or len(set(classes)) < 3


def _axis_parity(n, x, u, z, ell):
    count = 0
    for j in range(3):
        for sign in (-1, 1):
            if sign*u[j] != -x:
                continue
            v = tuple(sign*n if i == j else 0 for i in range(3))
            count += all((a-b) % ell == 0 for a, b in zip(v, z))
    return bool(count % 2)


def packed_parity_shortcut(norm, modulus, residue):
    """Return the exact parity for a recognized mirror-fixed class, else None."""
    data = _unpack(norm, modulus, residue)
    if data is not None and _mirror(modulus, residue):
        return _axis_parity(*data)
    return None


def _conic_count(n, x, u, z, ell, stats):
    a = x*x + 1
    # Leave a nonzero coefficient in the other two positions; minimizing
    # their squared length gives the shortest one-dimensional scan.
    d, j = min((a-c*c, i) for i, c in enumerate(u) if c*c < a)
    k, l = (i for i in range(3) if i != j)
    radius = isqrt(n*n*d)
    center = -n*x*u[j]
    lo, hi = -((radius-center) // a), (center+radius) // a
    count = 0
    for vj in range(lo + (z[j]-lo) % ell, hi+1, ell):
        if stats is not None:
            stats["coordinate_tests"] += 1
        dot = -n*x - u[j]*vj
        square = d*(n*n-vj*vj) - dot*dot
        if square < 0:
            continue
        root = isqrt(square)
        if root*root != square:
            continue
        for wedge in ((root, -root) if root else (0,)):
            vk = u[k]*dot - u[l]*wedge
            vl = u[l]*dot + u[k]*wedge
            if vk % d == 0 and vl % d == 0:
                count += (vk//d-z[k]) % ell == 0 and (vl//d-z[l]) % ell == 0
    return count


def native_ternary_parity(norm, modulus, residue, *, stats=None):
    """Answer one exact parity. Use scalar enumeration for an unrelated query."""
    data = _unpack(norm, modulus, residue)
    if stats is not None:
        stats["coordinate_tests"] = 0
    if data is None:
        if stats is not None:
            stats["method"] = "scalar"
        return brute_sphere_parity(norm, modulus, residue)
    if _mirror(modulus, residue):
        if stats is not None:
            stats["method"] = "mirror"
        return _axis_parity(*data)
    if stats is not None:
        stats["method"] = "conic"
    return bool(_conic_count(*data, stats) % 2)
