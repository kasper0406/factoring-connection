"""Direct tests for a public rational plane on a mixed packed ternary sphere."""

from math import gcd, isqrt

if __package__:
    from .packed_point_rotations import _cross
    from .packed_ternary_decision import _packing_parameters
else:
    from packed_point_rotations import _cross
    from packed_ternary_decision import _packing_parameters


def _dot(u, v):
    return sum(a*b for a, b in zip(u, v))


def _validate(n, x, u, normal):
    if len(normal) != 3 or not any(normal) or len(u) != 3 or _dot(u, u) != x*x + 1:
        raise ValueError("expected a ternary first vector of norm x²+1 and nonzero normal")
    return _packing_parameters(n, x)[0]


def homogeneous_plane_factor(n, x, u, normal):
    """A proper split is guaranteed if this plane contains a productive packed point."""
    scale = _validate(n, x, u, normal)
    content = gcd(*normal)
    w = tuple(a // content for a in normal)
    for candidate in (scale, _dot(u, w), _dot(w, w)):
        divisor = gcd(candidate, n)
        if 1 < divisor < n:
            return divisor
    return None


def plane_partners(n, x, u, w, height=0):
    """Return all v with w·(u−Qv)=height, or None when the entire conic qualifies."""
    scale = _validate(n, x, u, w)
    a, b, c = x*x + 1, -n*x, n*n
    e, rem = divmod(_dot(w, u) - height, scale)
    if rem:
        return ()
    axis = _cross(u, w)
    det = _dot(axis, axis)
    pairing, norm = _dot(u, w), _dot(w, w)
    if not det:
        return None if a*e == pairing*b else ()
    discriminant = det*c - norm*b*b + 2*b*e*pairing - a*e*e
    if discriminant < 0 or isqrt(discriminant)**2 != discriminant:
        return ()
    first, second = b*norm - e*pairing, a*e - b*pairing
    partners = []
    for t in sorted({isqrt(discriminant), -isqrt(discriminant)}):
        vector = [first*i + second*j + t*k for i, j, k in zip(u, w, axis)]
        if all(entry % det == 0 for entry in vector):
            partners.append(tuple(entry // det for entry in vector))
    return tuple(partners)
