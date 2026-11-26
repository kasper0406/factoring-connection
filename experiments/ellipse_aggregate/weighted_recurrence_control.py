"""Reverse-reduction calibration. This explicitly takes the two prime factors."""

from math import isqrt

from weighted_recurrence import fourths, squares


def factored_observables(br, p, q):
    """Evaluate the full cap including both principal and mixed contributions."""
    n = br.n
    if n != p * q or p == q:
        raise ValueError("expected the two distinct input prime factors")
    b0 = br.delta % 2
    # x/y in the (x,y) integral coordinates, w=2x+b0*y.
    xp = (br.c - b0 * br.ell) * pow(2 * br.ell, -1, p) % p
    xq = (-br.c - b0 * br.ell) * pow(2 * br.ell, -1, q) % q
    s = xp + p * ((xq - xp) * pow(p, -1, q) % q)
    a = n
    b = 2 * s + b0
    d, rem = divmod(s * s + b0 * s + (br.delta + b0) // 4, n)
    assert rem == 0
    yu, yv = 0, 1
    reductions = 0
    while True:
        reductions += 1
        shift = (a - b) // (2 * a)
        yv += shift * yu
        d += shift * (b + shift * a)
        b += 2 * shift * a
        if a <= d:
            break
        a, d = d, a
        yu, yv = yv, yu
    assert abs(b) <= a <= d and 4 * a * d - b * b == br.delta

    k = isqrt(br.cap)
    trace = br.c * k * (k + 1) // 2
    raw = [k, squares(k), fourths(k), br.ell**2 * squares(k), br.ell**2 * fourths(k)]
    m = isqrt(br.cap // a)
    trace += isqrt(4 * n * a - br.delta * yu**2) * m * (m + 1) // 2
    ray = (
        m,
        a * squares(m),
        a * a * fourths(m),
        yu**2 * squares(m),
        a * yu**2 * fourths(m),
    )
    raw = [r + t for r, t in zip(raw, ray, strict=True)]
    for u in (-1, 0, 1):
        j = a * u * u + b * u + d
        if j <= br.cap:
            y2 = (u * yu + yv) ** 2
            trace += isqrt(4 * n * j - br.delta * y2)
            extra = (1, j, j * j, y2, j * y2)
            raw = [r + t for r, t in zip(raw, extra, strict=True)]
    return tuple(r % n for r in raw), trace % n, reductions
