"""Exact rational-tangent batches for w^2 + D*y^2 = H. Integer arithmetic only."""

from math import isqrt


def cube_root(n):
    lo, hi = 0, 1 << ((n.bit_length() + 2) // 3)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid**3 <= n:
            lo = mid
        else:
            hi = mid - 1
    return lo


def slope(A, C, H, x, K):
    """Rational slope u/v with v<=K and |v*f'(x)-u| <= 9/(8K)."""
    num = A * x
    scale = 8 * K * K * (num + 1)
    den = isqrt(C * (H - A * x * x) * scale * scale)
    a, b = num * scale, den
    p0, p1, q0, q1 = 0, 1, 1, 0
    while b:
        k, r = divmod(a, b)
        p, q = k * p1 + p0, k * q1 + q0
        if q > K:
            break
        p0, p1, q0, q1 = p1, p, q1, q
        a, b = b, r
    return -p1, q1


def graph_size(A, C, H):
    return isqrt(H // (2 * A)), max(1, cube_root(isqrt(C * H) // (4 * A)))


def arc_cost(br, j):
    cost = 0
    for A, C in ((br.delta, 1), (1, br.delta)):
        last, K = graph_size(A, C, 4 * br.n * j)
        cost += last + 1 if K == 1 else 5 * (last // K + 1)
    return cost


def graph_points(A, C, H, work):
    """(x,z)>=0, A*x^2+C*z^2=H, A*x^2<=H/2."""
    last, K = graph_size(A, C, H)
    if K == 1:
        for x in range(last + 1):
            t, rem = divmod(H - A * x * x, C)
            work["direct_rows"] += 1
            if not rem:
                z = isqrt(t)
                if z * z == t:
                    yield x, z
        return
    for x0 in range(0, last + 1, K):
        work["arc_blocks"] += 1
        x1 = min(last, x0 + K - 1)
        u, v = slope(A, C, H, x0, K)
        center = isqrt(v * v * (H - A * x0 * x0) // C) - u * x0
        den = A * v * v + C * u * u
        for k in range(center - 2, center + 3):
            work["line_tests"] += 1
            rad = H * den - A * C * k * k
            if rad < 0:
                continue
            r = isqrt(rad)
            if r * r != rad:
                continue
            for vr in (v * r, -v * r) if r else (0,):
                x, rem = divmod(-C * u * k + vr, den)
                if rem or not x0 <= x <= x1:
                    continue
                z, rem = divmod(u * x + k, v)
                if not rem and z >= 0:
                    assert A * x * x + C * z * z == H
                    yield x, z


def arc_shell(br, j, lo, hi, work):
    H = 4 * br.n * j
    D = br.delta
    points = set(graph_points(D, 1, H, work))  # (y,w)
    points.update((y, w) for w, y in graph_points(1, D, H, work))
    ans = []
    for y0, w0 in points:
        if (w0 - D * y0) % 2:
            continue
        for y in (y0, -y0) if y0 else (0,):
            for w in (w0, -w0) if w0 else (0,):
                L = (br.c * y - br.ell * w) // 2
                if lo <= L <= hi:
                    ans.append(L)
    return sorted(ans)
