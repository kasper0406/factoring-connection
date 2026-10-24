"""Point counts of elliptic curves over F_l (brute force; small l only). Diagnostic/test data for E1."""
from typing import Dict, Tuple

# Cremona labels -> [a1, a2, a3, a4, a6]
CREMONA: Dict[str, Tuple[int, int, int, int, int]] = {
    "11a1": (0, -1, 1, -10, -20),
    "14a1": (1, 0, 1, 4, -6),
    "15a1": (1, 1, 1, -10, -10),
    "17a1": (1, -1, 1, -1, -14),
    "19a1": (0, 1, 1, -9, -15),
    "21a1": (1, 0, 0, -4, -1),
    "26a1": (1, 0, 1, -5, -8),
    "26b1": (1, -1, 1, -3, 3),
    "33a1": (1, 1, 0, -11, 0),
}


def discriminant(a):
    a1, a2, a3, a4, a6 = a
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    return -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6


def count_points(a, l: int) -> int:
    """#E(F_l) including the point at infinity (l prime of good reduction)."""
    a1, a2, a3, a4, a6 = a
    n = 1
    for x in range(l):
        rhs = (x * x * x + a2 * x * x + a4 * x + a6) % l
        for y in range(l):
            if (y * y + a1 * x * y + a3 * y - rhs) % l == 0:
                n += 1
    return n


def ap(a, l: int) -> int:
    return l + 1 - count_points(a, l)


def count_points_fl2(a, l: int) -> int:
    """#E(F_{l^2}) by brute force; F_{l^2} = F_l[w]/(w^2 - nr) for a non-residue nr (odd l) or w^2 + w + 1 (l = 2)."""
    a1, a2, a3, a4, a6 = a
    if l == 2:
        # w^2 = w + 1 (x^2 + x + 1 irreducible over F_2)
        def mul(x, y):
            (u, v), (s, t) = x, y  # u + v w
            return ((u * s + v * t) % 2, (u * t + v * s + v * t) % 2)
    else:
        nr = next(z for z in range(2, l) if pow(z, (l - 1) // 2, l) == l - 1)

        def mul(x, y):
            (u, v), (s, t) = x, y
            return ((u * s + nr * v * t) % l, (u * t + v * s) % l)

    def add(*xs):
        return (sum(x[0] for x in xs) % l, sum(x[1] for x in xs) % l)

    def c(k):
        return (k % l, 0)

    els = [(u, v) for u in range(l) for v in range(l)]
    n = 1
    for x in els:
        x2 = mul(x, x)
        rhs = add(mul(x2, x), mul(c(a2), x2), mul(c(a4), x), c(a6))
        for y in els:
            lhs = add(mul(y, y), mul(mul(c(a1), x), y), mul(c(a3), y))
            if lhs == rhs:
                n += 1
    return n
