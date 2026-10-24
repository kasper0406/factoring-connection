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
