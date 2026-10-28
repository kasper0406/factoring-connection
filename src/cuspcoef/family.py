"""PUBLIC DATA (no instance information): the fixed curve family F of the cuspcoef pre-registration.

F: y^2 = x^3 + A x + B, A in [-160, 160] minus {0}, B in [1, 160], 4A^3 + 27B^2 != 0, ordered by (max(|A|, B), A, B),
one curve per j-invariant (quadratic twists change a_N only by (d/N)), the 13 CM j-invariants removed.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from typing import List, Tuple

CM_J = {0, 1728, -3375, 8000, -32768, 54000, 287496, -884736, -12288000, 16581375, -884736000,
        -147197952000, -262537412640768000}


def j_invariant(A: int, B: int) -> Fraction:
    d = 4 * A ** 3 + 27 * B ** 2
    return Fraction(1728 * 4 * A ** 3, d)


def disc(A: int, B: int) -> int:
    return -16 * (4 * A ** 3 + 27 * B ** 2)


@lru_cache(maxsize=1)
def family(M: int = 160) -> Tuple[Tuple[int, int], ...]:
    L = [(A, B) for A in range(-M, M + 1) if A != 0 for B in range(1, M + 1) if 4 * A ** 3 + 27 * B ** 2 != 0]
    L.sort(key=lambda ab: (max(abs(ab[0]), ab[1]), ab[0], ab[1]))
    seen = set()
    out: List[Tuple[int, int]] = []
    for A, B in L:
        j = j_invariant(A, B)
        if j in seen or (j.denominator == 1 and j.numerator in CM_J):
            continue
        seen.add(j)
        out.append((A, B))
    return tuple(out)
