"""HARNESS (uses the factors). Brute-force identity checks and the local-data registry for E2.

Never imported by trichotomy.py.
"""
from __future__ import annotations

from math import gcd
from typing import Dict, Iterable, List, Tuple

from .trichotomy import QuasiPoly


def omega_small(n: int) -> int:
    """Number of distinct prime factors of a SMALL positive integer (harness only, trial division)."""
    k, d = 0, 2
    while d * d <= n:
        if n % d == 0:
            k += 1
            while n % d == 0:
                n //= d
        d += 1
    return k + (1 if n > 1 else 0)


class _UF:
    def __init__(self, n: int):
        self.p = list(range(n))

    def find(self, x: int) -> int:
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def count_conjugacy_classes(elements: List[tuple], gens: List[tuple], mul, inv) -> int:
    """Orbits of conjugation by a generating set (= conjugacy classes of the group generated)."""
    idx = {e: i for i, e in enumerate(elements)}
    uf = _UF(len(elements))
    ginv = [inv(g) for g in gens]
    for e, i in idx.items():
        for g, gi in zip(gens, ginv):
            uf.union(i, idx[mul(mul(g, e), gi)])
    return len({uf.find(i) for i in range(len(elements))})


def _m2(N):
    def mul(A, B):
        a, b, c, d = A
        e, f, g, h = B
        return ((a * e + b * g) % N, (a * f + b * h) % N, (c * e + d * g) % N, (c * f + d * h) % N)

    def inv(A):
        a, b, c, d = A
        det = (a * d - b * c) % N
        di = pow(det, -1, N)
        return ((d * di) % N, (-b * di) % N, (-c * di) % N, (a * di) % N)
    return mul, inv


def k_SL2(N: int) -> int:
    mul, inv = _m2(N)
    els = [(a, b, c, d) for a in range(N) for b in range(N) for c in range(N) for d in range(N)
           if (a * d - b * c) % N == 1]
    gens = [(0, N - 1, 1, 0), (1, 1, 0, 1)]
    return count_conjugacy_classes(els, gens, mul, inv)


def k_GL2(N: int) -> int:
    mul, inv = _m2(N)
    els = [(a, b, c, d) for a in range(N) for b in range(N) for c in range(N) for d in range(N)
           if gcd((a * d - b * c) % N, N) == 1]
    units = [u for u in range(1, N) if gcd(u, N) == 1]
    gens = [(0, N - 1, 1, 0), (1, 1, 0, 1)] + [(u, 0, 0, 1) for u in units]
    return count_conjugacy_classes(els, gens, mul, inv)


def k_Aff(N: int) -> int:
    units = [u for u in range(1, N) if gcd(u, N) == 1]
    els = [(a, b) for a in units for b in range(N)]

    def mul(g, h):          # (a,b)(c,d): x -> a(cx+d)+b
        return ((g[0] * h[0]) % N, (g[0] * h[1] + g[1]) % N)

    def inv(g):
        ai = pow(g[0], -1, N)
        return (ai, (-ai * g[1]) % N)
    gens = [(1, 1)] + [(u, 0) for u in units]
    return count_conjugacy_classes(els, gens, mul, inv)


def k_Heis(N: int) -> int:
    els = [(x, y, z) for x in range(N) for y in range(N) for z in range(N)]

    def mul(g, h):          # [[1,x,z],[0,1,y],[0,0,1]]
        return ((g[0] + h[0]) % N, (g[1] + h[1]) % N, (g[2] + h[2] + g[0] * h[1]) % N)

    def inv(g):
        return ((-g[0]) % N, (-g[1]) % N, (-g[2] + g[0] * g[1]) % N)
    return count_conjugacy_classes(els, [(1, 0, 0), (0, 1, 0)], mul, inv)


def gcdsum_SL2(N: int) -> int:
    return sum(3 ** omega_small(gcd((t * t - 4) % N, N) or N) for t in range(N))


def gcdsum_GL2(N: int) -> int:
    units = [u for u in range(1, N) if gcd(u, N) == 1]
    return sum(2 ** omega_small(gcd((t * t - 4 * d) % N, N) or N) for t in range(N) for d in units)


def pillai(N: int) -> int:
    return sum(gcd(g, N) for g in range(N))


def zero_product_count(N: int) -> int:
    return sum(1 for x in range(N) for y in range(N) if (x * y) % N == 0)


def subgroups_ZN2(N: int) -> int:
    pts = [(a, b) for a in range(N) for b in range(N)]
    seen = set()
    for i, g1 in enumerate(pts):
        for g2 in pts[i:]:
            H = frozenset(((x * g1[0] + y * g2[0]) % N, (x * g1[1] + y * g2[1]) % N)
                          for x in range(N) for y in range(N))
            seen.add(H)
    return len(seen)


def aff_degree_multiset(p: int, q: int) -> List[Tuple[int, int]]:
    """Predicted irreducible degrees of Aff(Z/pq) = Aff(F_p) x Aff(F_q): (degree, multiplicity)."""
    return [(1, (p - 1) * (q - 1)), (p - 1, q - 1), (q - 1, p - 1), ((p - 1) * (q - 1), 1)]


# ------------------------------------------------------------------ local-data registry for E2

def chi4_coeffs(a1, a3):
    return {1: a1, 3: a3}


REGISTRY: Dict[str, Tuple[QuasiPoly, str]] = {
    "phi":          (QuasiPoly("phi: r-1", 1, {0: (-1, 1)}), "FACT"),
    "psi":          (QuasiPoly("psi=sigma: r+1", 1, {0: (1, 1)}), "FACT"),
    "pillai":       (QuasiPoly("Pillai: 2r-1", 1, {0: (-1, 2)}), "FACT"),
    "kSL2":         (QuasiPoly("k(SL2): r+4", 1, {0: (4, 1)}), "FACT"),
    "kGL2_J2":      (QuasiPoly("k(GL2)=J2=#K3: r^2-1", 1, {0: (-1, 0, 1)}), "FACT"),
    "kHeis":        (QuasiPoly("k(Heis): r^2+r-1", 1, {0: (-1, 1, 1)}), "FACT"),
    "subgroups":    (QuasiPoly("#subgroups of Z_r^2: r+3", 1, {0: (3, 1)}), "FACT"),
    "unitcircle":   (QuasiPoly("unit circle: r - chi4(r)", 4, chi4_coeffs((-1, 1), (1, 1))), "FACT"),
    "K5":           (QuasiPoly("#K5: r^3-1", 1, {0: (-1, 0, 0, 1)}), "FACT"),
    "phi_sq":       (QuasiPoly("(r-1)^2", 1, {0: (1, -2, 1)}), "FACT"),
    "r2type":       (QuasiPoly("1+chi4(r)", 4, chi4_coeffs((2,), (0,))), "FIBER"),
    "nu3":          (QuasiPoly("1+(-3/r)", 3, {1: (2,), 2: (0,)}), "FIBER"),
    "r_times_r2":   (QuasiPoly("r(1+chi4(r))", 4, chi4_coeffs((0, 2), (0,))), "FIBER"),
    "eps_type":     (QuasiPoly("1 if r=1 mod 4 else 2", 4, chi4_coeffs((1,), (2,))), "FIBER"),
    "agg_chi4_r2":  (QuasiPoly("chi4(r) r^2", 4, chi4_coeffs((0, 0, 1), (0, 0, -1))), "AGG"),
}
