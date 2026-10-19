"""OPERATIONAL GF(2) linear algebra on Python-int bitsets (bit i = coordinate i).

Echelon form keyed by the HIGHEST set bit.  If the "eliminate-first" coordinates are placed at
the high bits, the echelon vectors whose pivot lies below a threshold t form a basis of
R ∩ F_2^{[0,t)} (the top bit of a combination of vectors with distinct pivots is the largest
pivot used, so nothing cancels at the top).
"""
from __future__ import annotations

from typing import Dict, Iterable, List


class Echelon:
    def __init__(self, vecs: Iterable[int] = ()):
        self.piv: Dict[int, int] = {}
        for v in vecs:
            self.add(v)

    def reduce(self, v: int) -> int:
        piv = self.piv
        while v:
            c = v.bit_length() - 1
            w = piv.get(c)
            if w is None:
                return v
            v ^= w
        return 0

    def add(self, v: int) -> bool:
        v = self.reduce(v)
        if v:
            self.piv[v.bit_length() - 1] = v
            return True
        return False

    @property
    def rank(self) -> int:
        return len(self.piv)

    def contains(self, v: int) -> bool:
        return self.reduce(v) == 0

    def basis(self) -> List[int]:
        return list(self.piv.values())

    def basis_below(self, t: int) -> List[int]:
        """Basis of span ∩ F_2^{[0, t)}."""
        return [v for c, v in self.piv.items() if c < t]


def rank(vecs: Iterable[int]) -> int:
    return Echelon(vecs).rank


def dot(a: int, b: int) -> int:
    return bin(a & b).count("1") & 1


def mask_of(coords: Iterable[int]) -> int:
    m = 0
    for c in coords:
        m |= 1 << c
    return m


def dim_in_coords(basis: List[int], mask: int) -> int:
    """dim(span(basis) ∩ F_2^{mask}) for a basis (linearly independent list)."""
    off = ~mask
    return len(basis) - rank(v & off for v in basis)


def restrict_rank(funcs: List[int], mask: int) -> int:
    """rank of the functionals restricted to the coordinates in mask."""
    return rank(f & mask for f in funcs)


def annihilated_dim(funcs: List[int], basis: List[int]) -> int:
    """dim{ sum c_i f_i : <., b> = 0 for all b in basis } for independent funcs f_i:
    = len(funcs) - rank of the |funcs| x |basis| pairing matrix."""
    rows = []
    for f in funcs:
        r = 0
        for j, b in enumerate(basis):
            if dot(f, b):
                r |= 1 << j
        rows.append(r)
    return len(funcs) - rank(rows)


def all_annihilate(funcs: List[int], basis: List[int]) -> bool:
    return all(dot(f, b) == 0 for f in funcs for b in basis)
