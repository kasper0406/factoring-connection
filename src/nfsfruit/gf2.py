"""OPERATIONAL: GF(2) left-kernel of a relation matrix, rows as Python int bitmasks."""
from __future__ import annotations
from typing import List, Tuple

def left_kernel(rows: List[int]) -> Tuple[int, List[int]]:
    """rows[i] = bitmask of column parities of relation i.
    Returns (rank, kernel_basis) where each kernel vector is a bitmask over relation
    indices; kernel vector j has a distinct top bit (its 'lead' relation)."""
    piv = {}          # pivot column -> (row_value, tag)
    kern = []
    for i, r in enumerate(rows):
        tag = 1 << i
        while r:
            c = r.bit_length() - 1
            if c in piv:
                pr, pt = piv[c]
                r ^= pr
                tag ^= pt
            else:
                piv[c] = (r, tag)
                break
        if r == 0:
            kern.append(tag)
    return len(piv), kern

def coords(S: int, basis: List[int]) -> int:
    """Express S in a triangular kernel basis (as returned by left_kernel).
    Returns bitmask over basis indices; raises ValueError if S not in span."""
    leads = sorted(((b.bit_length() - 1, j) for j, b in enumerate(basis)), reverse=True)
    out = 0
    for lead, j in leads:
        if (S >> lead) & 1:
            S ^= basis[j]
            out |= 1 << j
    if S:
        raise ValueError("not in span")
    return out

def rank_of(vectors: List[int]) -> int:
    piv = {}
    for r in vectors:
        while r:
            c = r.bit_length() - 1
            if c in piv:
                r ^= piv[c]
            else:
                piv[c] = r
                break
    return len(piv)
