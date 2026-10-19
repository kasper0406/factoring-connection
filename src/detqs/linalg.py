"""OPERATIONAL. Incremental exact Gaussian elimination over F_2 with Python-int
bitsets.  Rows are added one at a time; when a row reduces to zero, the
combination (a bitset over row indices) is a new independent left-kernel vector.
Deterministic (pivot = highest set column bit)."""
from __future__ import annotations

from typing import Dict, Optional, Tuple


class IncrementalF2:
    def __init__(self) -> None:
        self.piv: Dict[int, Tuple[int, int]] = {}  # column -> (vector, combo)
        self.nrows = 0

    def add(self, vec: int) -> Optional[int]:
        """Add a row; return the kernel combo (bitset of row indices) if dependent."""
        i = self.nrows
        self.nrows += 1
        combo = 1 << i
        while vec:
            c = vec.bit_length() - 1
            p = self.piv.get(c)
            if p is None:
                self.piv[c] = (vec, combo)
                return None
            vec ^= p[0]
            combo ^= p[1]
        return combo

    @property
    def rank(self) -> int:
        return len(self.piv)
