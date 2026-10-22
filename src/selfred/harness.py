"""HARNESS (uses the factorisation). Oracles and map outputs for testing the OPERATIONAL reductions."""
from __future__ import annotations

from math import prod
from typing import List

import sympy

from .rigidity import peval


def rad(M: int) -> int:
    return prod(sympy.factorint(M).keys()) if M > 1 else 1


def phi_mod(N: int, l: int) -> int:
    return int(sympy.totient(N)) % l


def rigid_map_output(p: int, q: int, A: List[int], B: List[int]) -> int:
    """What an exact algebraic self-reduction would output: N2 = A(p) B(q)."""
    return peval(A, p) * peval(B, q)
