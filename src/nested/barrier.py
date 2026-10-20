"""Exact helper computations for notes/nested_results.md (barrier propositions).  DIAGNOSTIC / analysis only.

Functions here are elementary and exact (Fractions / integers); floats only for display.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb, factorial
from typing import Dict, List

from src.thmA import arith


def is_smooth(x: int, B: int) -> bool:
    """x != 0 is B-smooth (trial division; for tests only)."""
    x = abs(x)
    assert x > 0
    for p in arith.primes_upto(B):
        while x % p == 0:
            x //= p
        if x == 1:
            return True
    return x == 1
