"""The cadobench test numbers: RSA-style balanced semiprimes from the shared seeded generator.

`balanced_semiprime(bits, seed)` (src/common/instances.py) with `bits` chosen so that every N of
that bit length has exactly `digits` decimal digits; seed 0 for every size. The factors are kept
here for the harness check only (factors_ok); nothing passes them to CADO.
"""
from __future__ import annotations

import os
import sys
from typing import Dict, List, Tuple

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from common.instances import balanced_semiprime  # noqa: E402

# digits -> bits; 2^(bits-1) >= 10^(digits-1) and 2^bits < 10^digits.
BITS = {60: 198, 70: 231, 80: 264, 90: 297, 100: 331, 110: 364}
SEED = 0


def _make() -> Dict[int, Tuple[int, int, int, int, int]]:
    out = {}
    for d, b in BITS.items():
        inst = balanced_semiprime(b, SEED)
        assert len(str(inst.N)) == d, (d, b)
        p, q = inst.factors
        out[d] = (inst.N, p, q, SEED, b)
    return out


NUMBERS = _make()


def factors_ok(digits: int, factors: List[int]) -> bool:
    """Harness check: CADO's output is exactly the hidden factor pair."""
    N, p, q, _, _ = NUMBERS[digits]
    return sorted(factors) == [p, q]
