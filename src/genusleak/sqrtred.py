"""genusleak square-root reduction (OPERATIONAL reduction; the oracle is supplied by the caller).

Reduction (Prop SQ): sample x = prod f_ell^{e_ell} with e_ell uniform in [0, N^2), ask the oracle for some z with
z^2 = x^2, and read a factor of N from the reduced form of z * x^{-1}, which lies in Cl[2].
"""
from __future__ import annotations

import random
from typing import Callable, Dict, List, Optional

from .forms import Form, ambiguous_split, compose, identity, inverse, is_ambiguous, power, prime_form, split_primes

__all__ = ["sample_class", "sqrt_reduction_query"]


def sample_class(D: int, gens: List[Form], N: int, rng: random.Random) -> Form:
    x = identity(D)
    for g in gens:
        x = compose(x, power(g, rng.randrange(N * N)))
    return x


def sqrt_reduction_query(N: int, D: int, gens: List[Form], oracle: Callable[[Form], Form],
                         rng: random.Random) -> Dict:
    x = sample_class(D, gens, N, rng)
    y = compose(x, x)
    z = oracle(y)
    assert compose(z, z) == y, "oracle returned a non-root"
    t = compose(z, inverse(x))
    assert is_ambiguous(t), "z/x must have order <= 2 (reduced ambiguous form)"
    return {"t": t, "split": ambiguous_split(N, t)}
