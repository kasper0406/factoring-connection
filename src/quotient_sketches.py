"""Track B quotient sketches (agent sketchBD): operational entry points.

Thin re-export of src/trackBD/sketch.py and src/trackBD/linalg.py (operational code only:
N and public parameters in, factors/failures out; every inversion is gcd-checked).
Analysis: notes/trackBD_analysis.md.  Candidate cards: notes/trackBD_candidates.md.
"""
from src.trackBD.sketch import (all_readouts, binomial_modulus, coef_gcds, content_gcd,  # noqa: F401
                                euclid_readout, lifted_derivative_divided, power_sketch,
                                scalar_F, scalar_G, sketch, sketch_G, try_inverse)
from src.trackBD.linalg import (berlekamp_QN, bm_mod_N, elim_profile, hankel_over_c,  # noqa: F401
                                hankel_proj, mult_matrix)


def sketch_factor_attempt(N: int, f, readouts=("coef", "content", "euclid")):
    """One operational attempt: sketch F_N mod (f, N) and return a proper divisor or None."""
    S = sketch(N, f)
    rd = all_readouts(N, f, S)
    return rd["divisors"][0] if rd["divisors"] else None
