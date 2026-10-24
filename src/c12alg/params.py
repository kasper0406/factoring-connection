"""c12alg: leading-order exponent bookkeeping for design 2L_DW(alpha), 3/8 < alpha < 2/5 (notes/c12alg_results.md §1).

Every quantity Q used in a comparison is represented by the leading term of log Q (or of a small positive power
of D, for precisions) as  c * D^e * lam^f  with lam = log log N, D = log N / lam.  A comparison "A << B" holds at
leading order iff (e_A, f_A) < (e_B, f_B) lexicographically; equal exponents are reported as "borderline".
This is a check of the note's arithmetic, not an independent proof (same status as uprime/params.py).
"""
from __future__ import annotations

from fractions import Fraction as Fr
from typing import Dict, List, Tuple

Term = Tuple[Fr, Fr]          # (D-exponent, lam-exponent)


def zeta0(al: Fr) -> Fr:
    """Height-margin exponent: eta0 * log N = D^zeta0 * lam, zeta0 = 3 alpha - 1/2."""
    return 3 * al - Fr(1, 2)


def comparisons(al: Fr) -> List[Tuple[str, Term, Term]]:
    """(name, smaller, larger): each must satisfy smaller << larger at leading order."""
    one = Fr(1)
    z0 = zeta0(al)
    s = 3 - 7 * al                      # u_2 exponent (parent values over B_2), precision exponent
    return [
        ("u2 -> infinity faster than lam^4 (BDH saving at B2)", (Fr(0), Fr(4)), (s, Fr(0))),
        ("level-1 algebraic size dominated by d log V1", (2 * al, one), (2 - 3 * al, one)),
        ("child bound below level-1 value", (4 * al - 1, one), (2 - 3 * al, one)),
        ("(SQ) automatic: u2 log u2 << log B", (s, one), (al, one)),
        ("N^eta0 beats B2^C", (4 * al - 1, one), (z0, one)),
        ("N^eta0 beats L[alpha]^C", (al, one), (z0, one)),
        ("eta0 << mu (slack lower order at level 2)", (z0, one), (2 * al, one)),
        ("eta0 << theta_F (slack lower order at level 1)", (z0, one), (2 - 3 * al, one)),
        ("log V1 << log M", (1 - al, one), (2 * al, one)),
        ("log V2 << log M (child fibres, FIB)", (4 * al - 1, one), (2 * al, one)),
        ("thin loss Xi (child) << BDH/BV saving at B", (3 * al - 1, Fr(1)), (al, Fr(-2))),
        ("thin loss Xi2 (parent) << BDH saving at B2", (2 - 5 * al, Fr(1)), (s, Fr(-2))),
        ("thin loss Xi2 inside Hildebrand range at B2", (2 - 5 * al, one), (4 * al - 1, one)),
        ("parent rational pair loss Xi_r = D^(1-2a) << BDH saving at B2", (1 - 2 * al, Fr(2)), (s, Fr(-2))),
        ("parent rational pair loss Xi_r inside Hildebrand range at B2", (1 - 2 * al, Fr(2)), (4 * al - 1, one)),
        ("box size-spread loss D^(1-2a) << BV saving at B (child)", (1 - 2 * al, Fr(2)), (al, Fr(-2))),
        ("filter-threshold spread D^(2a-1/2) << BV saving at B (child)", (2 * al - Fr(1, 2), Fr(2)), (al, Fr(-2))),
        ("eta independent of polylog: W/mu1 <= 3 u2 (pointwise ratio; structural, recorded)", (Fr(0), Fr(0)), (Fr(0), Fr(1))),
        ("first-order child error D^-alpha << eta = D^-(3-7alpha)", (-al, Fr(9)), (-s, Fr(-9))),
        ("sqrt(A B) << eta: (3-6a)/2 > 3-7a", (-(3 - 6 * al) / 2, Fr(9)), (-s, Fr(-9))),
        ("recursion: FACTOR_U on parent values is L^o(1)", (Fr(2, 5) * (2 - 3 * al), one), (al, one)),
        ("IR rejection exp(-D^(4a-1)) << gamma2", (Fr(0), one), (4 * al - 1, Fr(0))),
        ("content rejection 2^-d << gamma2", (Fr(0), one), (1 - 2 * al, Fr(0))),
        ("k < (d^2+1)/4 (mu_opt branch)", (al, Fr(0)), (2 - 4 * al, Fr(0))),
        ("level floor inactive: mu << d nu", (2 * al - 1, Fr(0)), (1 - 3 * al, Fr(0))),
        ("parent pair moduli << Psi(X_F)", (1 - al, one), (2 - 3 * al, one)),
        ("u_child -> infinity (Harper range, HT)", (Fr(0), Fr(4)), (al, Fr(0))),
    ]


def margins(al: Fr) -> Dict[str, str]:
    out = {}
    for name, a, b in comparisons(al):
        out[name] = "ok" if a < b else ("borderline" if a == b else "FAIL")
    return out


def all_ok(al: Fr) -> bool:
    return all(v == "ok" for v in margins(al).values())


def sizes(al: Fr, d0: Fr = Fr(1), kappa: Fr = Fr(1, 8), a: Fr = Fr(3, 5), b2: Fr = Fr(2)) -> Dict[str, Tuple]:
    """Leading constants (c, D-exponent) of log-sizes in units of lam, for illustration."""
    return {
        "log Y (special-q)": (1 / (2 * kappa), 1 - al),
        "log M": (1 / d0, 2 * al),
        "log X_r (level-1 rational)": (1 / d0, 2 * al),
        "log X_F (level-1 algebraic)": (d0 / (4 * kappa), 2 - 3 * al),
        "log B": (a, al),
        "log B2": (b2, 4 * al - 1),
        "u2 (level-1 alg over B2)": (d0 / (4 * kappa * b2), 3 - 7 * al),
        "log X_F2 (child algebraic)": (3 / d0 + d0 * b2, 2 * al),
        "sigma1 (L[alpha] box exponent)": (al / (a * d0), 0),
        "sigma2 (child box exponent)": (al / a * (4 / d0 + d0 * b2), 0),
    }


def cost_constant(al: Fr, a: Fr, kappa: Fr, d0: Fr, b2: Fr, C_gamma: Fr) -> float:
    """c = max(2a, a + s + max(a, sigma1 + sigma2), (a + s)/(1 - 2^-alpha)), s = 2 kappa (1 + C_gamma) (uprime §6)."""
    s = 2 * kappa * (1 + C_gamma)
    sig1 = al / (a * d0)
    sig2 = al / a * (4 / d0 + d0 * b2)
    return max(float(2 * a), float(a + s + max(a, sig1 + sig2)), float(a + s) / (1 - 2 ** (-float(al))))


def constraint_A1(al: Fr, a: Fr, d0: Fr) -> bool:
    """(A1) a > sigma_2(q) for q near B (b(q) -> 0): a^2 > 4 alpha / d0.  Makes Delta' < q in every child box."""
    return a * a > 4 * al / d0
