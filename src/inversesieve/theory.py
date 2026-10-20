"""Exact exponent bookkeeping for the inverse-sieve route (fractions.Fraction; no floats in the ledger).

Ledger (dchworst Thm W, coverage38 Thm U / U2): a single band at lambda = J^theta costs
  X = N^{(1+theta)/(6+theta)},  Y = P^{2/(6+theta)},  J ~ Y,  lambda ~ Y^theta = N^{theta/(6+theta)}.
Coppersmith container: <= P^{1/2} = N^{1/4} integer candidates, polynomial time.
"""
from __future__ import annotations

from fractions import Fraction as Fr


def exponent(theta: Fr) -> Fr:
    """log_N X for a single-band design whose coverage holds at lambda = J^theta (U2)."""
    return (1 + theta) / (6 + theta)


def ledger(theta: Fr = Fr(1, 5)) -> dict:
    """All exponents as powers of P = N^{1/2} (exact)."""
    x_N = exponent(theta)
    y_P = Fr(2) / (6 + theta)                    # Y = P^{y}
    lam_P = theta * y_P                          # lambda = J^theta ~ Y^theta
    X_P = 2 * x_N                                # X = P^{2x}
    mass_frac_P = X_P + Fr(1, 2) - 1             # (X * P^{1/2}) / P  = P^{...}: container mass fraction
    h2_P = lam_P                                 # H^2 = lambda Y / J ~ lambda (up to log)
    return {"theta": theta, "X_N": x_N, "X_P": X_P, "Y_P": y_P, "lambda_P": lam_P,
            "container_mass_fraction_P": mass_frac_P, "H2_P": h2_P,
            "band_primes_per_coppersmith_modulus": (1 / y_P),   # log P / log Y: <= this many band primes in M <= P
            }


def poisson_uncovered_pairs_log(theta: Fr, n: float) -> float:
    """log of the Poisson-model expected number of uncovered pairs |S_n| e^{-lambda} (float, report only)."""
    import math
    P = 2.0 ** (n / 2)
    lam = P ** float(ledger(theta)["lambda_P"])
    return 2 * math.log(P) - 2 * math.log(n) - lam


def moment_order_needed(theta: Fr = Fr(1, 5)) -> Fr:
    """Smallest k such that a Poisson-quality k-th moment bound, exceptional density lambda^{-k/2},
    reaches the container mass fraction P^{-(1/2 - X_P)}: k >= 2*(1/2 - X_P)/lambda_P."""
    L = ledger(theta)
    need = -L["container_mass_fraction_P"]       # positive exponent
    return 2 * need / L["lambda_P"]


def local_container_bound_exponent(m_P: Fr, delta_P: Fr) -> Fr:
    """Theorem L: a container rule seeing only N mod M (M = P^{m}) and N to precision Delta = P^{d}
    has average container MASS >= c * min(P, Delta/M) / n^2.  Returns log_P of min(P, Delta/M)."""
    return min(Fr(1), delta_P - m_P)


def max_forgetful_modulus_exponent(theta: Fr = Fr(1, 5)) -> Fr:
    """With full size knowledge to precision Delta = P^2 (only n known) the mass bound P^{2-m} beats the
    container budget X P^{1/2} iff m < 2 - X_P - 1/2.  Returns that threshold for m (as a power of P)."""
    L = ledger(theta)
    return 2 - L["X_P"] - Fr(1, 2)
