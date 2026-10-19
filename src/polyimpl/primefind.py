"""Exact exponent arithmetic for prime construction via large prime factors in short intervals.

Lemma PF (notes/polyimpl_results.md section 3). Input theorem T(beta, theta): every interval
[X, X + X^beta] contains n with P+(n) > X^theta. With a deterministic factoring algorithm of
cost n^alpha, a prime p > z = X^theta is found in time z^{e+o(1)}, where e is the minimum of:
  e_fac   = (beta + alpha) / theta                 (factor every integer of the interval);
  e_strip = max(beta, 1 - theta) / theta, valid if theta > 1/2
            (batch-strip all primes <= 2X^{1-theta} by a remainder tree; the cofactor is P+(n)).
The baseline to beat is 1/2 (Polymath4 / Tao-Croot-Helfgott; Lagarias-Odlyzko).
"""
from __future__ import annotations

from fractions import Fraction as F


def exponents(beta: F, theta: F, alpha: F) -> dict:
    e_fac = (beta + alpha) / theta
    e_strip = (max(beta, 1 - theta) / theta) if theta > F(1, 2) else None
    best = e_fac if e_strip is None else min(e_fac, e_strip)
    return {"e_fac": e_fac, "e_strip": e_strip, "best": best}


def factoring_matters(beta: F, theta: F) -> bool:
    """True iff some alpha >= 0 gives best < 1/2 while e_strip (no factoring) does not.

    best < 1/2 for some alpha >= 0 needs beta/theta < 1/2, i.e. theta > 2 beta. Factoring
    is unnecessary when e_strip < 1/2, i.e. theta > 1/2 and max(beta, 1 - theta) < theta/2.
    """
    any_alpha = beta / theta < F(1, 2)
    strip_ok = theta > F(1, 2) and max(beta, 1 - theta) < theta / 2
    return any_alpha and not strip_ok


# Checked input theorems (beta, theta); epsilon losses suppressed in the limit.
KNOWN = {
    "Baker-Harman-Pintz 2001 (prime in X^0.525)": (F(21, 40), F(1)),
    "Merikoski 2021 (x^{1/2} log^c, P > x^{18/19})": (F(1, 2), F(18, 19)),
    "Jia-Liu 2000 (x^{1/2+eps}, P > x^{25/26})": (F(1, 2), F(25, 26)),
    "R. Li 2026 preprint (x^{1/2+eps}, P > x^{35/36})": (F(1, 2), F(35, 36)),
    "Sylvester-Schur (any beta: P > interval length)": (F(1, 4), F(1, 4)),
}


def table(alphas=(F(0), F(1, 6), F(1, 5))):
    rows = []
    for name, (b, t) in KNOWN.items():
        rows.append((name, b, t, [exponents(b, t, a)["best"] for a in alphas],
                     factoring_matters(b, t)))
    return rows


if __name__ == "__main__":
    for name, b, t, es, fm in table():
        print(f"{name}: beta={b} theta={t} best e at alpha=0,1/6,1/5: "
              f"{[str(e) for e in es]} = {[round(float(e), 4) for e in es]}; "
              f"factoring can matter: {fm}")
