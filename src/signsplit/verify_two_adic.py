"""HARNESS-free check of the user's mod-64 identity (and its 2^k generalisation).

Claim (user, parallel session): with 8 | M, r0 -> r0 + 4 shifts the offset C0 = tau
mod 64 by +-16.  General form checked here: for k >= 3, M = 2^k, any odd N, any
(a, b) with the lock b r^2 == a N (mod 2^k) (r^2 == 1 mod 8 for every odd r, so
the lock lattice is common to the 4 roots {+-r, +-r + 2^{k-1}}):

    tau(r + 2^{k-1}) - tau(r) == +-2^{2k-2} (mod 2^{2k})   if a is odd,
                              in {0, 2^{2k-1}}             if a is even.

Run: .venv/bin/python -m src.signsplit.verify_two_adic
"""
from __future__ import annotations

from collections import Counter

from .identity import tau


def run(nmax: int = 400, kmax: int = 6, amax: int = 9):
    stats = Counter()
    for k in range(3, kmax + 1):
        Mk = 1 << k
        h = Mk >> 1
        for N in range(1, nmax, 2):
            for r in range(1, Mk, 2):
                rinv2 = pow(r * r, -1, Mk)
                for a in range(-amax, amax + 1):
                    b0 = (a * N * rinv2) % Mk
                    for b in (b0, b0 - Mk, b0 + Mk):
                        if (a, b) == (0, 0):
                            continue
                        d = (tau(a, b, (r + h) % Mk, N, Mk) - tau(a, b, r, N, Mk)) % (Mk * Mk)
                        q = 1 << (2 * k - 2)
                        if a % 2:
                            ok = d in (q, 3 * q)
                            stats[("odd_a", k, ok)] += 1
                        else:
                            ok = d in (0, 2 * q)
                            stats[("even_a", k, ok)] += 1
    return stats


if __name__ == "__main__":
    for key, v in sorted(run().items()):
        print(key, v)
