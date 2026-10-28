"""OPERATIONAL: E1 toy -- a fruitless-by-construction QS relation matrix plus h Dixon-descended rows.

Base rows: QS relations x^2 - N = smooth with Jacobi symbol J_N(x) = +1.  For Blum N with
N = 1 (mod 8), nfsfruit Cor J1 gives beta_i = [J_N(x_i) = -1] = 0 for every base row, so every
base-only dependency is trivial (Phi == 0): the worst case for a Lee--Venkatesan-type argument.
J_N is computed from N alone (gmpy2.jacobi), so this filter is operational.

Descent rows: r uniform in [1, N-1] (exact Theorem DD; lemmaM used [2, N-2]), R = r^2 mod N,
then lemmaM's `split_and_descend(R)` (a function of R only).  Hidden root r enters only X_i.
"""
from __future__ import annotations

import random
from math import gcd, isqrt
from typing import List

import gmpy2

from ..lemmaM import descent as dsc


class QSHybrid(dsc.DescentSystem):
    def add_qs_jplus(self, scan: int, limit: int) -> int:
        N, x0 = self.N, self.x0
        added = 0
        for t in range(scan):
            if added >= limit:
                break
            x = x0 + t
            if gmpy2.jacobi(x, N) != 1:
                continue
            Q = x * x - N
            self.counters["qs_tests"] += 1
            f = self._fbfac(Q)
            if f is None:
                continue
            self.X.append(x % N)
            self.E.append(f)
            self.rows.append(self._row(f))
            self.kind.append("Q")
            self.hidden_r.append(None)
            added += 1
        return added

    def add_descent_exact(self, rng: random.Random, max_trials: int) -> bool:
        N = self.N
        for _ in range(max_trials):
            self.counters["r_trials"] += 1
            r = rng.randrange(1, N)          # uniform on [1, N-1]
            g = gcd(r, N)
            if g != 1:
                self.counters["accidental"] += 1
                raise dsc.FactorFound(g)
            out = self.split_and_descend(r * r % N)
            if out is None:
                continue
            W, E = out
            self.X.append(r * W % N)
            self.E.append(E)
            self.rows.append(self._row(E))
            self.kind.append("D")
            self.hidden_r.append(r)
            self.counters["desc_accepted"] += 1
            return True
        return False

    def sub_congruence(self, idx: List[int], S: int):
        """Congruence for a kernel vector S over the sub-list idx of relation indices."""
        T = 0
        for j, i in enumerate(idx):
            if (S >> j) & 1:
                T |= 1 << i
        return self.congruence(T)


def fruitful(X: int, Y: int, N: int) -> bool:
    g = gcd((X - Y) % N, N)
    return 1 < g < N
