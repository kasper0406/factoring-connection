"""OPERATIONAL: everywhere randomized reduction  M = #E(Z/N)  =>  factor of N, for ANY fixed curve E (cx2 item 1,
GPT-6.1-sol adversarial round, phase 25; re-implemented and re-checked here).

One trial: uniform monic quintic F over Z/N, R = (Z/N)[z]/(F), uniform x in R, delta = f(x),
P = (delta x, delta^2) on Y^2 = X^3 + a delta^2 X + b delta^3 over R; compute [M]P with affine arithmetic over R.
Every non-invertible element met is examined (elimination on the 5x5 multiplication matrix with gcd checks on
every pivot candidate and updated entry, plus a coefficient gcd; on the splitting event equivalent to a norm gcd) with
N is output.  On the event "F = (linear)(irreducible quartic) mod p and F irreducible mod q" (probability -> 1/20),
the F_p component is killed by M with probability ~1/2 (delta a square there) while the F_{q^5} component is never
annihilated (exponent >= q^{5/2} - 1 > M), so a split is forced with probability >= 1/80 - o(1).

Inputs: N, (a, b), M, rng.  No factorisation, no primality test, no harness import.
"""
from __future__ import annotations

import random
from math import gcd
from typing import Dict, List, Optional, Tuple

DEG = 5


class Found(Exception):
    def __init__(self, g: int):
        self.g = g


class Obstruction(Exception):
    pass


class Ring:
    def __init__(self, N: int, F: List[int]):
        self.N = N
        self.F = F  # monic, F[0..4] low coefficients, z^5 = -sum F[i] z^i

    def mul(self, u: List[int], v: List[int]) -> List[int]:
        N = self.N
        prod = [0] * (2 * DEG - 1)
        for i, ui in enumerate(u):
            if ui:
                for j, vj in enumerate(v):
                    prod[i + j] += ui * vj
        for k in range(2 * DEG - 2, DEG - 1, -1):
            c = prod[k] % N
            if c:
                for i in range(DEG):
                    prod[k - DEG + i] -= c * self.F[i]
        return [x % N for x in prod[:DEG]]

    def add(self, u, v):
        return [(x + y) % self.N for x, y in zip(u, v)]

    def sub(self, u, v):
        return [(x - y) % self.N for x, y in zip(u, v)]

    def const(self, c: int) -> List[int]:
        return [c % self.N] + [0] * (DEG - 1)

    def is_zero(self, u) -> bool:
        return all(x == 0 for x in u)

    def _check_entry(self, x: int) -> None:
        g = gcd(x, self.N)
        if 1 < g < self.N:
            raise Found(g)

    def inv(self, w: List[int]) -> List[int]:
        """Inverse of w, or raise Found(factor) / Obstruction (non-unit with no proper gcd visible)."""
        N = self.N
        for x in w:
            self._check_entry(x)
        # multiplication matrix columns: w * z^j
        cols = []
        e = w[:]
        for j in range(DEG):
            cols.append(e)
            e = self.mul(e, [0, 1, 0, 0, 0])
        A = [[cols[j][i] for j in range(DEG)] + [1 if i == 0 else 0] for i in range(DEG)]  # solve A y = e_0
        for c in range(DEG):
            piv = None
            for r in range(c, DEG):
                x = A[r][c] % N
                self._check_entry(x)
                if x and gcd(x, N) == 1:
                    piv = r
                    break
            if piv is None:
                # no unit pivot: w is a zero divisor; try norm / coefficient gcds
                cg = 0
                for x in w:
                    cg = gcd(cg, x)
                g = gcd(cg, N)
                if 1 < g < N:
                    raise Found(g)
                raise Obstruction()
            A[c], A[piv] = A[piv], A[c]
            ip = pow(A[c][c], -1, N)
            A[c] = [x * ip % N for x in A[c]]
            for r in range(DEG):
                if r != c and A[r][c] % N:
                    f = A[r][c]
                    A[r] = [(x - f * y) % N for x, y in zip(A[r], A[c])]
                    for x in A[r]:
                        self._check_entry(x)
        return [A[i][DEG] for i in range(DEG)]

    def is_unit_or_zero(self, w) -> str:
        """'zero', 'unit', or raise Found / Obstruction."""
        if self.is_zero(w):
            return "zero"
        self.inv(w)
        return "unit"


def ec_add(R: Ring, P, Q, A: List[int]):
    if P is None:
        return Q
    if Q is None:
        return P
    (x1, y1), (x2, y2) = P, Q
    dx = R.sub(x2, x1)
    if R.is_zero(dx):
        s = R.add(y1, y2)
        if R.is_zero(s):
            return None
        inv2y = R.inv(s)  # = 2 y1 when y1 = y2 componentwise; raises on mixed cases
        x1sq = R.mul(x1, x1)
        num = R.add(R.mul(R.const(3), x1sq), A)
        lam = R.mul(num, inv2y)
    else:
        lam = R.mul(R.sub(y2, y1), R.inv(dx))
    x3 = R.sub(R.sub(R.mul(lam, lam), x1), x2)
    y3 = R.sub(R.mul(lam, R.sub(x1, x3)), y1)
    return (x3, y3)


def ec_mul(R: Ring, k: int, P, A):
    Rp = None
    for bit in bin(k)[2:]:
        Rp = ec_add(R, Rp, Rp, A)
        if bit == "1":
            Rp = ec_add(R, Rp, P, A)
    return Rp


def ext_trial(N: int, a: int, b: int, M: int, rng: random.Random) -> Tuple[Optional[int], str]:
    F = [rng.randrange(N) for _ in range(DEG)]
    R = Ring(N, F)
    x = [rng.randrange(N) for _ in range(DEG)]
    try:
        delta = R.add(R.add(R.mul(R.mul(x, x), x), R.mul(R.const(a), x)), R.const(b))
        if R.is_zero(delta):
            return None, "delta_zero"
        R.inv(delta)  # delta must be a unit (else factor or obstruction)
        d2 = R.mul(delta, delta)
        A = R.mul(R.const(a), d2)
        P = (R.mul(delta, x), d2)
        ec_mul(R, M, P, A)
    except Found as f:
        return f.g, "split"
    except Obstruction:
        return None, "obstruction"
    return None, "no_split"


def ext_factor(N: int, a: int, b: int, M: int, rng: random.Random, max_trials: int = 400,
               trial_div: int = 101) -> Dict:
    for r in range(2, trial_div):
        if N % r == 0 and 1 < r < N:
            return {"factor": r, "trials": 0, "outcomes": {"trial_division": 1}}
    outcomes: Dict[str, int] = {}
    for t in range(1, max_trials + 1):
        g, out = ext_trial(N, a, b, M, rng)
        outcomes[out] = outcomes.get(out, 0) + 1
        if g is not None:
            assert 1 < g < N and N % g == 0
            return {"factor": g, "trials": t, "outcomes": outcomes}
    return {"factor": None, "trials": max_trials, "outcomes": outcomes}
