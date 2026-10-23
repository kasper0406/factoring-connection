"""navigate core (agent `navigate`, phase 25).

OPERATIONAL helpers (see only N, a target a, the factor base and a seed):
``primes_upto``, ``gauss_split`` (copied from src/smoothbuild/core.py, rule 7), ``SmoothTester`` (batch-gcd style
smoothness test), ``factor_over``, ``verify_relation`` and the three E1 algorithms ``run_s0``, ``run_g1``, ``run_j2``.

No routine here factors N or uses p, q, phi(N) or a square root of the target. Smoothness is decided by the classical
test  "z is y-smooth  iff  P^(2^e) == 0 (mod z)"  with P the product of the factor base and 2^e >= log2 z; factor
exponents are then read off by trial division over the factor base only.
"""
from __future__ import annotations

import math
import random
from typing import Dict, List, Optional, Sequence, Tuple

import gmpy2


def primes_upto(n: int) -> List[int]:
    """All primes <= n (sieve of Eratosthenes)."""
    if n < 2:
        return []
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i:: i] = bytearray(len(range(i * i, n + 1, i)))
    return [i for i in range(n + 1) if s[i]]


def gauss_split(t: int, N: int) -> Tuple[int, int]:
    """Rational reconstruction (copied from smoothbuild): (u, v) with u = t*v (mod N), |u| < sqrt(N), |v| <~ sqrt(N)."""
    t %= N
    r0, r1 = N, t
    s0, s1 = 0, 1
    while r1 * r1 >= N:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1
        s0, s1 = s1, s0 - qq * s1
    return r1, s1


class SmoothTester:
    """y-smoothness test over a fixed factor base, with an operation counter."""

    def __init__(self, fb: Sequence[int]):
        self.fb = list(fb)
        P = gmpy2.mpz(1)
        for p in self.fb:
            P *= p
        self.P = P
        self.tests = 0

    def is_smooth(self, z: int) -> bool:
        self.tests += 1
        z = abs(int(z))
        if z == 0:
            return False
        if z == 1:
            return True
        r = self.P % z
        e = max(1, (z.bit_length() - 1).bit_length())
        for _ in range(e + 1):
            if r == 0:
                return True
            r = r * r % z
        return r == 0


def factor_over(z: int, fb: Sequence[int]) -> Optional[Dict[int, int]]:
    """Exponents of |z| over fb (None if not fb-smooth). Diagnostic/verification helper (trial division by fb only)."""
    z = abs(int(z))
    if z == 0:
        return None
    out: Dict[int, int] = {}
    for p in fb:
        if z == 1:
            break
        if z % p == 0:
            k = 0
            while z % p == 0:
                z //= p
                k += 1
            out[p] = k
    return out if z == 1 else None


def verify_relation(u: int, v: int, a: int, N: int, fb: Sequence[int], height_factor: float = 3.0) -> bool:
    """Exact check of a NAV output: u = a*v (mod N), |u|,|v| fb-smooth, gcd(v, N) = 1, ln|u|+ln|v| <= 3 ln N."""
    if u == 0 or v == 0:
        return False
    if math.gcd(v, N) != 1:
        return False
    if (u - a * v) % N != 0:
        return False
    if factor_over(u, fb) is None or factor_over(v, fb) is None:
        return False
    return math.log(abs(u)) + math.log(abs(v)) <= height_factor * math.log(N) + 1e-9


def reduced_key(u: int, v: int) -> Tuple[int, int]:
    """Canonical form of the rational u/v (sign on the numerator, lowest terms): identifies identical relations."""
    g = math.gcd(u, v)
    u, v = u // g, v // g
    if v < 0:
        u, v = -u, -v
    return (u, v)


# ---------------------------------------------------------------------------------------------------------------
# E1 algorithms (operational). Each returns a dict with outputs and counters. ``null_cb`` is an optional harness
# callback receiving (|u'|, |v'|) for each split whose pieces were tested (used only for the size-matched null).


def _try_split(t: int, N: int, tester: SmoothTester, budget: int, null_cb=None):
    """Split t, test pieces (respecting the test budget). Returns (u', v', status) with status in
    {'ok', 'fail', 'budget', 'degenerate'}."""
    u1, v1 = gauss_split(t, N)
    if u1 == 0 or v1 == 0 or math.gcd(v1, N) != 1:
        return u1, v1, "degenerate"
    if tester.tests >= budget:
        return u1, v1, "budget"
    if null_cb is not None:
        null_cb(abs(u1), abs(v1))
    if not tester.is_smooth(u1):
        return u1, v1, "fail"
    if tester.tests >= budget:
        return u1, v1, "budget"
    if not tester.is_smooth(v1):
        return u1, v1, "fail"
    return u1, v1, "ok"


def run_s0(N: int, a: int, fb: Sequence[int], seed, budget: int = 400, null_cb=None) -> Dict[str, object]:
    rng = random.Random(f"navigate-s0|{N}|{a}|{seed}")
    tester = SmoothTester(fb)
    outs: List[Tuple[int, int]] = []
    splits = 0
    sizes: List[float] = []
    while tester.tests < budget:
        w = rng.choice(fb) * rng.choice(fb)
        t = a * w % N
        u1, v1, st = _try_split(t, N, tester, budget, null_cb)
        splits += 1
        if st in ("ok", "fail"):
            sizes.append(math.log(abs(u1)) + math.log(abs(v1)))
        if st == "ok":
            outs.append((u1, w * v1))
        elif st == "budget":
            break
    return {"outputs": outs, "tests": tester.tests, "splits": splits, "tested_logsizes": sizes}


def run_g1(N: int, a: int, fb: Sequence[int], seed, budget: int = 400, c: int = 16, null_cb=None,
           cap: float = 1.9) -> Dict[str, object]:
    rng = random.Random(f"navigate-g1|{N}|{a}|{seed}")
    tester = SmoothTester(fb)
    outs: List[Tuple[int, int]] = []
    splits = 0
    lnN, lny = math.log(N), math.log(max(fb))
    t, w = a % N, 1
    sizes: List[float] = []
    while tester.tests < budget:
        if math.log(w) + lny > cap * lnN:
            t, w = a % N, 1
        best = None
        for _ in range(c):
            ell = rng.choice(fb)
            u1, v1 = gauss_split(t * ell % N, N)
            splits += 1
            if u1 == 0 or v1 == 0:
                continue
            sc = math.log(abs(u1)) + math.log(abs(v1))
            if best is None or sc < best[0]:
                best = (sc, ell)
        if best is None:
            continue
        ell = best[1]
        t, w = t * ell % N, w * ell
        u1, v1, st = _try_split(t, N, tester, budget, null_cb)
        if st in ("ok", "fail"):
            sizes.append(math.log(abs(u1)) + math.log(abs(v1)))
        if st == "ok":
            outs.append((u1, w * v1))
        elif st == "budget":
            break
    return {"outputs": outs, "tests": tester.tests, "splits": splits, "tested_logsizes": sizes}


def run_j2(N: int, a: int, fb: Sequence[int], seed, budget: int = 400, c: int = 8, m: int = 16,
           cap: float = 1.9) -> Dict[str, object]:
    rng = random.Random(f"navigate-j2|{N}|{a}|{seed}")
    tester = SmoothTester(fb)
    outs: List[Tuple[int, int]] = []
    splits = 0
    lnN, lny = math.log(N), math.log(max(fb))
    t, w = a % N, 1
    done = False
    while not done:
        if math.log(w) + lny > cap * lnN:
            t, w = a % N, 1
        scores = []
        for _ in range(c):
            ell = rng.choice(fb)
            tl = t * ell % N
            s = 0
            for _ in range(m):
                if tester.tests >= budget:
                    done = True
                    break
                w2 = rng.choice(fb) * rng.choice(fb)
                u1, v1, st = _try_split(tl * w2 % N, N, tester, budget)
                splits += 1
                if st == "ok":
                    outs.append((u1, w * ell * w2 * v1))
                    s += 1
                elif st == "budget":
                    done = True
                    break
            scores.append((s, ell))
            if done:
                break
        if done or not scores:
            break
        top = max(s for s, _ in scores)
        ell = rng.choice([e for s, e in scores if s == top])
        t, w = t * ell % N, w * ell
    return {"outputs": outs, "tests": tester.tests, "splits": splits}
