"""Track A (d): child selectors along a binary section path.  OPERATIONAL.

Inputs: N, public config, seed.  No factors are used.  Every value computed
modulo N is gcd-checked; nontrivial gcds are returned as discoveries with a
classification (see `Discovery`).

Path state: residue r modulo M = 2^k such that (we hope) S_{r,M}(F_N) != 0.
A selector chooses bit b, moving to r' = r + b*M modulo 2M.

Selectors
  RandomSelector            b uniformly random (cost 0) -- null baseline.
  ExactScalarSelector       evaluates both children at random scalars c in R_N
                            with one computation of (1+X)^N mod (X^{2M} - c)
                            (rank 2M); a nonzero value CERTIFIES that child is
                            nonzero.  Retries `tries` scalars, then evaluation in
                            R_N[Z]/(Z^s - a) (rank 2M*s), then falls back to random.
  BudgetSelector(D)         ExactScalar while 2M <= D, random afterwards: the
                            best one can do in the isolation model with algebras
                            of rank <= D (Theorem A1: for 2M > D the two child
                            classes are indistinguishable by such observables).

`descend()` runs a selector from level 0, optionally finishing by gcd scanning
of r + t*M once N/M is small (`finish_ratio`).
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field
from math import gcd
from typing import Dict, List, Optional

from ..common.metrics import Metrics
from .sections import section_values, section_values_ext


@dataclass
class Discovery:
    divisor: int
    level: int
    source: str          # 'section_value_gcd' | 'finish_scan' | 'nonunit'
    kind: str            # set by the harness: STRUCTURAL | ACCIDENTAL (UNCLASSIFIED here)
    residue: int = -1    # evaluated class residue (section_value_gcd only)
    modulus: int = 0     # evaluated class modulus


@dataclass
class LevelLog:
    level: int
    M: int
    bit: int
    certified: bool          # selector certified the chosen child nonzero
    evals: int               # number of algebra evaluations
    rank: int                # rank of each evaluation algebra
    q_mul_dimsum: int        # cumulative cost counter after this level
    values_gcd_hits: int


class RandomSelector:
    name = "random"

    def choose(self, N, r, M, rng, m, disc):
        return rng.getrandbits(1), False, 0, 0


def _check_values(N, vals, level, disc: List[Discovery], idx):
    M2 = len(vals)
    hits = 0
    for i in idx:
        g = gcd(vals[i], N)
        if 1 < g < N:
            hits += 1
            # A nonzero section value divisible by exactly one prime component
            # happens when the class is hit by only one prime's Lucas support
            # (notes/trackA_analysis.md, Prop. C3); the harness classifies it.
            disc.append(Discovery(g, level, "section_value_gcd", "UNCLASSIFIED", i, M2))
    return hits


class ExactScalarSelector:
    name = "exact_scalar"

    def __init__(self, tries: int = 3, ext_s: int = 2, check_gcd: bool = True):
        self.tries = tries
        self.ext_s = ext_s
        self.check_gcd = check_gcd

    def choose(self, N, r, M, rng, m, disc):
        level = M.bit_length() - 1
        evals = 0
        for _ in range(self.tries):
            c = rng.randrange(1, N)
            vals = section_values(N, 2 * M, c, m)
            evals += 1
            if self.check_gcd:
                _check_values(N, vals, level, disc, (r, r + M))
            if vals[r]:
                return 0, True, evals, 2 * M
            if vals[r + M]:
                return 1, True, evals, 2 * M
        a = rng.randrange(1, N)
        T = section_values_ext(N, 2 * M, self.ext_s, a, m)
        evals += 1
        if any(T[r]):
            return 0, True, evals, 2 * M * self.ext_s
        if any(T[r + M]):
            return 1, True, evals, 2 * M * self.ext_s
        return rng.getrandbits(1), False, evals, 2 * M


class BudgetSelector:
    def __init__(self, D: int, tries: int = 3):
        self.D = D
        self.inner = ExactScalarSelector(tries=tries)
        self.name = f"budget_{D}"

    def choose(self, N, r, M, rng, m, disc):
        if 2 * M <= self.D:
            return self.inner.choose(N, r, M, rng, m, disc)
        return rng.getrandbits(1), False, 0, 0


def descend(N: int, selector, seed: int = 0, max_level: Optional[int] = None,
            finish_ratio: Optional[float] = None, stop_on_discovery: bool = False,
            m: Optional[Metrics] = None):
    """Run `selector` down the binary section tree of F_N.

    finish_ratio: when N // M <= finish_ratio * M, scan i = r + t*M (t < N/M+1)
    computing gcd(i, N) (cost N/M gcds) and stop.  Returns (r, M, logs, discoveries).
    """
    m = m if m is not None else Metrics()
    rng = random.Random(f"descend|{N}|{selector.name}|{seed}")
    max_level = max_level if max_level is not None else N.bit_length() + 1
    r, M, logs, disc = 0, 1, [], []
    for k in range(max_level):
        if finish_ratio is not None and N // M <= finish_ratio * M:
            for t in range(N // M + 1):
                i = r + t * M
                m.inc("gcd")
                g = gcd(i, N)
                if 1 < g < N:
                    disc.append(Discovery(g, k, "finish_scan", "STRUCTURAL"))
                    break
            break
        bit, cert, evals, rank = selector.choose(N, r, M, rng, m, disc)
        hits = sum(1 for d in disc if d.level == k)
        r = r + bit * M
        M *= 2
        logs.append(LevelLog(k + 1, M, bit, cert, evals, rank, m.get("q_mul_dimsum"), hits))
        if stop_on_discovery and disc:
            break
        if M > N:
            break
    return r, M, logs, disc
