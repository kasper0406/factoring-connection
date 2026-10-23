"""E3 (HARNESS): first square dependency among uniform random integers, all primes allowed (CGPT 2012 check).

Pre-registered in notes/smoothbuild_prereg.md (E3). Run:
    .venv/bin/python -m src.smoothbuild.e3_cgpt
Writes experiments/raw/smoothbuild_E3.jsonl.
"""
from __future__ import annotations

import json
import os
import random
import sys
import time
from typing import Dict, List, Optional

import flint

OUT = "experiments/raw/smoothbuild_E3.jsonl"
SIZES = [32, 40, 48]
SEEDS = list(range(12))
PI_LIMIT_BITS = 26
FHAT_SAMPLES = 200_000


def pi_table(limit: int, points: List[int]) -> Dict[int, int]:
    """Exact pi(y) for each y in points (y <= limit), by a bytearray sieve."""
    s = bytearray([1]) * (limit + 1)
    s[0] = s[1] = 0
    i = 2
    while i * i <= limit:
        if s[i]:
            s[i * i:: i] = bytearray(len(range(i * i, limit + 1, i)))
        i += 1
    out = {}
    for y in points:
        out[y] = s[: y + 1].count(1)
    return out


def factor_support(z: int):
    """Odd-exponent prime support (frozenset) and largest prime factor of z >= 1 (HARNESS: flint factor)."""
    if z == 1:
        return frozenset(), 1
    fs = flint.fmpz(z).factor()
    odd = frozenset(int(p) for p, e in fs if e % 2)
    return odd, max(int(p) for p, _e in fs)


class Elim:
    """Incremental GF(2) elimination on sets of primes, pivot = largest prime."""

    def __init__(self, track: bool = False) -> None:
        self.piv: Dict[int, tuple] = {}
        self.track = track

    def add(self, row: frozenset, idx: int) -> Optional[int]:
        """Insert; return the combination bitmask if the row reduces to zero, else None."""
        r = set(row)
        comb = (1 << idx) if self.track else 0
        while r:
            m = max(r)
            if m in self.piv:
                pr, pc = self.piv[m]
                r ^= pr
                if self.track:
                    comb ^= pc
            else:
                self.piv[m] = (frozenset(r), comb)
                return None
        return comb if self.track else 1


def main(sizes=SIZES, seeds=SEEDS) -> None:
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fh = open(OUT, "a")
    for bits in sizes:
        X = 1 << bits
        t0 = time.time()
        grid = sorted({int(2 ** (k / 2)) for k in range(8, 2 * bits - 4 + 1)})
        grid_pi = [y for y in grid if y <= (1 << PI_LIMIT_BITS)]
        pis = pi_table(1 << PI_LIMIT_BITS, grid_pi)
        rng = random.Random(f"smoothbuild-E3-fhat|{bits}|1000")
        cnt = {y: 0 for y in grid}
        for _ in range(FHAT_SAMPLES):
            _odd, mx = factor_support(rng.randint(1, X))
            for y in grid:
                if mx <= y:
                    cnt[y] += 1
        J = {y: pis[y] * FHAT_SAMPLES / cnt[y] for y in grid_pi if cnt[y] > 0}
        ystar = min(J, key=J.get)
        J0 = J[ystar]
        print(bits, "J0hat", round(J0), "y*", ystar, round(time.time() - t0, 1), flush=True)
        cap = int(6 * J0)
        for seed in seeds:
            rng = random.Random(f"smoothbuild-E3|{bits}|{seed}")
            full = Elim(track=True)
            per_y = {y: Elim() for y in grid}
            Ty: Dict[int, Optional[int]] = {y: None for y in grid}
            T = None
            dep_max_prime = None
            rows: List[int] = []
            for i in range(cap):
                z = rng.randint(1, X)
                odd, mx = factor_support(z)
                rows.append(mx)
                if T is None:
                    c = full.add(odd, i)
                    if c is not None:
                        T = i + 1
                        members = [j for j in range(i + 1) if c >> j & 1]
                        dep_max_prime = max(rows[j] for j in members)
                        dep_size = len(members)
                for y in grid:
                    if Ty[y] is None and mx <= y:
                        if per_y[y].add(odd, i) is not None:
                            Ty[y] = i + 1
                if T is not None and all(Ty[y] is not None for y in grid if y <= ystar * 64):
                    break
            fin = [v for v in Ty.values() if v is not None]
            Tbest = min(fin) if fin else None
            ybest = min((y for y in grid if Ty[y] == Tbest), default=None) if Tbest else None
            row = {"exp": "E3", "bits": bits, "seed": seed, "T": T, "J0hat": J0, "ystar": ystar,
                   "T_over_J0": (T / J0) if T else None, "Tbest": Tbest, "ybest": ybest,
                   "Tbest_over_T": (Tbest / T) if (T and Tbest) else None,
                   "dep_max_prime": dep_max_prime, "dep_size": dep_size if T else None,
                   "n_drawn": len(rows), "cap": cap,
                   "Ty": {str(y): v for y, v in Ty.items() if y <= ystar * 1024}}
            fh.write(json.dumps(row) + "\n")
            fh.flush()
            print(bits, seed, "T", T, "T/J0", round(T / J0, 3) if T else None, "Tbest", Tbest, ybest,
                  "depmax", dep_max_prime, round(time.time() - t0, 1), flush=True)
    fh.close()


if __name__ == "__main__":
    args = [int(a) for a in sys.argv[1:]]
    main(sizes=args or SIZES)
