"""Public level-59 Manin-symbol probes on N = t² + 3t + 17."""

import json
import sys
from argparse import ArgumentParser
from math import isqrt
from pathlib import Path

from cubic59 import public_parameters

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.newfamily.modsym import _idx, _rep

LEVEL = 59
SHIFTS = (-2, -1, 0, 1, 2)
# For these small primes the irreducibility test for X³ + 2X + 1 gives
# the odd-prime trace in the two-dimensional characteristic-two S3 representation.
CALIBRATION = (3, 5, 7, 11, 13)


def word(a, n):
    """Bottom rows of oriented convergent edges, reduced in P¹(F_59)."""
    if n == 0:
        return ()
    a %= n
    previous, current = 1, 0
    sign = -1
    rows = []
    while n:
        q, r = divmod(a, n)
        previous, current = current, (q * current + previous) % LEVEL
        rows.append(_idx(sign * current, previous, LEVEL))
        a, n, sign = n, r, -sign
    return tuple(rows)


def mask(indices):
    out = 0
    for i in indices:
        out ^= 1 << i
    return out


class Quotient:
    def __init__(self):
        self.pivots = {}

    def reduce(self, value):
        for i in sorted(self.pivots, reverse=True):
            if value & (1 << i):
                value ^= self.pivots[i]
        return value

    def impose(self, value):
        value = self.reduce(value)
        if value:
            self.pivots[value.bit_length() - 1] = value

    @property
    def dimension(self):
        return LEVEL + 1 - len(self.pivots)


def hecke_endpoint(a, n, p):
    """Sum over p+1 double cosets. Intended for small known primes only."""
    value = mask(word(p * a, n))
    for j in range(p):
        value ^= mask(word(a + j * n, p * n))
    return value


def hecke_edges(p):
    at_zero = hecke_endpoint(0, 1, p)
    # (i,1) is the bottom row of [[1,0],[i,1]], endpoints 1/i and 0.
    images = [at_zero ^ hecke_endpoint(1, i, p) for i in range(LEVEL)]
    return images + [at_zero]


def cubic_trace(p):
    # For a good prime, a 3-cycle has trace one and the two other types have trace zero.
    return int(all((x**3 + 2 * x + 1) % p for x in range(p)))


class Manin59:
    def __init__(self):
        self.full = Quotient()
        for i in range(LEVEL + 1):
            c, d = _rep(i, LEVEL)
            s = _idx(d, -c, LEVEL)
            u, v = _idx(d, -c - d, LEVEL), _idx(-c - d, c, LEVEL)
            self.full.impose(mask((i, s)))
            self.full.impose(mask((i, u, v)))
        self.cm = Quotient()
        self.cm.pivots = dict(self.full.pivots)
        for p in CALIBRATION:
            val = cubic_trace(p)
            for i, image in enumerate(hecke_edges(p)):
                self.cm.impose(image ^ (val << i))

    def signature(self, n):
        t = public_parameters(n)["shanks_parameter"]
        root = isqrt(t)
        offsets = list(SHIFTS) + [root, -root]
        words = [word(t + s, n) for s in offsets]
        return {
            "N": n,
            "t": t,
            "offsets": offsets,
            "path_lengths": [len(v) for v in words],
            "full": [self.full.reduce(mask(v)) for v in words],
            "cm": [self.cm.reduce(mask(v)) for v in words],
        }

    def complete_sum(self, n):
        """Exact parity from all N cusps; the count of Euclidean edges includes the first edge."""
        public_parameters(n)
        images = [self.cm.reduce(1 << i) for i in range(LEVEL + 1)]
        winding = images[LEVEL]
        if not winding:
            raise ArithmeticError("the winding symbol vanished")
        if n.bit_length() > 63:
            total, edges = 0, 0
            for b in range(n):
                path = word(b, n)
                total ^= self.cm.reduce(mask(path))
                edges += len(path)
        else:
            import numpy as np

            # Chunk the exact continued fractions. Each live column is one cusp.
            table = np.zeros((LEVEL, LEVEL), dtype=np.uint64)
            for i in range(LEVEL):
                for j in range(LEVEL):
                    if i or j:
                        table[i, j] = images[_idx(i, j, LEVEL)]
            total, edges = 0, 0
            for start in range(0, n, 65536):
                a = np.arange(start, min(n, start + 65536), dtype=np.int64)
                d = np.full_like(a, n)
                prev, cur = np.ones_like(a), np.zeros_like(a)
                sign = -1
                while len(d):
                    q, rem = np.divmod(a, d)
                    prev, cur = cur, ((q % LEVEL) * cur + prev) % LEVEL
                    total ^= int(
                        np.bitwise_xor.reduce(table[(sign * cur) % LEVEL, prev])
                    )
                    edges += len(d)
                    keep = rem != 0
                    a, d, prev, cur = d[keep], rem[keep], prev[keep], cur[keep]
                    sign = -sign
        if total not in (0, winding):
            raise ValueError("sum outside the promised two outcomes")
        return {
            "N": n,
            "bit": int(total == 0),
            "sum": total,
            "winding": winding,
            "cusps": n,
            "edges": edges,
        }


if __name__ == "__main__":
    p = ArgumentParser(description=__doc__)
    p.add_argument("N", type=int)
    p.add_argument(
        "--complete", action="store_true", help="sum all N cusps to evaluate the bit"
    )
    args = p.parse_args()
    model = Manin59()
    print(
        json.dumps(
            model.complete_sum(args.N) if args.complete else model.signature(args.N),
            indent=2,
        )
    )
