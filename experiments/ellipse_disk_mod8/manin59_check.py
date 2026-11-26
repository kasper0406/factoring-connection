"""Independent matrix, complete-sum, and matched endpoint checks."""

import json
from argparse import ArgumentParser
from math import lcm
from pathlib import Path

from flint import fmpz_mod_poly_ctx
from manin59 import (
    CALIBRATION,
    LEVEL,
    SHIFTS,
    Manin59,
    cubic_trace,
    hecke_edges,
    mask,
    word,
)
from sympy import isprime, primerange

from src.newfamily.modsym import _idx, _rep, heilbronn

SMALL = [(127, 251), (373, 677), (3947, 7351)]
# First joint three-outcome collision found along t = 5385 + 947835 k.
COLLISION = [(3947, 7351), (59799213017, 108564958261), (65475548071, 99316393997)]


def apply(images, value):
    result = 0
    while value:
        bit = value & -value
        result ^= images[bit.bit_length() - 1]
        value ^= bit
    return result


def check_matrices(m):
    for p in (2, *CALIBRATION):
        h = heilbronn(p)
        for i, image in enumerate(hecke_edges(p)):
            u, v = _rep(i, LEVEL)
            other = mask(_idx(u * a + v * c, u * b + v * d, LEVEL) for a, b, c, d in h)
            assert m.full.reduce(other ^ image) == 0
    at2 = hecke_edges(2)
    winding = m.cm.reduce(mask(word(0, 1)))
    assert (m.full.dimension, m.cm.dimension, winding) == (11, 4, 1)
    assert all(not m.cm.reduce(apply(at2, v) ^ (1 << i)) for i, v in enumerate(at2))
    assert all(not m.cm.reduce(apply(at2, row)) for row in m.cm.pivots.values())
    held = [p for p in primerange(17, 101) if p != LEVEL]
    for p in held:
        b = cubic_trace(p)
        for i, v in enumerate(hecke_edges(p)):
            assert not m.cm.reduce(v ^ (b << i))
    return {
        "relative_dimension": m.full.dimension,
        "quotient_dimension": m.cm.dimension,
        "winding": winding,
        "heilbronn_columns": (1 + len(CALIBRATION)) * (LEVEL + 1),
        "withheld_primes": held,
    }


def root_count(p, q, t):
    assert p < q < 2 * p and isprime(p) and isprime(q)
    coeff = [t - 5, -8, -2 * t, 0, 1]
    return len(fmpz_mod_poly_ctx(p)(coeff).roots()) * len(
        fmpz_mod_poly_ctx(q)(coeff).roots()
    )


def check_complete(m):
    cases = []
    for p, q in SMALL:
        out = m.complete_sum(p * q)
        t = m.signature(p * q)["t"]
        count = root_count(p, q, t)
        assert out["bit"] == count % 2
        if out["N"] == 31877:
            literal = 0
            for b in range(p * q):
                literal ^= mask(word(b, p * q))
            assert m.cm.reduce(literal) == out["sum"]
        cases.append({**out, "quartic_roots_control": count})
    assert {row["quartic_roots_control"] for row in cases} == {0, 1, 16}
    return cases


def check_collisions(m):
    period = LEVEL * lcm(*(s * s - 3 * s + 17 for s in SHIFTS))
    rows = []
    fixed_word = None
    for p, q in COLLISION:
        n = p * q
        sig = m.signature(n)
        t = sig["t"]
        assert (t - 5385) % period == 0
        fixed = tuple(word(t + s, n) for s in SHIFTS)
        for s, w in zip(SHIFTS, fixed, strict=True):
            assert w == word(t + period + s, (t + period) ** 2 + 3 * (t + period) + 17)
        if fixed_word is not None:
            assert fixed == fixed_word
        fixed_word = fixed
        rows.append(
            {
                "N": n,
                "t": t,
                "progression_index": (t - 5385) // period,
                "control_factors": [p, q],
                "quartic_roots": root_count(p, q, t),
                "cm": sig["cm"],
                "full": sig["full"],
            }
        )
    assert len({tuple(r["cm"]) for r in rows}) == 1
    assert len({tuple(r["full"][: len(SHIFTS)]) for r in rows}) == 1
    assert {r["quartic_roots"] for r in rows} == {0, 1, 16}
    return {"period_for_five_fixed_shifts": period, "inputs": rows}


if __name__ == "__main__":
    p = ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    model = Manin59()
    results = {
        "matrices": check_matrices(model),
        "complete": check_complete(model),
        "collision": check_collisions(model),
    }
    data = json.dumps(results, indent=2) + "\n"
    if args.output:
        args.output.write_text(data)
    else:
        print(data, end="")
