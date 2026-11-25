"""Shared public inputs and exact verification for the M2 guide comparison."""
import json
import math
import sys
from pathlib import Path

from sympy import primerange

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ARCHIVE = REPO / "experiments/hatch/hatch_h1c"
sys.path.insert(0, str(ARCHIVE))
import experiment as hatch

K = hatch.K
SIZES = (64, 80, 96)
REPS = 192
GRID = (256, 384, 512, 768, 1024, 1536, 2048, 3072, 4096, 6144, 8192)
LOW = 2 * K * K
HIGH = 2 * LOW


def exact_factors(value, primes):
    original = value
    factors = {}
    for p in primes:
        if value == 1:
            break
        if value % p:
            continue
        e = 0
        while value % p == 0:
            value //= p
            e += 1
        factors[p] = e
    assert value == 1 and math.prod(p**e for p, e in factors.items()) == original
    return factors


def matched_divisor(factors, excluded):
    divisors = [1]
    for p, e in factors.items():
        if excluded % p == 0:
            continue
        # A prime dividing a primitive norm and prime to 2aD is split.
        new = []
        for d in divisors:
            for _ in range(e):
                d *= p
                if d > HIGH:
                    break
                if d > LOW:
                    return d
                new.append(d)
        divisors.extend(new)
    return None


def exclusion(form):
    a, b, c = form
    return abs(2 * a * (b*b - 4*a*c))


def truth(form, primes):
    smooth_points, *_ = hatch.full_lattice_sieve(form, primes)
    excluded = exclusion(form)
    matches = {}
    for point in smooth_points:
        factors = exact_factors(hatch.norm(form, *point), primes)
        if matched_divisor(factors, excluded):
            matches[point] = max(factors)
    return len(smooth_points), matches


def archived():
    for line in (ARCHIVE / "runs.jsonl").read_text().splitlines():
        row = json.loads(line)
        if row["bits"] in SIZES:
            yield row["bits"], row["rep"], row["N"], tuple(row["m2"]["form"])


def calibrate():
    totals = {bits: [0] * len(GRID) for bits in SIZES}
    primes = list(primerange(max(GRID) + 1))
    for bits, _, _, form in archived():
        _, matches = truth(form, primes)
        for i, bound in enumerate(GRID):
            totals[bits][i] += sum(p <= bound for p in matches.values())
    bounds = {bits: min(zip(GRID, counts), key=lambda bc: (abs(6 - bc[1]), bc[0]))[0]
              for bits, counts in totals.items()}
    result = dict(grid=GRID, totals=totals, bounds=bounds)
    (HERE / "calibration.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    calibrate()
