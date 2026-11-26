"""Controls for the Cardano symbol and norm-one square probes at discriminant -59."""

import json
import sys
from argparse import ArgumentParser
from collections import Counter
from itertools import islice, product
from pathlib import Path
from random import Random

from class_trace_check import balanced_cases, frobenius_gcd, mixed_class
from cubic59 import CyclicCubic, public_parameters

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.newfamily.diag import cubic_symbol_prime


def local_data(r, t):
    algebra = CyclicCubic(r, t)
    roots = [int(root) for root, multiplicity in algebra.f.roots() if multiplicity == 1]
    assert len(roots) in (0, 3)
    return algebra, roots


def is_local_square(u, local):
    algebra, roots = local
    r = algebra.n
    value = algebra.ring(list(map(int, u)))
    if roots:
        return all(pow(int(value(root)), (r - 1) // 2, r) == 1 for root in roots)
    return value.pow_mod((r**3 - 1) // 2, algebra.f) == 1


def tiny_exact():
    """Enumerate every norm-one unit at the two smallest good test primes used here."""
    rows = []
    for r, c in [(5, 1), (17, 3)]:
        assert (c * c + 59) % r == 0
        t = (c - 3) * pow(2, -1, r) % r
        algebra, roots = local_data(r, t)
        signs, squares = Counter(), 0
        for coords in product(range(r), repeat=3):
            u = algebra.ring(list(coords))
            if algebra.norm(u) == 1:
                squares += is_local_square(u, (algebra, roots))
                pair, _ = algebra.shifted_symbols(u)
                signs[",".join(map(str, pair))] += 1
        order = (r - 1) ** 2 if roots else r * r + r + 1
        assert sum(signs.values()) == order and squares == (
            order // 4 if roots else order
        )
        rows.append(
            {
                "prime": r,
                "root_count": len(roots),
                "order": order,
                "squares": squares,
                "shifted_signs": dict(sorted(signs.items())),
            }
        )
    return rows


def check_tier(start, size, draws):
    labels, signs, square_counts = (
        Counter(),
        [Counter(), Counter()],
        [Counter(), Counter()],
    )
    symbols, examples, factors = Counter(), {}, set()
    quartic_profiles = Counter()
    quartic_examples = {}
    halving = [Counter(), Counter()]
    attempts = 0
    min_bits, max_bits = 10**9, 0
    for p, q, c in islice(balanced_cases(59, start), size):
        n = p * q
        params = public_parameters(n)
        assert params["c"] == c
        t = params["shanks_parameter"]
        algebra = CyclicCubic(n, t)
        locals_ = [local_data(r, t) for r in (p, q)]
        original_degrees = [frobenius_gcd(r, r)[0] for r in (p, q)]
        assert original_degrees == [len(roots) for _, roots in locals_]
        bit = int(mixed_class(p, q, c) != (1, 1, 15))
        assert original_degrees == [0 if bit else 3] * 2
        quartic_roots = []
        for local, roots in locals_:
            quartic = local.ring(params["halving_quartic"])
            count = len(quartic.roots())
            square = is_local_square(local.ring.gen(), (local, roots))
            assert count == (1 if bit else 4 if square else 0)
            if bit:
                prime_half = local.public_half_probe(params["halving_quartic"])
                assert prime_half["half_gcd"] == prime_half["quartic_gcd"] == local.n
            quartic_roots.append(count)
        global_count = quartic_roots[0] * quartic_roots[1]
        quartic_profiles[global_count] += 1
        quartic_examples.setdefault(global_count, {"N": n, "t": t, "bit": bit})
        half = algebra.public_half_probe(params["halving_quartic"])
        g = half["quartic_gcd"]
        assert g in (1, p, q, n) and g % half["half_gcd"] == 0
        if g not in (1, n):
            factors.add(n)
        if g == n:
            assert all(quartic_roots)
        halving[bit]["root" if g == n else "factor" if g > 1 else "miss"] += 1
        labels[bit] += 1
        local_symbols = [cubic_symbol_prime(*params["beta"], r) for r in (p, q)]
        assert [v != 0 for v in local_symbols] == [bool(bit)] * 2
        assert sum(local_symbols) % 3 == params["cubic_jacobi_exponent"] == 0
        symbol_key = ",".join(map(str, local_symbols))
        symbols[symbol_key] += 1
        examples.setdefault(
            symbol_key, {"N": n, "c": c, "factor_control": [p, q], "bit": bit}
        )
        min_bits, max_bits = (
            min(min_bits, n.bit_length()),
            max(max_bits, n.bit_length()),
        )
        z = algebra.ring.gen()
        s1, s2 = algebra.conjugates(z)
        assert s2.compose_mod(algebra.sigma_z, algebra.f) == z
        assert (
            algebra.f.compose_mod(s1, algebra.f) == 0
            and algebra.f.compose_mod(s2, algebra.f) == 0
        )
        rng = Random(n)
        for _ in range(draws):
            attempts += 1
            a = algebra.ring([rng.randrange(n) for _ in range(3)])
            u, divisor = algebra.sample_from(a)
            if divisor is not None:
                if divisor < n:
                    assert divisor in (p, q)
                    factors.add(n)
                continue
            assert algebra.norm(u) == 1
            square = all(is_local_square(u, d) for d in locals_)
            if bit:
                assert square
            square_counts[bit]["square" if square else "nonsquare"] += 1
            pair, ds = algebra.shifted_symbols(u)
            assert all(d in (p, q) for d in ds)
            if ds:
                factors.add(n)
            signs[bit][",".join(map(str, pair))] += 1
    assert sum(labels.values()) == size and all(labels.values())
    return {
        "min_prime": start,
        "modulus_bits": [min_bits, max_bits],
        "inputs": size,
        "attempted_draws": attempts,
        "bit_counts": dict(labels),
        "local_cubic_exponents": dict(symbols),
        "quartic_root_counts": dict(quartic_profiles),
        "quartic_examples": quartic_examples,
        "public_half_by_bit": halving,
        "examples": examples,
        "square_oracle_counts_by_bit": square_counts,
        "shifted_signs_by_bit": signs,
        "inputs_with_gcd_factor": len(factors),
    }


def run():
    assert all((x**4 - 8 * x - 5) % 17 for x in range(17))
    return {
        "tiny_exact": tiny_exact(),
        "tiers": [
            check_tier(start, size, 16)
            for start, size in [(101, 128), (2**18, 64), (2**30, 64), (2**40, 64)]
        ],
    }


if __name__ == "__main__":
    p = ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    data = json.dumps(run(), indent=2) + "\n"
    if args.output:
        args.output.write_text(data)
    else:
        print(data, end="")
