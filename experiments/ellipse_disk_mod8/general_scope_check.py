"""Cross-check the general disk-to-factoring reduction with independent arithmetic."""

import json
from collections import Counter
from math import gcd, isqrt, prod
from pathlib import Path
from random import Random

from general_scope import Request, chart, enumerated_oracle, factor, rows, split
from sympy import factorint, nextprime


def signed_roots(factors):
    roots, modulus = [0], 1
    for p, e in factors.items():
        q = p**e
        inverse = pow(modulus, -1, q)
        roots = [
            z + modulus * (((sign - z) * inverse) % q)
            for z in roots
            for sign in (-1, 1)
        ]
        modulus *= q
    return roots


def span(lo, hi, residue, stride):
    return (hi - residue) // stride - (lo - 1 - residue) // stride


def crt_oracle(n, req, roots, axis_stride):
    a, b, d = req.form
    disc = b * b - 4 * a * d
    assert gcd(a, b, d) == 1 and gcd(req.root, n) == 1
    assert (req.root**2 - disc) % n == 0
    assert req.cutoff <= isqrt(-disc // 3)
    inverse = pow(2 * a, -1, n)
    total = 0
    for v, lo, hi in rows(n, req):
        assert abs(v) <= 1
        if v == 0:
            total += span(lo, hi, 0, axis_stride) - 1
        else:
            total += sum(
                span(lo, hi, ((req.root * z - b * v) * inverse) % n, n) for z in roots
            )
    return total % 4


def checked_split(factors, direct, stats, controls):
    n = prod(p**e for p, e in factors.items())
    roots = signed_roots(factors)
    axis_stride = prod(p ** ((e + 1) // 2) for p, e in factors.items())
    assert len({z for z in roots if 0 < z < n and z % 2 == 0 and z != n - 1}) % 2 == 1
    if direct:
        assert sorted(roots) == [x for x in range(n) if x * x % n == 1]
    log = []

    def oracle(n, req):
        answer = crt_oracle(n, req, roots, axis_stride)
        if direct:
            assert answer == enumerated_oracle(n, req)
        log.append(answer)
        return answer

    result = split(n, oracle)
    m, r = 2, 0
    for main, axes in zip(log[::2], log[1::2]):
        m *= 2
        public = int((n - 1) % m == r)
        parity = sum(z % m == r and z != n - 1 for z in roots) % 2
        assert (main - axes - 2 * public) % 4 == 2 * parity
        if axes and "axis_control" not in controls:
            controls["axis_control"] = {
                "N": n,
                "m": m,
                "r": r,
                "main": main,
                "baseline": axes,
                "public_pairs": public,
                "hidden_parity": parity,
            }
        if not parity:
            r += m // 2
    assert r == result["root"] and r in roots and 1 < r < n - 1
    a, b = result["factors"]
    assert a * b == n and gcd(a, b) == 1 and a > 1 and b < n
    assert result["queries"] == 2 * ((n - 1).bit_length() - 1)
    stats["inputs"] += 1
    stats["queries"] += len(log)
    stats["max_queries"] = max(stats["max_queries"], len(log))
    stats["max_bits"] = max(stats["max_bits"], n.bit_length())
    stats["with_repeated_prime"] += max(factors.values()) > 1
    stats["at_least_three_distinct"] += len(factors) >= 3
    return result


def narrow_silence(n, p):
    c = isqrt(4 * n)
    delta = 4 * n - c * c
    if delta <= 4:
        c -= 1
        delta = 4 * n - c * c
    cap = isqrt(delta // 3)
    b = delta % 2
    q = (delta + b * b) // 4
    assert n % p == 0 and gcd(c, n) == 1 and delta > 4 * p * p * cap
    form = Request((1, b, q), cap, c)
    hits = {
        (x, y)
        for y, lo, hi in rows(n, form)
        for x in range(lo, hi + 1)
        if (x or y) and (x * x + b * x * y + q * y * y) % n == 0
    }
    v = ((c - b) // 2, 1)
    vb = (v[0] + b, -1)
    public = {
        (k * x, k * y)
        for x, y in (v, vb)
        for k in range(-isqrt(cap), isqrt(cap) + 1)
        if k
    }
    assert hits == public
    return {
        "N": n,
        "small_factor": p,
        "c": c,
        "delta": delta,
        "cap": cap,
        "signed_points": len(hits),
    }


def generate():
    controls, small, large = {}, Counter(), Counter()
    lookup = {}
    for n in range(3, 2001, 2):
        f = {int(p): int(e) for p, e in factorint(n).items()}
        if len(f) < 2:
            continue
        ans = checked_split(f, True, small, controls)
        if len(f) == 2 and max(f.values()) == 1:
            bit = ans["root"] % 4
            key = (n % 840, tuple(sorted(p % 8 for p in f)))
            if (
                key in lookup
                and lookup[key]["root"] % 4 != bit
                and "local_control" not in controls
            ):
                controls["local_control"] = [lookup[key], ans]
            lookup[key] = ans
    complete = Counter()
    for n in range(1, 2001):
        expected = dict(factorint(n))

        def native(modulus, req):
            complete["queries"] += 1
            return enumerated_oracle(modulus, req)

        assert factor(n, native) == expected
        complete["inputs"] += 1
        complete["odd_prime_powers"] += n % 2 == 1 and len(expected) == 1
        complete["even_inputs"] += n % 2 == 0
    historical = Counter()
    base = Path(__file__).resolve().parents[1]
    cases = json.loads((base / "ellipse_row_batches" / "public_cases.json").read_text())
    collision = []
    for n, *_ in cases:
        result = checked_split(dict(factorint(n)), True, historical, controls)
        if n in (11009, 28829):
            collision.append(result)
    rng = Random(353)
    for i in range(72):
        low = int(nextprime(rng.randrange(7, 3000)))
        bits = rng.randrange(80, 186)
        high = int(nextprime((1 << (bits - 1)) | rng.getrandbits(bits - 1)))
        factors = {low: 1 + i % 3, high: 1 + (i % 8 == 0)}
        for j in range(i % 4):
            p = int(nextprime(rng.randrange(4000 + 5000 * j, 6000 + 5000 * j)))
            factors[p] = 1 + (i % 5 == 0)
        checked_split(factors, False, large, controls)
    f = {3: 1, 10007: 1}
    work = Counter()
    native_unbalanced = split(30021, lambda n, q: enumerated_oracle(n, q, work))
    assert native_unbalanced["factors"] == sorted(f)
    # The baseline can have a different primitive normalization even at the first child.
    n = 147
    controls["example_forms"] = [chart(n, 4, 0, k).__dict__ for k in (False, True)]
    return {
        "original_near_square_silence": [
            narrow_silence(10713, 3),
            narrow_silence(652055, 5),
        ],
        "complete_public_factorizations": dict(complete),
        "historical_public_sample": dict(historical),
        "collision_pair": collision,
        "exact_enumeration_to_2000": dict(small),
        "seeded_factored_oracle_only": dict(large),
        "seed": 353,
        "native_unbalanced": {**native_unbalanced, **work},
        **controls,
    }


if __name__ == "__main__":
    result = generate()
    data = json.dumps(result, indent=2, sort_keys=True) + "\n"
    Path(__file__).with_name("general_scope_results.json").write_text(data)
    print(data, end="")
