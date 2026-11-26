"""Independent checks on grouped traces, matched inputs, and a fixed cubic probe."""

import json
from argparse import ArgumentParser
from collections import Counter
from itertools import islice
from math import gcd, isqrt
from pathlib import Path

from class_trace import (
    ClassTraces,
    UnsupportedClassGroup,
    branch_from_n,
    compose,
    evaluate,
    evaluate_branch,
    field_bits,
)
from flint import fmpz_mod_poly_ctx
from row_batches import Branch
from sympy import factorint, isprime, primerange

from genusleak.forms import prime_form, reduce_form, sqrt_mod_prime

HERE = Path(__file__).resolve().parent


def balanced_cases(delta):
    """Control generation only: p is known here; readers get N alone."""
    start = max(101, delta // 2)
    for p in primerange(start, start + 1_500_000):
        if pow(-delta, (p - 1) // 2, p) != 1:
            continue
        root = sqrt_mod_prime(-delta, p)
        lo = isqrt(4 * p * p - delta) + 1
        hi = isqrt(8 * p * p - delta)
        for r in (root, p - root):
            odd = r if r % 2 else r + p
            c = lo + (odd - lo) % (2 * p)
            if c <= hi:
                q, rem = divmod(c * c + delta, 4 * p)
                if rem == 0 and p < q <= 2 * p and isprime(q):
                    assert isqrt(4 * p * q) == c
                    yield p, q, c


def mixed_class(p, q, c):
    f = reduce_form((p, c, q))
    return compose(f, f)


def norm_readout(model, n):
    counts, work = model.norm_profile(n)
    assert counts[0] + 2 * sum(counts[1:]) == 4
    return model.transform(counts), dict(work)


def other_case(model, cutoffs, primary, calibrations):
    lo, hi = cutoffs
    estimates = {f: i for i, f in enumerate(model.forms)}
    primitive = {i for i, (_, order, _) in enumerate(model.orbits) if order == model.h}
    small = set(range(len(model.orbits))) - primitive
    for p, q, c in balanced_cases(model.delta):
        e = estimates[mixed_class(p, q, c)]
        u = [model.unordered(k, e) for k, _, _ in model.orbits]
        if any(primary[i] != u[i] for i in small):
            continue
        if all(primary[i] == u[i] for i in primitive):
            continue
        answers = [model.contract(labels[0], w, u) for labels, w in calibrations]
        if answers[0] != answers[1]:
            continue
        measured, work = norm_readout(model, p * q)
        assert u == measured
        return {
            "N": p * q,
            "control_factors": [p, q],
            "c": c,
            "same_proper_character_traces": True,
            "absolute_cutoffs": [lo, hi],
            "band_parities": answers,
            "proper_trace_values": [field_bits(u[i]) for i in sorted(small)],
            "readout": work,
        }
    raise RuntimeError("no matched control in the search range")


def local_ideal_labels(model, cutoff):
    """Independent coefficient construction by factoring only the small control levels."""
    logs = {f: i for i, f in enumerate(model.forms)}
    labels = [0] * model.h
    for t in range(cutoff // 4 + 1, cutoff + 1):
        if t % 2 == 0 or isqrt(t) ** 2 == t:
            continue
        exponents = [0]
        for p, v in factorint(t).items():
            if pow(-model.delta, (p - 1) // 2, p) != 1:
                exponents = []
                break
            step = v * logs[prime_form(-model.delta, p)]
            exponents = [
                (e + sign * step) % model.h for e in exponents for sign in (-1, 1)
            ]
        for e in exponents:
            labels[e] ^= 1
    return labels


def check_fixture(case):
    n = case["N"]
    br, shortcut = branch_from_n(n)
    assert shortcut is None
    cutoffs = [r["B"] for r in case["cutoffs"]]
    expected = [r["parity"] for r in case["cutoffs"]]
    run = evaluate(n, cutoffs)
    assert [r["band_parity"] for r in run["cutoffs"]] == expected == [0, 1]
    m = ClassTraces(br.delta)
    u, _ = norm_readout(m, n)
    p = case["cutoffs"][-1]["factor_witness"]
    assert 1 < p < n and n % p == 0
    e_control = m.forms.index(mixed_class(p, n // p, br.c))
    assert u == [m.unordered(k, e_control) for k, _, _ in m.orbits]
    calibrations = [m.calibration(b) for b in cutoffs]
    for b, (labels, w) in zip(cutoffs, calibrations, strict=True):
        ideal_labels = local_ideal_labels(m, b)
        for e in range(m.h):
            measured = m.contract(
                labels[0], w, [m.unordered(k, e) for k, _, _ in m.orbits]
            )
            assert measured == ideal_labels[e]
    differences = [
        a ^ b for a, b in zip(calibrations[0][0], calibrations[1][0], strict=True)
    ]
    jump = {e for e in range(m.h) if differences[min(e, (-e) % m.h)]}
    stabilizer = [s for s in range(m.h) if {(j + s) % m.h for j in jump} == jump]
    assert stabilizer == [0]
    spectral = m.transform(differences)
    rank = differences[0] + sum(
        2 * d for (_, _, d), v in zip(m.orbits, spectral, strict=True) if v
    )
    return {
        **run,
        "all_class_checks": len(cutoffs) * m.h,
        "jump_size": len(jump),
        "jump_translation_stabilizer_size": len(stabilizer),
        "jump_translate_rank": rank,
        "paired_control": other_case(m, cutoffs, u, calibrations),
    }


def frobenius_gcd(n, exponent):
    """Safe Euclidean probe over Z/nZ; a nonunit leading coefficient is observable."""
    ring = fmpz_mod_poly_ctx(n)
    x = ring.gen()
    f = x**3 + 2 * x + 1
    a, b = f, x.pow_mod(exponent, f) - x
    while b:
        lead = int(b.leading_coefficient())
        d = gcd(lead, n)
        if d != 1:
            return None, d
        b = b * pow(lead, -1, n)
        a, b = b, a % b
    return a.degree(), None


def cubic_controls():
    m = ClassTraces(59)
    assert m.h == 3
    hist = Counter()
    examples = {}
    for p, q, c in islice(balanced_cases(59), 128):
        n = p * q
        u, _ = norm_readout(m, n)
        bit = int(u[0])
        deg_p, _ = frobenius_gcd(p, p)
        deg_q, _ = frobenius_gcd(q, q)
        assert deg_p == deg_q == (0 if bit else 3)
        assert bit == int(mixed_class(p, q, c) != m.forms[0])
        deg, divisor = frobenius_gcd(n, n)
        if divisor:
            assert 1 < divisor < n and n % divisor == 0
        key = f"trace_{bit}_probe_{deg if deg is not None else 'factor'}"
        hist[key] += 1
        examples.setdefault(
            key,
            {"N": n, "c": c, "control_factors": [p, q], "divisor_from_probe": divisor},
        )
    assert (
        sum(hist.values()) == 128
        and hist["trace_0_probe_0"]
        and hist["trace_1_probe_0"]
    )
    return {"cases": sum(hist.values()), "histogram": hist, "first_examples": examples}


def sample_controls():
    rows = json.loads(
        (HERE.parent / "ellipse_row_batches/public_cases.json").read_text()
    )
    stats = Counter()
    for args in rows:
        br = Branch(*args)
        if br.delta <= 3 or br.delta % 8 != 3 or not isprime(br.delta):
            stats["outside_prime_discriminant_family"] += 1
            continue
        try:
            run = evaluate_branch(br, range(br.cap + 1))
        except UnsupportedClassGroup as err:
            stats[str(err)] += 1
            continue
        # Literal positive-quadrant disk. Only primitive integral basis coordinates
        # can contribute to the nonsquare band for these fundamental discriminants.
        norms = []
        radius = 4 * br.n * br.cap
        for y in range(1, isqrt(radius // br.delta) + 1):
            for w in range(1, isqrt(radius - br.delta * y * y) + 1):
                t, rem = divmod(w * w + br.delta * y * y, 4 * br.n)
                if rem == 0 and gcd((w - y) // 2, y) == 1 and isqrt(t) ** 2 != t:
                    norms.append(t)
        for r in run["cutoffs"]:
            b = r["B"]
            assert r["band_parity"] == sum(b < 4 * t and t <= b for t in norms) % 2
            stats["cutoff_checks"] += 1
            stats["positive_cutoffs"] += r["band_parity"]
        stats["eligible_branches"] += 1
    return {"sample_branches": len(rows), **stats}


def run():
    rows = json.loads((HERE / "conductor_band_results.json").read_text())["cases"]
    target = {138468065327, 39866375420117}
    checks = [check_fixture(c) for c in rows if c["N"] in target]
    assert len(checks) == 2
    return {
        "prime_discriminants": checks,
        "small_sample": sample_controls(),
        "order_three": cubic_controls(),
    }


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    opt = parser.parse_args()
    data = json.dumps(run(), indent=2) + "\n"
    if opt.output:
        opt.output.write_text(data)
    else:
        print(data, end="")
