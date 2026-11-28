"""End-to-end scalar oracle checks and independent frame enumeration."""

from collections import Counter
from math import gcd, isqrt, prod
from random import Random
from time import perf_counter

from sympy import jacobi_symbol

if __package__:
    from .check_mixed_ternary_embeddings import _all_frames
    from .packed_ternary_decision import (
        brute_sphere_count,
        brute_sphere_decision,
        brute_sphere_parity,
        recover_packed_factor,
        recover_packed_symmetric_factor,
    )
else:
    from check_mixed_ternary_embeddings import _all_frames
    from packed_ternary_decision import (
        brute_sphere_count,
        brute_sphere_decision,
        brute_sphere_parity,
        recover_packed_factor,
        recover_packed_symmetric_factor,
    )


def check_recovery(n):
    # Input-only selection. Small fixtures make a direct search for x cheap.
    x = next(k for k in range(n) if jacobi_symbol(-(k * k + 1), n) == -1)
    start = perf_counter()
    result = recover_packed_factor(n, x, brute_sphere_decision)
    elapsed = perf_counter() - start
    assert 1 < result.factor < n and n % result.factor == 0
    short_start = perf_counter()
    short = recover_packed_symmetric_factor(n, x, brute_sphere_decision)
    short_time = perf_counter() - short_start
    assert short.factor == result.factor and short.crt_modulus > n**6

    # Grade the claimed scalar answers from a disjoint calculation. This
    # table is constructed after the real supplier has completed the recovery.
    frames = [(u, v) for u, v in _all_frames(n, x) if u == tuple(sorted(map(abs, u)))]
    squares = {t * t for _, v in frames for t in v if t and t * t != n * n}
    symmetric = set()
    for _, v in frames:
        s, t, w = (c * c for c in v)
        symmetric.add((s * t + t * w + w * s, s * t * w))
    assert symmetric == {(0, 0), short.coefficients}
    assert 1 <= len(squares) <= 3
    assert result.coordinate_product == prod(squares)
    scale = result.packing_scale
    m = result.squared_radius
    assert m == (x + scale * n) ** 2 + 1 and scale // 7 < 8 * (x + 1) * (n + 1) <= scale
    assert result.crt_modulus > n**6
    assert prod(r.prime for r in result.readings) > n**24

    order = []
    yes_counts = {}
    category = set()
    allowed = {}
    ambiguous = []
    for row in result.readings:
        ell = row.prime
        assert gcd(ell, 210 * m * n * (x * x + 1)) == 1
        modulus = scale * ell
        allowed[modulus] = Counter(
            tuple((i - scale * j) % modulus for i, j in zip(u, v)) for u, v in frames
        )
        observed = {t % ell for t in squares} - {0, n * n % ell}
        assert observed == row.squares
        if len(observed) == len(squares):
            assert all(count == 1 for count in allowed[modulus].values())
        ambiguous.extend(
            (ell, residue, count) for residue, count in allowed[modulus].items() if count > 1
        )
        if len({t % ell for t in squares}) < len(squares):
            category.add("collision")
        if any(t % ell == 0 for t in squares):
            category.add("zero")
        if any((t - n * n) % ell == 0 for t in squares):
            category.add("axis")

    def grade_scalar_answer(norm, modulus, residue):
        assert norm == m
        assert gcd(modulus, *residue) == 1
        assert (sum(c * c for c in residue) - m) % modulus == 0
        # These conditions provide nonsingular local lifts at every divisor
        # of the odd modulus, even for queries whose exact answer is false.
        order.append((modulus, residue))
        answer = brute_sphere_decision(norm, modulus, residue)
        assert answer == (residue in allowed[modulus])
        assert brute_sphere_count(norm, modulus, residue) == allowed[modulus][residue]
        if answer:
            yes_counts[modulus] = yes_counts.get(modulus, 0) + 1
        return answer

    assert recover_packed_factor(n, x, grade_scalar_answer) == result
    assert len(order) == len(set(order)) == sum(row.queries for row in result.readings)
    assert all(yes_counts[scale * row.prime] == row.yes for row in result.readings)
    assert sum(yes_counts.values()) < len(order)
    assert recover_packed_factor(n, x, brute_sphere_parity).factor == result.factor

    def singleton_only(norm, modulus, residue):
        assert norm == m
        return allowed[modulus][residue] == 1

    assert recover_packed_factor(n, x, singleton_only).factor == result.factor

    def all_no(norm, modulus, residue):
        assert norm == m and next(replay_order) == (modulus, residue)
        return False

    replay_order = iter(order)
    try:
        recover_packed_factor(n, x, all_no)
    except ValueError as error:
        assert "inconsistent" in str(error)
    else:
        raise AssertionError("accepted an inconsistent all-negative oracle")
    assert next(replay_order, None) is None
    short_queries = sum(row.queries for row in short.readings)
    print(
        f"N={n} x={x} factor={result.factor} primes={len(result.readings)} "
        f"queries={len(order)} yes={sum(yes_counts.values())} "
        f"losses={','.join(sorted(category)) or 'none'} seconds={elapsed:.3f}; "
        f"short={len(short.readings)} primes/{short_queries} queries/{short_time:.3f}s; "
        f"ambiguous={ambiguous[:1]}"
    )
    return category


def check_scalar_supplier():
    rng = Random(41027)
    checks = 0
    for m in range(61):
        r = isqrt(m)
        points = [
            (i, j, k)
            for i in range(-r, r + 1)
            for j in range(-r, r + 1)
            for k in range(-r, r + 1)
            if i * i + j * j + k * k == m
        ]
        for modulus in (1, 2, 3, 5, 9):
            residues = Counter(tuple(t % modulus for t in w) for w in points)
            samples = [
                tuple(rng.randrange(-2 * modulus, 2 * modulus) for _ in range(3)) for _ in range(5)
            ]
            for z in [*residues, *samples]:
                assert brute_sphere_decision(m, modulus, z) == (
                    tuple(t % modulus for t in z) in residues
                )
                assert brute_sphere_count(m, modulus, z) == residues[tuple(t % modulus for t in z)]
                checks += 1
    return checks


if __name__ == "__main__":
    print("Independent scalar decisions:", check_scalar_supplier())
    reductions = set()
    for n in (15, 35, 65, 77, 91, 183):
        reductions.update(check_recovery(n))
    assert reductions == {"axis", "collision", "zero"}
