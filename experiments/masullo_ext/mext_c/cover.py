"""Run from scratch with the repository Python, with bytecode writes disabled."""
import itertools
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "experiments/subthird/codex"), str(ROOT)]
import experiment as ex
from flint import acb
from sympy import integer_nthroot, primerange

ALPHABET = list(itertools.product(range(-2, 3), repeat=2))


def pair_units():
    units = set()
    for a, b in ALPHABET:
        norm = a*a - 2*b*b
        if norm == 0:
            continue
        for c, d in ALPHABET:
            u, v = a*c-2*b*d, a*d-b*c
            if u % norm or v % norm:
                continue
            u, v = u//norm, v//norm
            if abs(u*u - 2*v*v) == 1:
                units.add((u,v))
    return sorted(units)


def unit_key(z, units):
    images = []
    for u, v in units:
        w = tuple(t for a, b in zip(z[::2], z[1::2])
                  for t in (a*u + 2*b*v, a*v + b*u))
        if max(map(abs, w)) <= 2:
            images.append(w)
    return min(images)


def add_independent(basis, z):
    if ex.kranks(basis + [z]) > len(basis):
        basis.append(z)


def prefix(row, points):
    multiples = [[v*j for j in range(-2, 3)] for v in row]
    ordered = []
    for z in points:
        value = abs(sum((tab[j+2] for tab, j in zip(multiples, z)), acb(0)))
        ordered.append((value.lower(), value.upper(), z))
    ordered.sort(key=lambda item: item[0])
    provisional = []
    for lower, upper, z in ordered:
        add_independent(provisional, z)
        if len(provisional) == 3:
            threshold = lower
            break
    else:
        raise ValueError("Coefficient box has deficient rank")
    basis = []
    for lower, upper, z in ordered:
        if lower >= threshold:
            break
        if upper < threshold:
            add_independent(basis, z)
    assert len(basis) < 3
    return basis


def run(bits, mode, units):
    N = ex.fixture(bits, 19)
    _, b, P, Q = ex.select(N)
    f, g = P*P-2, Q*Q-2
    U, V = [ex.expand(X, X*X-2, ex.generators(X, mode, 3)) for X in (P,Q)]
    residues = [int(u(b)) % N for u in U]
    assert residues == [int(v(b)) % N for v in V]
    points = [z for z in itertools.product(range(-2,3), repeat=6)
              if next((j for j in z if j), 0) > 0 and
              math.gcd(sum(t*c for t, c in zip(z,residues)), N) == 1]
    rows, _ = ex.rows_at(ex.coeffs(f), list(map(ex.coeffs,U)),
                         ex.coeffs(g), list(map(ex.coeffs,V)), 192)
    spaces = [prefix(row, points) for row in rows]
    cap = 3**32 * 2**((bits+15)//16)
    fourth, exact = integer_nthroot(bits*(bits-1).bit_length()**3, 4)
    B = 2**(int(fourth) + int(not exact))
    primes = list(primerange(B+1))
    primorial = math.prod(primes)
    fields = [ex.Ideals(X,primes) for X in (f,g)]
    capped_vectors, keys, uncovered_keys = 0, set(), set()
    covered_rows, other_rows, columns = [], [], {}
    for z in points:
        a, w = ex.element(U,z), ex.element(V,z)
        v1, v2 = abs(int(f.resultant(a))), abs(int(g.resultant(w)))
        value = v1*v2
        if not 0 < value <= cap:
            continue
        capped_vectors += 2
        key = unit_key(z, units)
        if key in keys:
            continue
        keys.add(key)
        covered = any(ex.kranks(basis+[z]) == len(basis) for basis in spaces)
        if not covered:
            uncovered_keys.add(key)
        if ex.smooth(value,primorial):
            mask = 0
            for side, (alpha, norm) in enumerate(((a,v1), (w,v2))):
                for prime in fields[side].parity(alpha,norm,primes):
                    col = columns.setdefault((side,prime),len(columns))
                    mask |= 1 << col
            (covered_rows if covered else other_rows).append(mask)
    covered_rank = ex.rank(covered_rows)
    full_rank = ex.rank(covered_rows + other_rows)
    return dict(bits=bits, N=str(N), mode=mode, B=B, cap=str(cap),
                valid_vectors=2*len(points), capped_vectors=capped_vectors,
                capped_unit_classes=len(keys), uncovered_classes=len(uncovered_keys),
                smooth_classes=len(covered_rows)+len(other_rows),
                smooth_uncovered=len(other_rows), covered_rank=covered_rank,
                smooth_rank=full_rank, incremental_rank=full_rank-covered_rank,
                prefix_ranks=[len(s) for s in spaces])


if __name__ == "__main__":
    assert ex.n == 8
    units = pair_units()
    for bits, mode in itertools.product((72,96), ("divided","power")):
        print(json.dumps(run(bits,mode,units)), flush=True)
