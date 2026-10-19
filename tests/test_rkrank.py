import random
from fractions import Fraction
from itertools import product
from math import gcd, isqrt

from src.rkrank import gf2, sampler as smp
from src.rkrank.harness import peel_core
from src.rkrank.toynfs import ToyField, gen_functionals, real_root_intervals, below_root, peval
from src.nfsfruit.instances import composite


def _span(vecs):
    s = {0}
    for v in vecs:
        s |= {x ^ v for x in s}
    return s


def test_echelon_basis_below_bruteforce():
    rng = random.Random(1)
    for _ in range(40):
        n, k, t = 9, rng.randrange(1, 8), rng.randrange(1, 9)
        vecs = [rng.getrandbits(n) for _ in range(k)]
        E = gf2.Echelon(vecs)
        S = _span(vecs)
        assert E.rank == len(S).bit_length() - 1
        low = {x for x in S if x < (1 << t)}
        assert _span(E.basis_below(t)) == low


def test_dim_in_coords_and_annihilated_dim():
    rng = random.Random(2)
    for _ in range(40):
        n = 8
        vecs = [rng.getrandbits(n) for _ in range(rng.randrange(1, 7))]
        B = gf2.Echelon(vecs).basis()
        mask = rng.getrandbits(n)
        S = _span(B)
        assert 2 ** gf2.dim_in_coords(B, mask) == sum(1 for x in S if x & ~mask == 0)
        funcs = gf2.Echelon([rng.getrandbits(n) for _ in range(4)]).basis()
        ann = [x for x in _span(funcs) if all(gf2.dot(x, b) == 0 for b in B)]
        assert 2 ** gf2.annihilated_dim(funcs, B) == len(ann)


def test_real_root_isolation_and_sign():
    rng = random.Random(3)
    f = [-7, 3, -2, 1]          # x^3 - 2x^2 + 3x - 7
    ivs = real_root_intervals(f)
    for lo, hi in ivs:
        assert peval(f, lo) * peval(f, hi) < 0
    for _ in range(200):
        a, b = rng.randrange(-60, 61), rng.randrange(1, 30)
        for lo, hi in ivs:
            # exact reference: bisect to width < 1/(10^6 b)
            L, H = lo, hi
            while H - L > Fraction(1, 10 ** 6 * b):
                M = (L + H) / 2
                if (peval(f, M) > 0) == (peval(f, H) > 0):
                    H = M
                else:
                    L = M
            q = Fraction(a, b)
            if q < L or q > H:
                assert below_root(a, b, f, (lo, hi)) == (q < L)


def _small_field():
    N, fs = composite("blum", 44, 5)
    return N, fs, ToyField(N, 1500, 12, 0)


def test_rows_parities_and_W_in_lambda_half():
    N, fs, T = _small_field()
    rows, pairs, _ = T.sieve(3000, 600, int(2.0 * T.ncols))
    assert len(rows) > T.ncols
    for r, (a, b) in zip(rows[:50], pairs[:50]):
        t = a - b * T.m
        ft = T.fb.factor(t)
        for p in T.ps:
            assert (r >> T.col[("r", p)]) & 1 == ft.get(p, 0) & 1
        assert r & 1 == (t < 0)
    E = gf2.Echelon(rows)
    Wb = E.basis_below(T.D_rat)
    gen = gen_functionals(T, fs)
    assert gf2.all_annihilate(gen, Wb)            # W ⊆ Λ/2 (harness check)
    high = ((1 << T.ncols) - 1) ^ ((1 << T.D_rat) - 1)
    assert len(Wb) == E.rank - gf2.rank(r & high for r in rows)


def test_peel_core_rows_not_in_any_dependency():
    rng = random.Random(4)
    for _ in range(30):
        n, k = 7, 9
        rows = [rng.getrandbits(n) for _ in range(k)]
        core = set(peel_core(rows, (1 << n) - 1))
        for S in product([0, 1], repeat=k):
            if any(S):
                acc = 0
                for i, s in enumerate(S):
                    if s:
                        acc ^= rows[i]
                if acc == 0:
                    assert {i for i, s in enumerate(S) if s} <= core


def test_cell_points_complete_and_sampler_valid():
    N = 10007 * 10009
    rng = random.Random(5)
    h = smp.h_two_thirds(N)
    T, H = isqrt(h) // 2, h // 2
    for _ in range(20):
        x1 = rng.randrange(1, N)
        got = sorted(smp.cell_points(N, x1, T, H))
        ref = []
        for t in range(-T, T + 1):
            s = (x1 * x1 + 2 * x1 * t) % N
            if s > N // 2:
                s -= N
            if abs(s) <= H:
                ref.append(t)
        if got:                                   # skewed cells may be skipped (returns [])
            assert got == ref
    for x, z in smp.small_squares(N, 200, 1, h):
        assert (x * x - z) % N == 0 and abs(z) <= h
