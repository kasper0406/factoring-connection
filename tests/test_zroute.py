"""Tests for src/zroute (correctness of the deterministic Z route pieces; toy sizes)."""
import math
import random
from math import gcd

import flint
import pytest

from src.common.instances import balanced_semiprime
from src.detGRH.core import FactorBase, fb_bound, grh_bound
from src.zroute import core


def test_gauss_reduce_is_shortest_and_in_lattice():
    rng = random.Random(1)
    for _ in range(200):
        N = rng.randrange(10 ** 6, 10 ** 7) | 1
        lam = rng.randrange(1, N)
        (u, v), (u2, v2) = core.gauss_reduce(N, lam)
        assert (u - lam * v) % N == 0 and (u2 - lam * v2) % N == 0
        assert abs(u * v2 - u2 * v) == N                     # a basis
        # shortest: no lattice vector with small coefficients is shorter
        n0 = u * u + v * v
        for i in range(-3, 4):
            for j in range(-3, 4):
                if i or j:
                    a, b = i * u + j * u2, i * v + j * v2
                    assert a * a + b * b >= n0
        assert n0 <= 2 * N / math.sqrt(3) + 1


def test_base_m_field():
    N = balanced_semiprime(36, 3).N
    fld = core.base_m_field(N, 3)
    assert sum(c * fld.m ** i for i, c in enumerate(fld.f)) == N and fld.f[3] == 1
    # F(a, b) = b^3 f(a/b); psi(a - b alpha) = a - b m
    for a, b in [(5, 7), (-11, 3), (2, 1)]:
        assert fld.F(a, b) == sum(c * a ** i * b ** (3 - i) for i, c in enumerate(fld.f))


def test_xkz_relations_verify_and_rank_formula():
    N = balanced_semiprime(32, 5).N
    fb = FactorBase.make(N, fb_bound(N))
    fld = core.base_m_field(N, 3)
    ctr = {}
    rows = core.nfs_rows(fld, fb, 600, 120, 10 ** 9, ctr)[:400]
    zs = core.xkz(fld, rows, ctr)
    assert ctr.get("NONUNIT", 0) == 0 and ctr["xk_undecided"] == 0
    rels = core.rational_relations(fb, rows, zs, ctr)
    assert ctr["rel_verify_fail"] == 0
    for o in rels:
        assert core.ev_signed(fb, o) == 1
    # torsion dimension = dim Lambda_1 - unit rank seen
    assert ctr["torsion_vectors"] == ctr["dim_Lambda1"] - ctr["unit_rank_seen"]


def test_dominance_design_forces_rank_and_factors():
    inst = balanced_semiprime(36, 7)
    N = inst.N
    fb = FactorBase.make(N, fb_bound(N))
    rows, H, tries = core.dominance_rows(N, fb, 3000, {})
    assert len(rows) == len(fb.primes)
    assert core.is_diag_dominant(fb, rows)
    res = core.pin_and_exponent(fb, [r.rat for r in rows], grh_bound(N))
    assert res["full"] and res["rank"] == fb.n
    g, base = core.miller_split(N, res["D"])
    assert g is not None and N % g == 0 and 1 < g < N


def test_levy_desplanques_random():
    """Strict row diagonal dominance => nonsingular (the algebra used by Lemma DD1)."""
    rng = random.Random(2)
    for _ in range(50):
        n = rng.randrange(2, 9)
        M = [[rng.randrange(-3, 4) for _ in range(n)] for _ in range(n)]
        for i in range(n):
            M[i][i] = sum(abs(M[i][j]) for j in range(n) if j != i) + rng.randrange(1, 3)
        assert flint.fmpz_mat(M).det() != 0


def test_pin_exponent_kills_pinned_primes():
    N = balanced_semiprime(32, 9).N
    fb = FactorBase.make(N, fb_bound(N))
    ctr = {}
    rows = core.orbit_rows(N, fb, 2, N.bit_length() + 1, 3000, 200, ctr)
    res = core.pin_and_exponent(fb, [r.rat for r in rows], grh_bound(N))
    for q in fb.primes:
        if q <= grh_bound(N) and q not in res["unpinned"]:
            assert pow(q, res["D"], N) == 1


def test_hidden_labels_rank_toy():
    """Lemma HR toy: random large labels give full rank as soon as there are n rows."""
    N = balanced_semiprime(28, 1).N
    fb = FactorBase.make(N, fb_bound(N))
    rng = random.Random(3)
    rows = []
    while len(rows) < fb.n:
        x = [rng.randrange(N ** 2) for _ in range(fb.n)]
        r = core.split_row(N, fb, core.ev_fast(fb, x), x, "random")
        if r is not None:
            rows.append(r.rat)
    assert flint.fmpz_mat(rows).rank() == fb.n


def test_ev_fast_matches():
    N = balanced_semiprime(30, 2).N
    fb = FactorBase.make(N, 200)
    rng = random.Random(4)
    for _ in range(5):
        x = [rng.randrange(1000) for _ in range(fb.n)]
        assert core.ev_fast(fb, x) == core.ev_signed(fb, x)


def test_operational_core_has_no_randomness():
    import inspect
    src = inspect.getsource(core)
    assert "import random" not in src and "random." not in src.replace("random-exponent", "")


def test_wronskian_labels_force_rank():
    """Lemma WR (codex zroute_1 §3.4): labels x_{j,i} = (d/dt)^i t^(n+j) at an integer t above the Cauchy bound
    make (X - V)^T Z of full rank n for every accepted subset, kernel Z (rank >= n) and descent part V."""
    from math import factorial, perm
    rng = random.Random(5)
    n, P, H0, B0 = 3, 6, 3, 4
    A0 = H0 * ((n + P) ** n + P * B0)
    C = factorial(n) * (P + 1) ** n * A0 ** n
    t = C + 2
    def label(j):
        return [perm(n + j, i) * t ** (n + j - i) for i in range(n)]
    for _ in range(300):
        S = sorted(rng.sample(range(1, P + 1), rng.randrange(n, P + 1)))
        M = len(S)
        while True:
            Z = [[rng.randrange(-H0, H0 + 1) for _ in range(n)] for _ in range(M)]   # M x n, columns = kernel vectors
            if flint.fmpz_mat(Z).rank() == n:
                break
        V = [[rng.randrange(-B0, B0 + 1) for _ in range(n)] for _ in range(M)]
        X = [label(j) for j in S]
        R = [[X[i][r] - V[i][r] for r in range(n)] for i in range(M)]               # M x n rational parts
        O = flint.fmpz_mat(Z).transpose() * flint.fmpz_mat(R)                       # n x n outputs
        assert O.rank() == n
