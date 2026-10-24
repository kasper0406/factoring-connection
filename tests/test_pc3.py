import numpy as np
import pytest
from sympy import factorint

from src.pc3 import toy
from src.pc3 import dc


def test_tables_exact():
    X, y, Y = 200_000, 50, 300
    s, p = toy.build_tables(X, y, Y)
    rng = np.random.default_rng(0)
    for n in rng.integers(1, X + 1, size=1500):
        n = int(n)
        f = factorint(n)
        sm = n == 1 or max(f) <= y
        big = [q for q in f if q > y]
        par = len(big) == 1 and f[big[0]] == 1 and Y < big[0] <= 2 * Y
        assert s[n] == sm
        assert p[n] == par
    assert not s[0] and not p[0]


def test_label_basis_and_box():
    for q, r in [(701, 123), (503, 0), (997, 996)]:
        e1, e2 = toy.label_basis(q, r)
        assert abs(e1[0] * e2[1] - e1[1] * e2[0]) == q
        for e in (e1, e2):
            assert (e[0] - r * e[1]) % q == 0
        box = toy.child_box(e1, e2, 6)
        keys = set()
        for i, j, a, b in box:
            assert (a - r * b) % q == 0
            assert np.gcd(a, b) == 1
            assert (a, b) == toy.canon((a, b))
            keys.add((a, b))
        assert len(keys) == len(box)


def test_pair_lattice_is_unit_graph():
    # {(G_g(P'), G_g(P''))} = {n1 = lam n2 mod det}, lam a unit (note §2.3)
    from math import gcd
    rng = np.random.default_rng(1)
    for _ in range(30):
        while True:
            P1 = tuple(int(v) for v in rng.integers(-40, 41, 2))
            P2 = tuple(int(v) for v in rng.integers(-40, 41, 2))
            D = abs(P1[0] * P2[1] - P1[1] * P2[0])
            if D > 1 and gcd(*P1) == 1 and gcd(*P2) == 1:
                break
        g = rng.integers(-500, 501, size=(400, 2))
        v1 = g[:, 1] * P1[0] - g[:, 0] * P1[1]
        v2 = g[:, 1] * P2[0] - g[:, 0] * P2[1]
        lams = [l for l in range(D) if gcd(l, D) == 1 and np.all((v1 - l * v2) % D == 0)]
        assert len(lams) == 1


def test_local_model_basic():
    s, _ = toy.build_tables(300_000, 50, 300)
    tau = toy.tau_table(s, 50_000, 250_000, 50)
    assert abs(toy.local_model(1, tau, 50)) < 1e-12
    assert abs(toy.local_model(53, tau, 50) - (1 / (1 - 1 / 53) - 1)) < 1e-12
    assert toy.local_model(12, tau, 50) > 0


def test_prop_dc_inequality_random():
    # Prop DC: E[Z 1{N=0}] <= E[Z] (K phi + (sig^2 + sig)/(K-1)) for every joint law and K > 1
    rng = np.random.default_rng(2)
    for _ in range(300):
        n = int(rng.integers(5, 60))
        p = rng.dirichlet(np.ones(n))
        Z = rng.poisson(rng.uniform(0.5, 30), n).astype(float) * (rng.random(n) < 0.9)
        N = rng.poisson(rng.uniform(0.2, 10), n).astype(float)
        if (p * Z).sum() == 0 or (p * N).sum() == 0:
            continue
        lhs = (p * Z * (N == 0)).sum()
        for K in (1.5, 2.0, 5.0, 40.0):
            assert lhs <= dc.dc_bound(p, Z, N, K) * (1 + 1e-9) + 1e-12
        assert lhs <= dc.dc_bound_opt(p, Z, N) * (1 + 1e-9) + 1e-12
