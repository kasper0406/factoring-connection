import random
from math import gcd

import pytest

from src.genusleak import cayley, cf, forms, sqrtred


def _pari():
    import cypari2
    return cypari2.Pari()


def test_reduce_and_compose_match_pari():
    pari = _pari()
    rng = random.Random(1)
    for N in (10403, 1000003 * 1000033, 999983 * 1000037):
        D = -4 * N
        ells = forms.split_primes(D, 8)
        gs = [forms.prime_form(D, l) for l in ells]
        for _ in range(30):
            f, g = rng.choice(gs), rng.choice(gs)
            f = forms.power(f, rng.randrange(1, 50))
            g = forms.power(g, rng.randrange(1, 50))
            h = forms.compose(f, g)
            assert forms.disc(h) == D and forms.is_reduced(h)
            ph = pari.qfbred(pari.qfbcomp(pari.Qfb(*f), pari.Qfb(*g)))
            assert tuple(int(x) for x in pari.Vec(ph)[:3]) == h


def test_power_class_number_is_identity():
    pari = _pari()
    N = 10007 * 10009
    D = -4 * N
    h = int(pari.qfbclassno(D))
    for l in forms.split_primes(D, 5):
        assert forms.power(forms.prime_form(D, l), h) == forms.identity(D)


def test_ambiguous_forms_and_public():
    p, q = 10009, 30011  # N = 1 mod 4? check below
    N = p * q
    D = -4 * N
    pub = forms.public_ambiguous(N)
    assert (1, 0, N) in pub
    if N % 4 == 1:
        assert (2, 2, (N + 1) // 2) in pub
    f = forms.reduce_form((p, 0, q))
    assert forms.is_ambiguous(f) and forms.ambiguous_split(N, f) in (p, q)


def test_squfof_and_midpoint():
    N = 10007 * 10037
    s = cf.squfof_scan(N, 1)
    assert s["L"] > 0
    for i in s["proper"]:
        assert i % 2 == 0
    # midpoint: when it splits, the divisor is correct
    for k in (1, 2, 3, 5, 7):
        r = cf.midpoint_split(N, k)
        if r["split"]:
            assert N % r["split"] == 0 and 1 < r["split"] < N
    assert any(cf.midpoint_split(N, k)["split"] for k in (1, 2, 3, 5, 6, 7, 10, 11))


def test_free_walks_and_girth_lemma():
    assert cayley.free_closed_walks(1, 6) == [1, 0, 2, 0, 6, 0, 20]
    N = 1009 * 1013
    D, ells, gens, elems, idx, table = cayley.setup(N, 1, 4)
    W = cayley.closed_walks(table, 16)
    Wf = cayley.free_closed_walks(len(ells), 16)
    g0 = cayley.girth_bound(N, 1, ells)
    assert all(W[t] == Wf[t] for t in range(min(g0, 17)))
    assert all(W[t] >= Wf[t] for t in range(17))


def test_sqrt_reduction_toy():
    N = 1009 * 1013  # 1 mod 4
    assert N % 4 == 1
    D = -4 * N
    ells = forms.split_primes(D, 8)
    gens = [forms.prime_form(D, l) for l in ells]
    elems, idx, table = cayley.enumerate_subgroup(D, gens)
    roots = {}
    for z in elems:
        y = forms.compose(z, z)
        if y not in roots or z < roots[y]:
            roots[y] = z
    rng = random.Random(5)
    hits = 0
    for _ in range(40):
        out = sqrtred.sqrt_reduction_query(N, D, gens, lambda y: roots[y], rng)
        if out["split"]:
            assert out["split"] in (1009, 1013)
            hits += 1
    assert 8 <= hits <= 32


def test_lemma0_ambiguous_forms():
    """Lemma 0: for N > 12 k^2 every reduced ambiguous form of D = -4kN either has a | 4k (public) or yields
    a factor of N by gcd; and the useful ones are at least half of Cl[2]."""
    from src.genusleak import genus
    pari = _pari()
    for (p, q) in ((101, 113), (103, 131), (107, 109), (211, 307), (1009, 1013), (1019, 1031)):
        N = p * q
        for k in (1, 2, 3, 5, 6, 7, 10, 11):
            D = -4 * k * N
            amb = genus.ambiguous_reduced_forms(D)
            useful = 0
            for f in amb:
                s = forms.ambiguous_split(N, f)
                if s:
                    assert s in (p, q)
                    useful += 1
                else:
                    assert (4 * k) % f[0] == 0, (p, q, k, f)
            # number of reduced ambiguous forms = |Cl[2]| (one per class, D < -4)
            cyc = [int(c) for c in pari.quadclassunit(D)[1]]
            assert len(amb) == 2 ** sum(1 for c in cyc if c % 2 == 0)
            assert 2 * useful >= len(amb)


def test_qr_via_principal_genus():
    """Prop GEN (QR_N <= GEN-SQ), with a brute-force 'is a square in Cl(D)' oracle (independent of genus theory)."""
    from src.genusleak import genus
    import gmpy2
    rng = random.Random(3)
    for (p, q, k) in ((1009, 1013, 1), (1009, 1019, 3), (211, 307, 2)):
        N = p * q
        D = -4 * k * N
        ells = forms.split_primes(D, 12)
        gens = [forms.prime_form(D, l) for l in ells]
        elems, idx, table = cayley.enumerate_subgroup(D, gens)
        assert len(elems) == int(_pari().qfbclassno(D))
        squares = {forms.compose(z, z) for z in elems}
        oracle = lambda f: f in squares
        for a in range(2, 200):
            if gmpy2.jacobi(a, N) != 1 or gcd(a, N) != 1:
                continue
            ans, m, _ = genus.qr_via_principal_genus(N, k, a, oracle, rng)
            assert ans == (gmpy2.legendre(a, p) == 1), (p, q, k, a)


def test_operational_isolation():
    """Operational modules import neither the harness nor PARI/sympy factoring code."""
    import ast, pathlib
    root = pathlib.Path(__file__).resolve().parent.parent / "src" / "genusleak"
    for name in ("forms.py", "cf.py", "cayley.py", "sqrtred.py", "genus.py"):
        tree = ast.parse((root / name).read_text())
        mods = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                mods |= {a.name for a in node.names}
            elif isinstance(node, ast.ImportFrom):
                mods.add(node.module or "")
        assert not any(m and (m.startswith("cypari2") or m.startswith("sympy") or "harness" in m
                              or m.startswith("src.common")) for m in mods), (name, mods)
