"""Tests for src/tracedata (phase 25)."""
import ast
import os
import random

import numpy as np
import pytest
import sympy

from src.tracedata import decode as D
from src.tracedata import galois as G
from src.tracedata import harness as H
from src.tracedata import reductions as R
from src.tracedata.xe_descent import xe_factor

ROOT = os.path.join(os.path.dirname(__file__), "..", "src", "tracedata")
OPERATIONAL = ("reductions.py", "decode.py", "xe_descent.py", "galois.py")
FORBIDDEN_MODULES = {"src.tracedata.harness", "cypari2", "sympy", "src.ecount.harness", "src.common.instances"}
FORBIDDEN_NAMES = {"factorint", "factor", "isprime", "is_prime", "primefactors", "ellap", "ellcard"}


@pytest.mark.parametrize("fn", OPERATIONAL)
def test_operational_no_leakage(fn):
    tree = ast.parse(open(os.path.join(ROOT, fn)).read())
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                assert a.name not in FORBIDDEN_MODULES, (fn, a.name)
        if isinstance(node, ast.ImportFrom):
            assert node.module not in FORBIDDEN_MODULES, (fn, node.module)
            for a in node.names:
                assert a.name not in FORBIDDEN_NAMES, (fn, a.name)
        if isinstance(node, ast.Call):
            f = node.func
            name = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", None)
            assert name not in FORBIDDEN_NAMES, (fn, name)


def _pairs(lo, hi, k, seed):
    rng = random.Random(seed)
    P = list(sympy.primerange(lo, hi))
    out = []
    while len(out) < k:
        p, q = sorted(rng.sample(P, 2))
        out.append((p, q))
    return out


def test_a_pk_matches_pari():
    E = H.ell("C2_37a1")
    for p in sympy.primerange(41, 300):
        a = H.ap("C2_37a1", p)
        for k in (2, 3):
            assert H.a_pk(a, p, k) == int(H.pari.ellak(E, p ** k))


@pytest.mark.parametrize("name", ["C1_11a1", "C2_37a1", "C3_2_3", "CM1_32a2", "T7_26b1"])
def test_bc2_bc3_mix_tor(name):
    m = H.TORSION_M.get(name)
    for (p, q) in _pairs(10 ** 4, 10 ** 5, 25, name):
        N = p * q
        a, b = H.ap(name, p), H.ap(name, q)
        aN = a * b
        g = R.bc2_factor(N, aN, H.a_pk(a, p, 2) * H.a_pk(b, q, 2))
        if (a, b) != (0, 0):
            assert g in (p, q)
        g = R.bc3_factor(N, aN, H.a_pk(a, p, 3) * H.a_pk(b, q, 3))
        if aN != 0:
            assert g in (p, q)
        M = (p + 1 - a) * (q + 1 - b)
        assert R.mix_factor(N, sorted((a, b)), M) in (p, q)
        if m:
            assert R.s_mod_m_from_torsion(N, aN, m) == (p + q) % m


def test_bc2_both_supersingular_fails_cm():
    # CM1 at p = q = 3 mod 4: a_p = a_q = 0; (a_N, a_{N^2}) = (0, N) carries nothing
    p, q = 10007, 10039
    assert p % 4 == 3 and q % 4 == 3
    assert H.ap("CM1_32a2", p) == 0 and H.ap("CM1_32a2", q) == 0
    assert R.bc2_factor(p * q, 0, (-p) * (-q)) is None


def test_eta_weight3_and_cm2():
    c = H.eta4_6_coeffs(400)
    for p in sympy.primerange(3, 400):
        assert c[p] == H.cm_weight3_ap(p) == H.cm_weight3_ap_fast(p)
    for (p, q) in _pairs(10 ** 4, 10 ** 5, 60, "cm2"):
        if p % 4 != 1 or q % 4 != 1:
            continue
        N = p * q
        aN = H.ap("CM1_32a2", p) * H.ap("CM1_32a2", q)
        g3 = H.cm_weight3_ap(p) * H.cm_weight3_ap(q)
        assert R.linked_pair_factor(N, aN, g3, c=1, e=2) in (p, q)
        assert R.cm_split_bit(aN)


def test_cm_split_bit():
    for (p, q) in _pairs(10 ** 3, 10 ** 4, 80, "cmqr"):
        N = p * q
        if N % 4 != 1:
            continue
        aN = H.ap("CM1_32a2", p) * H.ap("CM1_32a2", q)
        assert R.cm_split_bit(aN) == (p % 4 == 1)


def test_class_counts_bruteforce():
    for l in (3, 5, 7):
        assert (G.class_counts(l) == G.class_counts_bruteforce(l)).all()
        assert G.class_counts(l).sum() == (l * l - 1) * (l * l - l)


def test_psi3_lines_match_pari():
    for name in ("C2_37a1", "C5_92"):
        for p in sympy.primerange(5, 1500):
            if H.disc(name) % p == 0:
                continue
            r = H.roots_psi3(name, p)
            if r == 2:
                assert p % 3 == 2
            if r in (1, 4):
                assert p % 3 == 1
    T = G.psi3_table()
    assert sum(T.values()) == 48
    assert {r for (_, r) in T} == {0, 1, 2, 4}


def test_info_and_advantages_sane():
    for l in (3, 5, 7):
        assert 0 < G.info_ordered(l) < 1
        assert 0 < G.info_product(l) < G.info_ordered(l)
    # cross-check adv_aN_beta(3) by brute force over GL2(F_3)^2
    els = list(G.gl2_elements(3))
    joint = {}
    for g in els:
        for h in els:
            dg = (g[0] * g[3] - g[1] * g[2]) % 3
            dh = (h[0] * h[3] - h[1] * h[2]) % 3
            if dg != dh:
                continue
            key = ((g[0] + g[3]) * (h[0] + h[3])) % 3
            joint.setdefault(key, [0, 0])[0 if dg == 1 else 1] += 1
    tot = sum(a + b for a, b in joint.values())
    adv = sum(max(v) for v in joint.values()) / tot - 0.5
    assert abs(adv - G.adv_aN_beta(3)) < 1e-12


def test_decoder_on_model_samples():
    rng = np.random.default_rng(1)
    for l in (5, 7, 11):
        n = G.class_counts(l).astype(float)
        for d in range(1, l):
            t = rng.choice(l, size=60 * l * l, p=n[:, d] / n[:, d].sum())
            assert D.rank_p(l, t)[0] == d


def test_crt_and_end_to_end_synthetic():
    # traces drawn from the model with the true residues of p: the end-to-end decoder must return p
    rng = np.random.default_rng(2)
    p, q = 1000003, 1000033
    N = p * q
    ells = [3, 5, 7, 11, 13, 17]
    # build an integer trace list whose residues mod each l follow the model for d = p mod l (CRT per index)
    K = 20000
    cols = []
    for l in ells:
        n = G.class_counts(l).astype(float)[:, p % l]
        cols.append(rng.choice(l, size=K, p=n / n.sum()))
    traces = [D._crt([int(c[i]) for c in cols], ells)[0] for i in range(K)]
    assert D.factor_from_ordered(N, traces, ells, top=2) == p


def test_xe_descent_splits_generic():
    a, b = 2, 3
    p, q = 10007, 10009
    Ep = H.pari.ellinit([0, 0, 0, a, b])
    M = int(H.pari.ellcard(Ep, p)) * int(H.pari.ellcard(Ep, q))
    res = xe_factor(p * q, a, b, M, random.Random(0), B=100, max_trials=20)
    assert res["factor"] in (p, q)


def test_scalar_twist_symmetry_products():
    # Theorem SC: the per-curve law of t_p t_q given (p, q) = (d, c/d) mod l is invariant under d -> lambda^2 d
    for l in (5, 7, 11, 13):
        L = G.prod_likelihood(l)
        for c in range(1, l):
            for d in range(1, l):
                for lam in range(2, l):
                    d2 = lam * lam * d % l
                    assert np.allclose(L[d, c * pow(d, -1, l) % l], L[d2, c * pow(d2, -1, l) % l])
        # identifiable classes: square classes if (c/l) = 1, everything if (c/l) = -1
        for c in range(1, l):
            cls = G.identifiable_classes(l, "product", c)
            assert len(cls) == (2 if G.legendre(c, l) == 1 else 1)
            assert len(G.identifiable_classes(l, "single", c)) == l - 1
            # unordered pairs and M mod l identify exactly {d, c/d}, i.e. s mod l
            want = sorted(sorted({d, c * pow(d, -1, l) % l}) for d in range(1, l))
            want = [list(x) for x in {tuple(w) for w in want}]
            for kind in ("unordered", "M"):
                got = sorted(sorted(x) for x in G.identifiable_classes(l, kind, c))
                assert got == sorted(want), (l, c, kind)


def test_cm_quartic_pair_and_family():
    H.CURVES.setdefault("q2", [0, 0, 0, -2, 0])
    Ds = list(range(2, 40))
    for D in Ds:
        H.CURVES.setdefault(f"q{D}", [0, 0, 0, -D, 0])
    n = 0
    for (p, q) in _pairs(10 ** 5, 10 ** 6, 200, "cmpair"):
        if p % 4 != 1 or q % 4 != 1:
            continue
        n += 1
        N = p * q
        A1 = H.ap("CM1_32a2", p) * H.ap("CM1_32a2", q)
        if p % 8 == 5 and q % 8 == 5:
            assert R.cm_quartic_pair_factor(N, A1, H.ap("q2", p) * H.ap("q2", q)) in (p, q)
        A = [H.ap(f"q{D}", p) * H.ap(f"q{D}", q) for D in Ds]
        assert R.cm_quartic_family_factor(N, A1, A) in (p, q)
    assert n > 5


def test_laws_normalised():
    for l in (5, 11):
        for kind in ("single", "ordered", "unordered", "product", "M"):
            for d in range(1, l):
                assert abs(G.law(l, kind, d, 3 % l).sum() - 1) < 1e-12, (l, kind, d)


def test_mix_degenerate_case_documented():
    # cx2: E: y^2 = x^3 - 7x - 2, p = 5, q = 7: a_p = a_q = 1, M = N; the quadratic is identically zero
    assert R.mix_factor(35, [1, 1], 35) is None


def test_ext_reduction_any_curve_including_XE():
    from src.tracedata.ext_reduction import ext_factor
    pari = H.pari
    cases = [(-13392, -1080432, 54443, 62819), (27331, 2, 151, 181), (-1296, 11664, 1000003, 1000033),
             (-1, 1, 999983, 1000037)]
    for a, b, p, q in cases:
        E = pari.ellinit([0, 0, 0, a, b])
        M = int(pari.ellcard(E, p)) * int(pari.ellcard(E, q))
        r = ext_factor(p * q, a, b, M, random.Random(3), max_trials=200)
        assert r["factor"] in (p, q)
