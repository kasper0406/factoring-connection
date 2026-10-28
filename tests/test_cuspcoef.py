import ast
import os
import random
from fractions import Fraction
from itertools import product
from math import isqrt, sqrt

import numpy as np
import pytest

from src.cuspcoef import decode as D
from src.cuspcoef.family import family, j_invariant, CM_J

ROOT = os.path.join(os.path.dirname(__file__), "..")
FORBIDDEN_MODULES = {"src.cuspcoef.harness", "cypari2", "sympy", "src.common.instances", "gmpy2", "flint"}
FORBIDDEN_NAMES = {"factorint", "factor", "isprime", "is_prime", "primefactors", "ellap", "ellcard"}


@pytest.mark.parametrize("fn", ["src/cuspcoef/decode.py", "src/cuspcoef/family.py"])
def test_no_leakage(fn):
    tree = ast.parse(open(os.path.join(ROOT, fn)).read())
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                assert a.name not in FORBIDDEN_MODULES, (fn, a.name)
        if isinstance(node, ast.ImportFrom):
            assert node.module not in FORBIDDEN_MODULES, (fn, node.module)
        if isinstance(node, ast.Call):
            f = node.func
            name = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", None)
            assert name not in FORBIDDEN_NAMES, (fn, name)


def test_class_counts_bruteforce():
    for l in (3, 5):
        n = np.zeros((l, l))
        for a, b, c, d in product(range(l), repeat=4):
            det = (a * d - b * c) % l
            if det:
                n[(a + d) % l, det] += 1
        assert np.array_equal(n, D.class_counts(l))


def test_family():
    F = family()
    assert len(F) >= 32000
    js = [j_invariant(A, B) for A, B in F[:2000]]
    assert len(set(js)) == len(js)
    assert all(not (j.denominator == 1 and j.numerator in CM_J) for j in js)


def test_candidate_pairs_and_spf():
    spf = D.SPF(10 ** 5)
    assert spf.decompose(2 * 2 * 3 * 7 * 7 * 97) == {2: 2, 3: 1, 7: 2, 97: 1}
    C = D.candidate_pairs(-84, 20, spf)
    assert set(C) == {(x, -84 // x) for x in range(-20, 21) if x and 84 % abs(x) == 0 and abs(84 // x) <= 20}


def test_decoding_primes_and_crt():
    for X in (2 ** 15, 2 ** 19, 2 ** 23):
        ells = D.decoding_primes(X)
        P = 1
        for l in ells:
            P *= l
        assert P > 2 * X and P // ells[-1] <= 2 * X
    x, m = D.crt([1, 2, 3], [3, 5, 7])
    assert m == 105 and x % 3 == 1 and x % 5 == 2 and x % 7 == 3


def test_decode_toy_end_to_end():
    from src.cuspcoef import harness as H
    N, p, q = H.instance(1000, 32)
    F = family()[:2000]
    A = [H.ap_AB(a, b, p) * H.ap_AB(a, b, q) for a, b in F]
    X = 1 << 15
    spf = D.SPF(8 * X)
    res = D.decode(N, 32, F, A, mode="full", spf=spf)
    assert res["p"] in (p, q)


# ------------------------------------------------------------------ Theorem DFB: the likelihood-ratio inequality

def _gekeler_law(p):
    from src.cuspcoef import harness as H
    T = isqrt(4 * p - 1)
    law = {t: Fraction(H.hurwitz(4 * p - t * t), 12 * 2 * p) for t in range(-T, T + 1)}
    assert sum(law.values()) == 1          # Kronecker-Hurwitz
    return law


def test_dfb_ratio_bound_exact():
    """Exact check of the pointwise inequality P(A = a) <= B_p B_q (2R+1)^2 / sqrt(N) * S(a) behind Theorem DFB,
    for A = a_p a_q under the Gekeler law, and S the law of x'y' with x', y' uniform on [-R, R]."""
    X = 512
    p, q = 607, 1009
    N = p * q
    R = isqrt(8 * X)
    Lp, Lq = _gekeler_law(p), _gekeler_law(q)
    Bp = max(Lp.values()) * sqrt(p)
    Bq = max(Lq.values()) * sqrt(q)
    PA = {}
    for x, u in Lp.items():
        for y, v in Lq.items():
            PA[x * y] = PA.get(x * y, 0) + u * v
    S = {}
    w = Fraction(1, (2 * R + 1) ** 2)
    for x in range(-R, R + 1):
        for y in range(-R, R + 1):
            S[x * y] = S.get(x * y, 0) + w
    C = float(Bp * Bq) * (2 * R + 1) ** 2 / sqrt(N)
    worst = max(float(PA[a] / S[a]) for a in PA)
    assert all(a in S for a in PA)
    assert worst <= C * (1 + 1e-12)
    assert C <= 32 * float(Bp * Bq) * 1.2


def test_cong_decoder_runs():
    rng = random.Random(1)
    N = 40009 * 49991
    A = [rng.randint(-400, 400) for _ in range(200)]
    res = D.decode(N, 32, family()[:200], A, mode="cong", max_combos=10)
    assert res["tried"] <= 10
