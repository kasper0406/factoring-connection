"""Tests for src/factorbits (phase 21)."""
import ast
import random
from math import gcd
from pathlib import Path

import pytest

from src.common.instances import balanced_semiprime
from src.factorbits import equiv as E
from src.factorbits import harness as H
from src.factorbits import pure as P

ROOT = Path(__file__).resolve().parents[1]


# ----------------------------------------------------------------------------- isolation
def test_equiv_is_isolated_from_factoring():
    """equiv.py (OPERATIONAL) imports neither sympy, flint nor the harness."""
    tree = ast.parse((ROOT / "src/factorbits/equiv.py").read_text())
    names = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            names |= {a.name for a in n.names}
        elif isinstance(n, ast.ImportFrom):
            names.add(n.module or "")
    assert not any(x.startswith(("sympy", "flint")) or "harness" in x for x in names), names


# ----------------------------------------------------------------------------- Prop E1 / E2
@pytest.mark.parametrize("seed", range(6))
def test_successor_ladder_exact_calls(seed):
    N = balanced_semiprime(48, seed).N
    c = E.Counter()
    v = E.factor_via_successor(N, H.succ_oracle, c)
    assert E.value(v) == N and len(v) == 2
    assert c.calls == bin(N).count("1") - 1


@pytest.mark.parametrize("w", [1, 3, 6])
def test_window_ladder_calls(w):
    N = balanced_semiprime(60, 7).N
    c = E.Counter()
    v = E.factor_via_addition_window(N, H.add_oracle, w, c)
    assert E.value(v) == N
    digits = []
    M = N
    while M:
        digits.append(M % (1 << w)); M >>= w
    assert c.calls == sum(1 for d in digits[:-1] if d)      # nonzero digits after the leading one
    assert c.calls <= -(-N.bit_length() // w)


@pytest.mark.parametrize("seed", range(10))
def test_one_call_prime_trick(seed):
    N = balanced_semiprime(64, seed).N
    c = E.Counter()
    v = E.factor_via_addition_prime(N, H.add_oracle, c, max_gap=10 ** 5)
    assert v is not None and E.value(v) == N and c.calls == 1


@pytest.mark.parametrize("seed", range(10))
def test_one_call_successor_multiplier(seed):
    N = balanced_semiprime(64, 100 + seed).N
    c = E.Counter()
    v = E.factor_via_successor_multiplier(N, H.succ_oracle, c)
    assert v is not None and E.value(v) == N and c.calls == 1


def test_addition_from_factoring_and_splitting():
    rng = random.Random(3)
    for _ in range(30):
        a, b = rng.randrange(2, 10 ** 12), rng.randrange(2, 10 ** 12)
        va, vb = H.fact(a), H.fact(b)
        c1, c2 = E.Counter(), E.Counter()
        assert E.add_via_factoring(va, vb, H.fact, c1) == H.fact(a + b) and c1.calls == 1
        v = E.add_via_splitting(va, vb, H.split_oracle, c2)
        assert v == H.fact(a + b)
        omega_big = sum(v.values())
        assert c2.calls == omega_big - 1 if omega_big > 1 else c2.calls == 0
        assert c2.calls <= (a + b).bit_length() - 1


# ----------------------------------------------------------------------------- Prop E3 partial outputs
@pytest.mark.parametrize("seed", range(4))
def test_support_oracle_suffices(seed):
    N = balanced_semiprime(48, 200 + seed).N
    c = E.Counter()
    v = E.factor_via_support(N, H.support_oracle, c)
    assert E.value(v) == N and c.calls == bin(N).count("1") - 1


@pytest.mark.parametrize("m", [2, 3])
def test_expmod_oracle_suffices_including_powers(m):
    # targets with large m-th power parts force the recursion
    for N in [balanced_semiprime(40, 300).N, (10007 ** m) * 3 * 65537, (2 ** 61 - 1) ** m * 5]:
        c = E.Counter()
        v = E.factor_via_expmod(N, H.make_expmod_oracle(m), m, c, trial=100)
        assert E.value(v) == N, (N, v)


# ----------------------------------------------------------------------------- Prop Z
@pytest.mark.parametrize("p", [1009, 7919, 65537])
def test_zech_dlog_equivalence(p):
    import sympy
    g = int(sympy.primitive_root(p))
    tab = H.make_dlog_table(p, g)
    zech = H.make_zech_oracle(p, g, tab)
    rng = random.Random(p)
    for _ in range(40):
        h = rng.randrange(1, p)
        c = E.Counter()
        assert E.dlog_via_zech(h, p, g, zech, c) == tab[h]
        assert c.calls == (bin(h).count("1") if h > 1 else 0)
        assert c.calls <= p.bit_length()
        k = rng.randrange(p - 1)
        if (1 + pow(g, k, p)) % p:
            c2 = E.Counter()
            assert E.zech_via_dlog(k, p, g, lambda t: tab[t], c2) == zech(k) and c2.calls == 1


# ----------------------------------------------------------------------------- pure trees
def _small_tree():
    M = P.smooth_exponent(20, 12)
    t1 = P.PowTest(2, M)
    t2 = P.PowTest(3, P.smooth_exponent(60, 12))
    t3 = P.ECMTest(7, 2, P.smooth_exponent(30, 12))
    t4 = P.ProductTest(3000, 3500, "x")
    # adaptive: different follow-up test depending on the first answer
    return P.node(t1, P.node(t3, P.leaf("a"), P.leaf("b")), P.node(t2, P.node(t4, P.leaf("c"), P.leaf("d")), P.leaf("e")))


def test_partition_lemma_exhaustive():
    """Run the instrumented tree mod N on all pairs of 12-bit primes: it splits iff the diagonal
    classes differ, and otherwise ends in the common class (Partition lemma)."""
    primes = P.primes_top2(12)
    tree = _small_tree()
    cls = P.classes(tree, primes)
    assert len(set(cls.values())) >= 3
    for i, x in enumerate(primes):
        for y in primes[i + 1:]:
            out = P.run_on_N(tree, x * y)
            if cls[x] == cls[y]:
                assert out == ("leaf", cls[x])
            else:
                assert out[0] == "split" and out[1] in (x, y)


def test_pt_bound_and_eps_cross_bound():
    primes = P.primes_top2(12)
    P1 = [x for x in primes if x % 4 == 1]
    P3 = [x for x in primes if x % 4 == 3]
    trees = [_small_tree(), P.nonadaptive([P.ProductTest(1, 60, "x2")]),
             P.nonadaptive([P.ProductTest(1, 2100, "x2")]), P.nonadaptive([P.ProductTest(3072, 3600)])]
    targets = [lambda x, y: (min(x, y) >> 1) & 1, lambda x, y: (min(x, y) >> 10) & 1,
               lambda x, y: (min(x, y) >> 9) & 1, lambda x, y: x % 4]
    for tree in trees:
        cls = P.classes(tree, primes)
        for F in targets:
            r = P.pt_bound_check(cls, primes, F)
            assert r["success"] <= r["bound"] + 1e-12
            assert r["p_split"] >= r["eta"] - 1e-12
        r = P.eps_cross_check(cls, P1, P3)
        assert r["success"] <= r["bound"] + 1e-12
    # the x^2+1 product over a full half-range decides eps exactly and splits every mixed N
    cls = P.classes(P.nonadaptive([P.ProductTest(1, 2100, "x2")]), primes)
    r = P.eps_cross_check(cls, P1, P3)
    assert abs(r["success"] - 1.0) < 1e-12 and abs(r["p_split_A"]) < 1e-12 and abs(r["p_split_B"] - 1) < 1e-12


def test_polynomial_self_reduction_lemma_examples():
    """Lemma PSR sanity: h(N) = N^2 + N is divisible by p (g = x), and the exponential map
    N -> 2^N - 1 carries the divisor 2^p - 1 (divisibility sequence), which no polynomial map can."""
    for (p, q) in [(11, 13), (101, 103), (1009, 2003)]:
        N = p * q
        assert (N * N + N) % p == 0
        assert (pow(2, N, 2 ** p - 1) - 1) % (2 ** p - 1) == 0
