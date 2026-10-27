"""Tests for src/nativemod (phase 29)."""
import ast
import os
import random
from fractions import Fraction

import gmpy2
import pytest

from src.nativemod import ecalg
from src.nativemod import harness as H
from src.nativemod import reductions as R

ROOT = os.path.join(os.path.dirname(__file__), "..", "src", "nativemod")
OPERATIONAL = ("reductions.py", "ecalg.py")
FORBIDDEN_MODULES = {"src.nativemod.harness", "cypari2", "sympy", "src.common.instances", "src.tracedata.harness",
                     "src.ecount.harness"}
FORBIDDEN_NAMES = {"factorint", "factor", "isprime", "is_prime", "is_probab_prime", "primefactors", "ellap",
                   "ellcard", "next_prime"}


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


def _primes(lo, hi):
    return [p for p in range(lo, hi) if gmpy2.is_prime(p)]


def _pairs(k, lo, hi, seed):
    ps = _primes(lo, hi)
    rng = random.Random(seed)
    out = []
    while len(out) < k:
        p, q = sorted(rng.sample(ps, 2))
        out.append((p, q))
    return out


# --- E1a: Eichler-Selberg Hurwitz moments --------------------------------------------------------------------
@pytest.mark.parametrize("N,p,q", [(15, 3, 5), (35, 5, 7), (77, 7, 11), (221, 13, 17), (391, 17, 23)])
def test_hurwitz_moments(N, p, q):
    assert H.hurwitz_moment(N, 2) == 2 * (N + q)
    for k in (4, 6, 8, 10, 14):
        assert H.hurwitz_moment(N, k) == -2 * (1 + p ** (k - 1))
    for k in (12, 16):
        assert H.hurwitz_moment(N, k) == -2 * H.trace_TN_level1(N, k) - 2 * (1 + p ** (k - 1))


# --- E1b: Hecke matrices on M_2(Gamma_0(M)) ----------------------------------------------------------------
@pytest.mark.parametrize("M", [11, 37])
def test_hecke_full_space(M):
    for p, q in [(3, 5), (5, 7), (13, 17)]:
        if M in (p, q):
            continue
        N = p * q
        TN, Tp, Tq = H.hecke_full(M, N), H.hecke_full(M, p), H.hecke_full(M, q)
        assert TN == Tp * Tq
        assert int(H.pari.subst(H.pari.charpoly(TN), H.pari("x"), H.sigma(N))) == 0
        I = H.pari.matid(int(H.pari.matsize(TN)[0]))
        assert TN * TN - H.hecke_full(M, N * N) - N * I == p * H.hecke_full(M, q * q) + q * H.hecke_full(M, p * p)


# --- E1c: theta series ---------------------------------------------------------------------------------------
def test_theta_level11_and_class_number_one():
    B = 200
    t1, t2 = H.theta_counts(H.L11_1, B), H.theta_counts(H.L11_2, B)
    f = H.coeffs_11a1(B)
    assert all(t1[n] - t2[n] == 4 * f[n] for n in range(1, B + 1))
    # Brandt/Eichler-weighted combination (post hoc E1c(ii')): (3 theta_1 + 2 theta_2)/5 = E_11
    assert all(3 * t1[n] + 2 * t2[n] == 12 * H.sigma(n) for n in range(1, B + 1) if n % 11)
    z4 = H.theta_counts([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], B, even=False)
    assert all(z4[n] == 8 * H.sigma(n) for n in range(1, B + 1, 2))
    d4 = H.theta_counts(H.D4_HURWITZ, B)
    assert all(d4[n] == 24 * H.sigma(n) for n in range(1, B + 1, 2))


def test_merel_count_matches_bruteforce():
    from src.geometry.modcurve import heilbronn_merel
    for n in (2, 3, 6, 15, 21, 35):
        assert H.merel_size(n) == len(heilbronn_merel(n))
    # the b = 0 part alone is sigma(n): factor-free definition, divisor content inside
    assert sum(n // a for a in range(1, 36) if 35 % a == 0) == H.sigma(35)


# --- reductions ------------------------------------------------------------------------------------------------
def _chi4(n):
    return 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)


def _one(n):
    return 1


@pytest.mark.parametrize("k,chi,chip", [(2, _one, _one), (4, _one, _one), (12, _one, _one), (3, _one, _chi4),
                                        (3, _chi4, _one), (2, _chi4, _chi4), (5, _one, _chi4)])
def test_eis_factor(k, chi, chip):
    for p, q in _pairs(15, 10 ** 5, 2 * 10 ** 5, k):
        N = p * q
        a = sum(chi(N // d) * chip(d) * d ** (k - 1) for d in (1, p, q, N))
        assert R.eis_factor(N, a, k, chi, chip) in (p, q)


def test_eis_weight1_is_bounded():
    # weight 1, (1, chi_-4): a(N) = r_2(N)/4 in {0, 4}: O(1) bits, no factor
    vals = set()
    for p, q in _pairs(30, 10 ** 4, 2 * 10 ** 4, 1):
        N = p * q
        vals.add(sum(_one(N // d) * _chi4(d) for d in (1, p, q, N)))
    assert vals <= {0, 4, -4, 2, -2}


def test_crt_trace_factor_and_promise_exception():
    ok = 0
    for p, q in _pairs(6, 2 ** 20, 2 ** 20 + 2 ** 18, 7):
        N = p * q
        ap, aq = H.ap("37a1", p), H.ap("37a1", q)
        t = (ap * q * pow(q, -1, p) + aq * p * pow(p, -1, q)) % N
        ok += R.crt_trace_factor(N, t) in (p, q)
    assert ok == 6
    # y^2 = x^3 - x, p = q = 3 mod 4: a_p = a_q = 0, t = 0 carries nothing
    p, q = 1048583, 1048639
    assert p % 4 == 3 and q % 4 == 3 and H.ap("32a2", p) == 0 and H.ap("32a2", q) == 0


def test_class_pair_factor():
    ls = _primes(3, 200)
    for p, q in _pairs(10, 10 ** 6, 2 * 10 ** 6, 3):
        N = p * q
        pairs = {l: (p % l, q % l) for l in ls if N % l}
        assert R.class_pair_factor(N, pairs) in (p, q)


# --- ecalg: the literal power map --------------------------------------------------------------------------
def test_power_map_trace_identity():
    psi = H.divpoly("37a1", 5)
    for p, q in [(1009, 1013), (10007, 10009), (65537, 65539)]:
        N = p * q
        g = ecalg.monic_division_poly(psi, 5, N)
        st = ecalg.statistics(N, g)
        gp = ecalg.monic_division_poly(psi, 5, p)
        assert st["TrY"] % p == ecalg.power_sums(gp, q + 1, p)[q]   # Tr(x^N) = s_q (mod p)
        gq = ecalg.monic_division_poly(psi, 5, q)
        assert st["TrY"] % q == ecalg.power_sums(gq, p + 1, q)[p]   # Tr(x^N) = s_p (mod q)


def test_power_map_not_additive():
    # codex's counterexample: (1 + 1)^35 != 1^35 + 1^35 in F_5
    assert pow(2, 35, 5) == 3 and (pow(1, 35, 5) + pow(1, 35, 5)) % 5 == 2


def test_power_sums_newton():
    g = [3, -2, 0, 1]   # x^3 - 2x + 3
    s = ecalg.power_sums(g, 6)
    import numpy as np
    roots = np.roots([1, 0, -2, 3])
    for k in range(6):
        assert abs(sum(r ** k for r in roots).real - s[k]) < 1e-6


def test_crt_trace_equal_traces_exception():
    # cx2 / proof critic: a_p = a_q != 0 gives t = a_p, and gcd(t - a_q, N) = N: no split (exception set {a_p = a_q})
    for curve, p, q in (("37a1", 1013, 1097), ("11a1", 1009, 1019)):
        ap, aq = H.ap(curve, p), H.ap(curve, q)
        assert ap == aq != 0
        t = (ap * q * pow(q, -1, p) + aq * p * pow(p, -1, q)) % (p * q)
        assert t == ap % (p * q)
        assert R.crt_trace_factor(p * q, t) is None


def test_full_gcd_checker_catches_coefficient_split():
    # cx2: N = 67*103 had an unchecked coefficient gcd in v1; the v1.1 checker raises Split on it or earlier
    psi = H.divpoly("37a1", 3)
    N = 67 * 103
    g = ecalg.monic_division_poly(psi, 3, N)
    with pytest.raises(ecalg.Split):
        ecalg.statistics(N, g, check_all=True)


def test_sigma_mod_N_is_rsa_type():
    # sigma_{k-1}(N) mod N = 1 + s^{k-1} (mod N): cross terms of (p+q)^{k-1} are multiples of pq
    for p, q in [(101, 103), (1009, 1013)]:
        N = p * q
        for k in (2, 4, 6):
            assert H.sigma(N, k - 1) % N == (1 + pow(p + q, k - 1, N)) % N


def test_mod5_leak_11a1():
    # proof critic: 11a1 has rational 5-torsion, so a(N) = (1+p)(1+q) (mod 5)
    for p, q in [(13, 17), (101, 103), (1009, 1013)]:
        assert (H.ap("11a1", p) * H.ap("11a1", q) - (1 + p) * (1 + q)) % 5 == 0
