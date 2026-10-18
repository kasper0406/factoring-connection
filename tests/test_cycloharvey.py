import random
from fractions import Fraction
from math import gcd

import flint
import pytest

from src.common.instances import balanced_semiprime
from src.cycloharvey import calibration as cal
from src.cycloharvey import covering as cov
from src.cycloharvey import lattice as lat
from src.cycloharvey import operational as op


# ---------------------------------------------------------------- calibration
def test_lp_exponents_exact():
    assert cal.exponent_lp("lehman")[0] == Fraction(1, 3)
    assert cal.exponent_lp("hittmeir")[0] == Fraction(2, 9)
    assert cal.exponent_lp("harvey")[0] == Fraction(1, 5)
    assert cal.exponent_lp("hh")[0] == Fraction(1, 5)
    assert cal.exponent_lp("bsgs_sum")[0] == Fraction(1, 4)


def test_ceil_sqrt_ratio():
    for num in range(0, 400):
        for den in (1, 2, 3, 7, 16):
            w = cal.ceil_sqrt_ratio(num, den)
            assert den * w * w >= num and (w == 0 or den * (w - 1) ** 2 < num)


def test_sum_ceil_sqrt_bruteforce():
    N = 10 ** 13 + 7
    for r in (1, 3, 40, 300):
        d = cal.divisor_counts(r)
        P = cal.prefix(d)
        for C in (16, 16 * r * r, 16 * r * r * 25):
            bf = sum(int(d[k]) * cal.ceil_sqrt_ratio(N, C * k) for k in range(1, r + 1))
            assert cal.sum_ceil_sqrt(N, C, r, P) == bf


def test_harvey_count_matches_pair_loop():
    N = 2 ** 40 + 2 ** 20 + 1
    r, m = 30, 200
    brute = sum(max(1, cal.ceil_sqrt_ratio(N, 16 * r * r * m * m * a * b))
                for a in range(1, r + 1) for b in range(1, r // a + 1))
    assert cal.cost_harvey(N, r, m)["giant"] == brute


# ---------------------------------------------------------------- lattices
def test_wronskian():
    for I in lat.all_supports():
        assert lat.wronskian_det(I) == lat.wronskian_formula(I) != 0


def test_kappa_prop33_and_general():
    N = 1000003 * 1000033
    for M in (5, 7, 11):
        for r in range(1, M):
            assert lat.class_lattice_index((-1, 1), N, r, M)["kappa"] == 3
            assert lat.class_lattice_index((-1, 1, 2), N, r, M)["kappa"] == 6


def test_prop33_lock_constant_mod_M2():
    # (a, b) with b r^2 == a N (mod M) => a q + b p constant mod M^2 on p == r (M)
    M, N = 7, 1000003 * 1000033
    r = 3
    for a in range(1, 6):
        b = (a * N * pow(r * r, -1, M)) % M
        vals = set()
        for x in range(60):
            pp = r + M * x
            qq = (N * pow(pp, -1, M * M)) % (M * M)
            vals.add((a * qq + b * pp) % (M * M))
        assert len(vals) == 1


def test_guess_neutrality_formula():
    for I in lat.all_supports():
        e0 = lat.support_stats(I)["exponent"]
        assert lat.exponent_with_guesses(I, Fraction(1, 7)) == e0
        assert lat.exponent_with_guesses(I, Fraction(1, 7), D_exponent=Fraction(1, 7)) == e0 + Fraction(1, 7)


def test_min_heuristic_exponent_is_lehman_support():
    rows = sorted(lat.all_supports(), key=lambda I: lat.support_stats(I)["exponent"])
    assert rows[0] == (-1, 1)
    assert lat.support_stats(rows[0])["exponent"] == Fraction(1, 5)
    assert all(lat.support_stats(I)["exponent"] > Fraction(1, 5) for I in rows[1:])


# ---------------------------------------------------------------- covering
def _brute_ok(t, x0, x1, W, M=1, r=0):
    vals = [cov.evaluate(t, p) for p in range(x0, x1 + 1) if (p - r) % M == 0]
    if not vals:
        return True
    mn, mx = min(vals), max(vals)
    U = -((-mn.numerator) // mn.denominator)
    top = mx.numerator // mx.denominator
    return top - U + 1 <= W


def test_oscillation_exact_vs_bruteforce():
    rng = random.Random(1)
    N = 10007 * 10009
    for _ in range(200):
        I = rng.choice([(-1, 1), (-2, -1, 1, 2), (1, 2), (-3, -1, 1, 3)])
        a = [rng.randint(-3, 3) for _ in I]
        if not any(a):
            continue
        t = cov.combine(cov.monomial_basis(I, N), a)
        x0 = rng.randint(7000, 10000)
        x1 = x0 + rng.randint(0, 60)
        W = rng.choice([1, 4, 50, 10 ** 6])
        for M, r in ((1, 0), (5, 2)):
            assert cov.oscillation_ok(t, x0, x1, W, M=M, r=r)[0] == _brute_ok(t, x0, x1, W, M, r)


def test_sweep_is_a_valid_cover():
    N = 2 ** 26 + 2 ** 13 + 1
    plo, phi = cov.balanced_range(N)
    for I in [(-1, 1), (-2, -1, 1, 2)]:
        basis = cov.monomial_basis(I, N)
        res = cov.sweep(basis, plo, phi, 8, record=True)
        nxt = plo
        for x0, x1, v, U in res["targets"]:
            assert x0 == nxt
            t = cov.combine(basis, v)
            assert _brute_ok(t, x0, x1, 8)
            nxt = x1 + 1
        assert nxt > phi


# ---------------------------------------------------------------- operational
@pytest.mark.parametrize("bits", [32, 40])
def test_harvey_search_factors(bits):
    for seed in range(3):
        inst = balanced_semiprime(bits, seed)
        N = inst.N
        b = cal.best_harvey(N)
        res = op.harvey_search(N, b["r"], b["m"], seed, early_exit=False)
        assert res["factor"] in (inst.labels["p"], N // inst.labels["p"])
        assert res["met"].get("giant") is None or res["met"]["giant"] == b["giant"]


def test_cyclo_inverse_and_handle_identity():
    inst = balanced_semiprime(32, 3)
    N, p = inst.N, inst.labels["p"]
    q = N // p
    for m in (5, 8, 12):
        R = op.Cyclo(N, m)
        rng = random.Random(m)
        a = R.elt([rng.randrange(N) for _ in range(R.D)])
        ai = R.inv(a)
        assert R.coeffs(R.mul(a, ai)) == [1] + [0] * (R.D - 1)
        F = {(0, 1): 3, (1, 0): -2, (1, 1): 1}
        T = op.handle(R, a, p % m, F)
        e = 3 * q - 2 * p + p * q
        Rp = op.Cyclo(p, m)
        f = 1
        while pow(p, f, m) != 1:
            f += 1
        lhs = Rp.coeffs(Rp.pow(Rp.elt(R.coeffs(a)), e % (p ** f - 1)))
        assert lhs == [x % p for x in R.coeffs(T)]


def test_multipoint_eval():
    N = 2 ** 61 - 1
    rng = random.Random(0)
    f = flint.nmod_poly([rng.randrange(N) for _ in range(20)], N)
    pts = [rng.randrange(N) for _ in range(13)]
    met = op.NullMetrics()
    assert op.multipoint_eval(f, pts, N, met) == [int(f(x)) for x in pts]


def test_cyclo_search_factors_small():
    inst = balanced_semiprime(32, 1)
    N = inst.N
    b = cal.best_harvey(N)
    res = op.cyclo_harvey_search(N, 5, max(2, b["r"]), b["m"], 0)
    assert res["factor"] in (inst.labels["p"], N // inst.labels["p"])
