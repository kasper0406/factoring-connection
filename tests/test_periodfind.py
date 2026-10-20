"""Tests for src/periodfind (agent periodfind)."""
import math
import random

import gmpy2
import numpy as np
import pytest

from src.periodfind import observables as ob
from src.periodfind import spectra as sp
from src.periodfind import generic as gen
from src.periodfind.sft import CountingOracle, sft, recover_order, cf_denominators

SMALL = [(23, 47), (29, 53), (37, 43), (11, 101)]


def _unit(rng, N):
    while True:
        a = rng.randrange(2, N)
        if math.gcd(a, N) == 1:
            return a


def test_orbit_and_order():
    p, q = 1009, 2003
    N = p * q
    a = 5
    info = sp.order_info(a, p, q)
    r = info["r"]
    assert pow(a, r, N) == 1
    for ell in set(gmpy2.mpz(x) for x in [2, 3, 7]):
        if r % ell == 0:
            assert pow(a, r // int(ell), N) != 1
    orb = sp.orbit(a, N, r)
    assert orb[0] == 1 and all(int(orb[t]) == pow(a, t, N) for t in (1, 17, r - 1))


@pytest.mark.parametrize("p,q", SMALL)
def test_transfer_identity(p, q):
    """Theorem T: F^(j) = sum over chi with chi(a) = e(j/r) of g^(chi)."""
    N = p * q
    rng = random.Random(p * q)
    for _ in range(2):
        a = _unit(rng, N)
        for name in ["P2_additive", "P3_interval", "P6_shifted_jacobi", "P8_lsb_square", "P12_hash"]:
            gfun = (lambda nm: (lambda y, N: ob.evaluate(nm, y, N, random.Random(1))))(name)
            if name == "P2_additive":
                gfun = lambda y, N: ob.obs_additive(y, N, 7)
            assert sp.transfer_check(N, p, q, a, gfun) < 1e-10


@pytest.mark.parametrize("p,q", SMALL + [(1009, 2003)])
def test_gauss_bound_additive(p, q):
    """Prediction W for additive characters: |F^(j)| <= sqrt(N)/r for all j != 0."""
    N = p * q
    rng = random.Random(q)
    for _ in range(3):
        a = _unit(rng, N)
        r = sp.order_info(a, p, q)["r"]
        h = _unit(rng, N)
        F = ob.obs_additive(sp.orbit(a, N, r), N, h)
        Fh = np.abs(np.fft.fft(F) / r)
        assert Fh[1:].max() <= math.sqrt(N) / r * (1 + 1e-9)


def test_shifted_jacobi_bound():
    p, q = 1009, 2003
    N = p * q
    rng = random.Random(3)
    for _ in range(3):
        a = _unit(rng, N)
        r = sp.order_info(a, p, q)["r"]
        F = ob.obs_shifted_jacobi(sp.orbit(a, N, r), N)
        Fh = np.abs(np.fft.fft(F) / r)
        assert Fh[1:].max() <= math.sqrt(N) / r * (1 + 1e-9)


def test_comb_flat():
    p, q = 1009, 2003
    N = p * q
    a = 3
    r = sp.order_info(a, p, q)["r"]
    orb = sp.orbit(a, N, r)
    F = ob.obs_comb(orb, N, int(orb[123]))
    P, half, norm = sp.spectrum(F)
    assert np.allclose(P, 1.0 / r)


def test_jacobi_exceptional_only():
    p, q = 1009, 2003
    N = p * q
    rng = random.Random(5)
    for _ in range(4):
        a = _unit(rng, N)
        r = sp.order_info(a, p, q)["r"]
        orb = sp.orbit(a, N, r)
        F = ob.obs_jacobi(orb, N)
        ja = gmpy2.jacobi(a, N)
        assert np.array_equal(F, np.where((np.arange(r) & 1) & (ja == -1), -1.0, 1.0))
        met = sp.metrics(sp.spectrum(F), r, ja)
        assert met["M_peak"] < 1e-20 and abs(met["exceptional_mass"] - 1) < 1e-9


def test_orbit_shortcuts_match_direct():
    """The harness shortcuts in run_experiments.panel_values equal the direct observables."""
    from src.periodfind.run_experiments import panel_values
    p, q = 1009, 2003
    N = p * q
    rng = random.Random(9)
    for _ in range(3):
        a = _unit(rng, N)
        r = sp.order_info(a, p, q)["r"]
        orb = sp.orbit(a, N, r)
        for name, direct in [("P5_jacobi", ob.obs_jacobi), ("P7_jacobi_interval", ob.obs_jacobi_interval),
                             ("P9_lsb_inverse", ob.obs_lsb_inverse), ("P11_oa", ob.obs_oa)]:
            assert np.array_equal(panel_values(name, a, N, r, orb, rng), direct(orb, N)), name


def test_density_bound():
    """Any 0/1 observable of orbit density rho has nonzero relative peaks <= rho/(1-rho)."""
    rng = np.random.default_rng(0)
    r = 4096
    for rho in (0.01, 0.1, 0.5):
        F = (rng.random(r) < rho).astype(float)
        F -= F.mean()
        P, _, _ = sp.spectrum(F)
        dens = (F > 0).mean()
        assert P[1:].max() <= dens / (1 - dens) + 1e-12


def test_generic_ring_theorem_random_trees():
    """Theorem GR-F: visible mass <= 8 eps_split for random zero-test decision trees."""
    p, q = 433, 599
    N = p * q
    rng = random.Random(11)
    a = 5
    info = sp.order_info(a, p, q)
    r, rp, rq, g = info["r"], info["r_p"], info["r_q"], info["g"]
    orb = sp.orbit(a, N, r)
    exps = [rp, rp // 2, rq, rq // 3 if rq % 3 == 0 else rq, r // 2, rng.randrange(2, N), 6, 12]
    nonconst = 0
    for _ in range(30):
        tree = gen.random_tree(3, rng, N, exps)
        res = gen.check_tree(tree, orb, N, p, q, r, g)
        assert res["ok"], res
        nonconst += res["nonconstant"]
    assert nonconst > 0


def test_sft_finds_planted_tone_and_recovers_order():
    p, q = 1009, 2003
    N = p * q
    a = 3
    r = sp.order_info(a, p, q)["r"]
    m = 2 * N.bit_length() + 2
    rr = np.uint64(r)
    tab = np.exp(2j * np.pi * np.arange(r) / r)
    oracle = CountingOracle(lambda x: tab[(x % rr).astype(np.int64)])
    res = sft(oracle, m, 0.25, seed=1)
    assert recover_order(res["out"], 1 << m, a, N) == r


def test_sft_on_comb_finds_nothing_useful():
    p, q = 1009, 2003
    N = p * q
    a = 3
    r = sp.order_info(a, p, q)["r"]
    orb = sp.orbit(a, N, r)
    tab = ob.obs_comb(orb, N, int(orb[7]))
    rr = np.uint64(r)
    oracle = CountingOracle(lambda x: tab[(x % rr).astype(np.int64)])
    m = 2 * N.bit_length() + 2
    res = sft(oracle, m, 0.25, seed=2)
    assert recover_order(res["out"], 1 << m, a, N) != r


def test_cf_denominators():
    Q = 1 << 40
    r = 1000003
    k = round(Q * 17 / r)
    assert r in cf_denominators(k, Q, 10 ** 7)


def test_generic_ring_nonsplit_nonconstant_branch():
    """A zero test y^{r/d} = 1 with d | g, gcd(d, r_p r_q / g^2) = 1 is non-split and non-constant:
    all its mass sits on (r/g)Z, consistent with GR-F (visible mass 0 = 8 * eps_split)."""
    rng = random.Random(21)
    found = 0
    for p, q in [(1009, 2003), (433, 599), (877, 1013), (1201, 1489)]:
        N = p * q
        for _ in range(40):
            a = _unit(rng, N)
            info = sp.order_info(a, p, q)
            r, rp, rq, g = info["r"], info["r_p"], info["r_q"], info["g"]
            al, be = rp // g, rq // g
            ds = [d for d in range(2, g + 1) if g % d == 0 and math.gcd(d, al * be) == 1]
            if not ds:
                continue
            d = ds[0]
            tree = gen.Node(m=r // d, c=1, yes=gen.Node(leaf=1.0 + 0j), no=gen.Node(leaf=-1.0 + 0j))
            res = gen.check_tree(tree, sp.orbit(a, N, r), N, p, q, r, g)
            assert res["eps_split"] == 0 and res["nonconstant"] and res["visible_mass"] < 1e-20
            found += 1
            break
    assert found >= 1
