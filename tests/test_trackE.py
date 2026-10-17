"""Tests for Track E (partial-information formats -> lattice decoders)."""
import itertools
import math
import random
import re
from math import gcd
from pathlib import Path

import pytest

from src.common.instances import balanced_semiprime
from src.trackE import oracles
from src.trackE.coppersmith import (AffineHint, balanced_promise, best_m, choose_params,
                                    decode_affine, decode_interval, decode_low_bits,
                                    decode_residue, rigorous_ok)
from src.trackE.formats import (crt, decode_N_mod_pm1, decode_random_bits, decode_residues,
                                decode_residues_noisy, decode_sum_mod_2k, decode_sum_mod_M,
                                hs_branch, labels_resolve, legendre, quad_root_cosets,
                                quad_roots_mod_prime, small_primes, sqrt_mod_prime)
from src.trackE.info_floor import census

ROOT = Path(__file__).resolve().parents[1]


def _fac(g, N):
    return g is not None and 1 < g < N and N % g == 0


# ---------------------------------------------------------------- leakage audit
@pytest.mark.parametrize("mod", ["coppersmith.py", "formats.py"])
def test_operational_modules_do_not_import_harness(mod):
    src = (ROOT / "src" / "trackE" / mod).read_text()
    code = "\n".join(l for l in src.splitlines() if not l.strip().startswith("#"))
    code = re.sub(r'""".*?"""', "", code, flags=re.S)
    for bad in ("oracles", "instances", "inst.", ".labels[", "factorint", ".factor(", "is_prime"):
        assert bad not in code, (mod, bad)


# ---------------------------------------------------------------- (a) / (b)
@pytest.mark.parametrize("n", [64, 128, 256])
def test_low_bits_decoder(n):
    for seed in range(3):
        inst = balanced_semiprime(n, seed)
        k = n // 4 + 10
        h = oracles.oracle_low_bits(inst, k)
        g = decode_low_bits(inst.N, h["k"], h["r"])
        assert _fac(g, inst.N)


def test_wrong_residue_fails_cleanly():
    inst = balanced_semiprime(128, 0)
    k = 50
    r = oracles.oracle_low_bits(inst, k, flip_positions=[3])["r"]
    assert decode_low_bits(inst.N, k, r, d=11) is None


def test_arbitrary_modulus_and_interval():
    inst = balanced_semiprime(128, 1)
    N = inst.N
    M = 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29 * 31 * 37 * 41    # ~ 2^47
    assert _fac(decode_residue(N, M, oracles.oracle_residue(inst, M)["r"]), N)
    hb = oracles.oracle_high_bits(inst, 40)
    assert _fac(decode_interval(N, hb["A"], hb["W"]), N)


def test_nonunit_modulus_is_returned_as_factor():
    inst = balanced_semiprime(64, 2)
    p = inst.labels["p"]
    g = decode_affine(inst.N, AffineHint(p * 3, 5, *balanced_promise(inst.N)))
    assert g == p


def test_rigorous_condition_implies_success():
    # whenever the exact HG+LLL inequality holds, decoding must succeed (theorem)
    for seed in range(3):
        inst = balanced_semiprime(128, seed)
        N = inst.N
        L, U = balanced_promise(N)
        for k in (40, 44, 50):
            hint = AffineHint(1 << k, inst.labels["p"] % (1 << k), L, U)
            lo, hi = hint.x_range(); xm = (lo + hi) // 2; Y = max(hi - xm, xm - lo)
            pr = choose_params(N, L, Y, dmax=40)
            assert pr is not None
            assert rigorous_ok(N, L, Y, *pr)
            assert _fac(decode_affine(N, hint, d=pr[0], m=pr[1]), N)


# ---------------------------------------------------------------- format (i)
def test_quad_root_cosets_bruteforce_and_v_formula():
    rng = random.Random(7)
    for _ in range(200):
        p = rng.randrange(3, 3000) | 1
        q = rng.randrange(3, 3000) | 1
        if p == q:
            continue
        N, k = p * q, rng.randrange(1, 12)
        s = (p + q) % (1 << k)
        brute = {x for x in range(1 << k) if (x * x - s * x + N) % (1 << k) == 0}
        cos = quad_root_cosets(s, N, 2, k)
        got = {x for a, e in cos for x in range(a % (1 << e), 1 << k, 1 << e)}
        assert got == brute
        v = ((p - q) & -(p - q)).bit_length() - 1
        if k > 2 * v:
            assert len(brute) == 2 ** (v + 1) and len(cos) == 2
            assert {e for _, e in cos} == {k - v}
            assert {a % (1 << (k - v)) for a, _ in cos} == {p % (1 << (k - v)), q % (1 << (k - v))}
        else:
            assert len(cos) == 1 and cos[0][1] == (k + 1) // 2


def test_sum_mod_2k_decoder():
    inst = balanced_semiprime(128, 3)
    k = 48
    o = oracles.oracle_sum_mod(inst, 1 << k)
    assert _fac(decode_sum_mod_2k(inst.N, o["s"], k), inst.N)


# ---------------------------------------------------------------- format (ii)/(iv)
def test_sqrt_and_roots_mod_prime():
    for l in (3, 5, 13, 17, 97, 1009, 65537):
        for a in range(1, min(l, 200)):
            r = sqrt_mod_prime(a, l)
            assert (r is None) == (legendre(a, l) == -1)
            if r is not None:
                assert r * r % l == a % l
    p, q = 1000003, 1000033
    for l in (101, 103, 107):
        assert set(quad_roots_mod_prime((p + q) % l, p * q, l)) == {p % l, q % l}


def test_legendre_labels_resolve_iff_N_nonresidue():
    # exhaustive: the two roots u, c/u have different Legendre symbols iff (c/l) = -1
    for l in (3, 5, 7, 11, 13, 17, 19, 23):
        for c in range(1, l):
            for u in range(1, l):
                v = c * pow(u, -1, l) % l
                if u == v:
                    continue
                assert (legendre(u, l) != legendre(v, l)) == (legendre(c, l) == -1)


def test_sum_mod_M_list_and_labels():
    inst = balanced_semiprime(64, 0)
    N = inst.N
    ps = [l for l in small_primes(200) if l > 2 and N % l][:8]
    o = oracles.oracle_sum_mod(inst, math.prod(ps))
    s_res = {l: o["s"] % l for l in ps}
    from src.common.metrics import Metrics
    mt = Metrics()
    g = decode_sum_mod_M(N, s_res, metrics=mt)
    assert _fac(g, N)
    amb = mt.fields["omega_ambiguous"]
    assert mt.fields["list_size"] == 2 ** amb
    # with Legendre labels: list shrinks by 2 for every l with (N/l) = -1
    lab = dict(zip(ps, oracles.oracle_legendre(inst, ps)["legendre"]))
    mt2 = Metrics()
    assert _fac(decode_sum_mod_M(N, s_res, labels=lab, metrics=mt2), N)
    n_res = sum(1 for l in ps if labels_resolve(N, l) and (inst.labels["p"] - N // inst.labels["p"]) % l)
    assert mt2.fields["list_size"] == 2 ** (amb - n_res)


# ---------------------------------------------------------------- format (iii)
def test_residues_and_noisy_residues():
    inst = balanced_semiprime(128, 2)
    N = inst.N
    ps = [l for l in small_primes(1000) if l > 2 and N % l][:40]
    o = oracles.oracle_residues(inst, ps[:14])
    assert _fac(decode_residues(N, o["primes"], o["residues"]), N)
    o2 = oracles.oracle_residues(inst, ps, wrong={0, 5, 9}, rng=random.Random(1))
    assert _fac(decode_residues_noisy(N, o2["primes"], o2["residues"], need_bits=44, trials=300,
                                      seed=0), N)


def test_crt():
    r, M = crt([2, 3, 4], [3, 5, 7])
    assert M == 105 and r % 3 == 2 and r % 5 == 3 and r % 7 == 4


# ---------------------------------------------------------------- (v), (vi)
def test_N_mod_pm1_pollard():
    for n in (64, 128):
        for seed in range(4):
            inst = balanced_semiprime(n, seed)
            assert _fac(decode_N_mod_pm1(inst.N, oracles.oracle_N_mod_pm1(inst)["e"]), inst.N)


def test_hs_branch_contains_truth_and_decodes():
    inst = balanced_semiprime(128, 4)
    N = inst.N
    p = inst.labels["p"]; q = N // p
    rng = random.Random(3)
    o = oracles.oracle_random_bits(inst, 0.6, 46, rng)
    cands = hs_branch(N, 46, o["known_p"], o["known_q"])
    assert (p % (1 << 46), q % (1 << 46)) in cands
    assert _fac(decode_random_bits(N, 46, o["known_p"], o["known_q"], depth=46), N)


# ---------------------------------------------------------------- information floor
def test_info_floor_census_prime():
    for l in (7, 11, 13):
        rows = {r["quantity"]: r for r in census(l)}
        assert rows["N mod m"]["determined_by_N_mod_m"]
        assert rows["(N/m) Legendre"]["determined_by_N_mod_m"]
        assert abs(rows["p mod m"]["H_min"] - math.log2(l - 1)) < 1e-3
        # symmetric quantities are functions of s given c (prime modulus)
        for q, r in rows.items():
            if r["symmetric"]:
                assert r["function_of_s_given_c"], q
        assert not rows["(p/m) Legendre"]["symmetric"]


def test_number_of_unordered_pairs():
    for l in (5, 7, 11, 13, 17):
        for c in range(1, l):
            pairs = {frozenset((u, c * pow(u, -1, l) % l)) for u in range(1, l)}
            assert len(pairs) == (l + legendre(c, l)) // 2
