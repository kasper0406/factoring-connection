"""Tests for src/nativeinv (agent nativeinv, phase 29)."""
import ast
import pathlib

import pytest

from src.common.instances import balanced_semiprime
from src.nativeinv import harness as H
from src.nativeinv.trichotomy import QuasiPoly, decode, global_tier, integer_roots, pair_tier


def test_integer_roots_exact():
    assert integer_roots([6, -5, 1], 0, 10) == [2, 3]
    assert integer_roots([-6, 11, -6, 1], -5, 10) == [1, 2, 3]
    assert integer_roots([1, -4, 4], -5, 5) == []          # (2x-1)^2: rational, not integer
    assert integer_roots([0, 0, 1], -3, 3) == [0]
    p, q = 1000003, 1000033
    assert integer_roots([p * q, -(p + q), 1], 2, p * q) == [p, q]


def test_operational_module_is_factor_free():
    src = pathlib.Path("src/nativeinv/trichotomy.py").read_text()
    tree = ast.parse(src)
    banned = {"sympy", "flint", "gmpy2", "harness", "instances"}
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            names = [a.name for a in node.names] + ([node.module] if getattr(node, "module", None) else [])
            for nm in names:
                assert not any(b in (nm or "") for b in banned), nm
        if isinstance(node, ast.Attribute):
            assert node.attr not in {"factorint", "factor", "isprime", "is_prime"}


def test_tiers_of_registry():
    for key, (F, tier) in H.REGISTRY.items():
        assert global_tier(F) == tier, key


@pytest.mark.parametrize("seed", range(5))
def test_decoder_fact_tier_splits(seed):
    inst = balanced_semiprime(64, seed)
    p, q = inst.labels["p"], inst.labels["q"]
    for key, (F, tier) in H.REGISTRY.items():
        V = F.value(p) * F.value(q)
        res = decode(inst.N, V, F)
        if pair_tier(F, p % F.m if F.m > 1 else 0, q % F.m if F.m > 1 else 0) == "FACT":
            assert res.split == (p, q), key
        else:
            assert res.split is None, key


def test_E1_identities_N15_N21():
    for p, q in [(3, 5), (3, 7)]:
        N = p * q
        assert H.k_SL2(N) == (p + 4) * (q + 4) == H.gcdsum_SL2(N)
        assert H.k_Aff(N) == N
        assert H.k_Heis(N) == (p * p + p - 1) * (q * q + q - 1)
        assert H.pillai(N) == H.zero_product_count(N) == (2 * p - 1) * (2 * q - 1)
        degs = H.aff_degree_multiset(p, q)
        assert sum(m * d * d for d, m in degs) == N * (p - 1) * (q - 1)
        assert sum(m for _, m in degs) == N
    assert H.k_GL2(15) == 8 * 24 == H.gcdsum_GL2(15)
    assert H.subgroups_ZN2(15) == 6 * 8
