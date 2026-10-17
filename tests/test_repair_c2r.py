"""Tests for the C-2R repair cycle (agent `repair`).  Small primes only."""
import re
from math import gcd
from pathlib import Path

import numpy as np

from src.repair import c2r_mixing as X
from src.repair.c2r_search import (brute_fold, criterion_fold, instance_primes, is_prime,
                                   params, search_prime, squares_mask, verify_certificate)

ROOT = Path(__file__).resolve().parents[1]


def test_instance_list_matches_audit():
    P = instance_primes()
    assert len(P) == 32 and P[:4] == [257, 263, 269, 271] and P[-1] == 24631
    assert all(is_prime(p) for p in P)


def test_reflection_criterion_equals_brute_force():
    for p in [11, 13, 19, 23, 29, 31, 43, 47]:
        m, _ = params(p)
        H = squares_mask(p)
        for E in range(2, p - 1, 2):
            if gcd(E, p - 1) != 2:
                continue
            for v in range(p):
                assert brute_fold(p, E, v) == criterion_fold(p, E, (-2 * v) % p, m, H)


def test_pair_sums_cover_all_translations_and_gcd_restriction():
    for p in [23, 29, 31, 37, 41, 43, 47, 53]:
        m, _ = params(p)
        best = m + 1
        for E in range(1, p - 1):
            for v in range(p):
                ex, d = brute_fold(p, E, v)
                if ex and gcd(E, p - 1) != 2:
                    assert d >= m + 1  # gcd 1 cannot beat ceil(L/2); gcd >= 3 is never exact
                if ex and gcd(E, p - 1) == 2:
                    best = min(best, d)
        assert search_prime(p)["best_d"][1] == best


def test_cross_check_against_audit_code():
    src = (ROOT / "notes" / "codex_r3_critic.md").read_text()
    code = re.search(r"```python\n(.*?)```", src, re.S).group(1).split("if __name__")[0]
    ns = {}
    exec(code, ns)  # audit's code, used only as an oracle
    for p in [101, 103, 107, 131, 151, 191, 199, 211, 227, 239]:
        a, b = ns["search"](p), search_prime(p)
        assert a["best"] == b["best_d"] and a["passes"] == b["passes"]


def test_certificate_p1031():
    r = verify_certificate()
    assert r["zero_set_is_Z_L"] and r["roots_equal_claim"] and r["E_is_(p+1)/2"]
    assert r["L"] == 33 and r["deg_Q"] == 11


def test_mixing_stage1_degree_matches_direct_evaluation():
    p = 103
    m, L = params(p)
    rng = np.random.default_rng(0)
    for _ in range(30):
        E1, E2 = sorted(rng.choice(np.arange(2, 65, 2), 2, replace=False).tolist())
        beta = int(rng.integers(1, p))
        W = np.array([(pow(j, E1, p) + beta * pow(j, E2, p)) % p for j in range(m + 1)])
        for d, c in X.exact_d_rowwise(W, p):
            if c is None:
                continue
            v = (-c * pow(2, -1, p)) % p
            vals = {pow((pow((x + m) % p, E1, p) + beta * pow((x + m) % p, E2, p) + v) % p, 2, p)
                    for x in [(-i) % p for i in range(L)]}
            assert len(vals) == d
