"""Tests for agent `classE` (phase 25)."""
import ast
import os
import random
from fractions import Fraction
from math import gcd, isqrt

import gmpy2

from src.classE import diag
from src.classE.fd import fd_C, fd_exact, fd_float, fd_scaled, floor_sums
from src.classE.structure import (convergents, cf, path_gcd_hit, remainders, sqrt_decomposition, two_halves)
from src.newfamily.recip import dedekind6

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _rp(rng, b):
    while True:
        x = rng.getrandbits(b) | (3 << (b - 2)) | 1
        if gmpy2.is_prime(x):
            return x


def _units(rng, N, k):
    out = []
    while len(out) < k:
        x = rng.randrange(1, N)
        if gcd(x, N) == 1:
            out.append(x)
    return out


def test_floor_sums_brute():
    rng = random.Random(1)
    for _ in range(400):
        n, m, a, b = rng.randrange(0, 50), rng.randrange(1, 40), rng.randrange(0, 90), rng.randrange(0, 90)
        f = sum((a * x + b) // m for x in range(n))
        g = sum(x * ((a * x + b) // m) for x in range(n))
        h = sum(((a * x + b) // m) ** 2 for x in range(n))
        assert floor_sums(n, a, b, m) == (f, g, h)


def test_fourier_dedekind_exact_vs_brute_and_float():
    rng = random.Random(2)
    for N in (5, 9, 15, 21, 35, 77, 101):
        for d in (1, 2, 3):
            for _ in range(4):
                a = _units(rng, N, d)
                n = rng.randrange(0, N)
                e = fd_exact(n, a, N)
                assert e == fd_exact(n, a, N, brute=True)
                assert abs(fd_float(n, a, N) - float(e)) < 1e-8


def test_d3_unshifted_is_three_dedekind_sums():
    """Prop FD3-0: s_0(a1,a2,a3;N) = (N-1)/(8N) - (1/2) sum_{i<j} s(a_i a_j^{-1}, N)."""
    rng = random.Random(3)
    for _ in range(150):
        N = rng.randrange(5, 3000)
        a = _units(rng, N, 3)
        ded = lambda x: Fraction(dedekind6(x % N, N), 6 * N)
        rhs = Fraction(N - 1, 8 * N) - Fraction(1, 2) * sum(
            ded(a[i] * pow(a[j], -1, N)) for i in range(3) for j in range(i + 1, 3))
        assert fd_exact(0, a, N) == rhs


def test_aggregate_congruences_d2_d3_shifted():
    """Prop FD-AGG (d = 2 proved; d = 3 verified): C/N^{d-1} mod N is a universal polynomial in m, a_j, a_j^{-1}."""
    rng = random.Random(4)
    for _ in range(120):
        N = rng.randrange(200, 4000) | 1
        if gcd(N, 6) != 1:
            continue
        a = _units(rng, N, 3)
        m = rng.randrange(0, N)
        inv = lambda x: pow(x, -1, N)
        # d = 2, shift 0: C/N = -(x + 1/x)/12, x = a1/a2
        C2 = fd_C(0, a[:2], N)
        x = a[0] * inv(a[1]) % N
        assert C2 % N == 0 and ((C2 // N) * 12 + x + inv(x)) % N == 0
        # d = 3 with shift m
        C3 = fd_C(m, a, N)
        assert C3 % (N * N) == 0
        ai = [inv(t) for t in a]
        R = sum(a[i] * ai[j] for i in range(3) for j in range(3) if i != j)
        T1 = sum(a[i] * ai[(i + 1) % 3] * ai[(i + 2) % 3] for i in range(3))
        S1 = sum(ai)
        S2 = ai[0] * ai[1] + ai[0] * ai[2] + ai[1] * ai[2]
        P = ai[0] * ai[1] * ai[2]
        pred = (2 * R - 4 * m * T1 - 12 * m * S1 + 12 * m * m * S2 - 8 * m ** 3 * P) * inv(48) % N
        assert (C3 // (N * N) - pred) % N == 0


def test_cf_identities_and_two_halves():
    rng = random.Random(5)
    for _ in range(300):
        N = rng.randrange(10 ** 5, 10 ** 9)
        a = rng.randrange(1, N)
        if gcd(a, N) != 1:
            continue
        qs = cf(a, N)
        H, S = convergents(qs)
        rs = remainders(a, N) + [0]
        assert S[-1] == N and rs[len(qs)] == 1
        for k in range(len(qs) + 1):
            assert a * S[k] - N * H[k] == (-1) ** k * rs[k + 1]
            if k >= 1:
                assert N == S[k] * rs[k] + S[k - 1] * rs[k + 1]
        th = two_halves(a, N)
        r = th["r"]
        assert (th["astar"] - (-1) ** (r - 1) * pow(a, -1, N)) % N == 0
        rev = cf(th["astar"], N) if th["astar"] > 0 else []
        # mirror formula: CF(astar/N) = [0; a_r, ..., a_1] (a_1 = 1 merges with the next quotient)
        if qs[0] >= 2 and r >= 2:
            assert rev == list(reversed(qs))


def test_sqrt_argument_structure():
    rng = random.Random(6)
    for _ in range(200):
        p, q = sorted((_rp(rng, 30), _rp(rng, 30)))
        if p == q:
            continue
        N = p * q
        d = sqrt_decomposition(N)
        r, sp = d["r"], d["sprime"]
        assert cf(r, N) == [r + d["t"]] + cf(sp, r)
        r2, u, v = diag.sqrt_uv(N, p, q)
        assert (sp + u * v) % r == 0


def test_prop_vis_necessity_and_planted():
    rng = random.Random(7)
    planted_hits = 0
    for t in range(3000):
        p, q = sorted((_rp(rng, 18), _rp(rng, 18)))
        if p == q:
            continue
        N = p * q
        if t % 2:
            a = rng.randrange(2, N)
        else:  # planted a = +-1 (mod q)
            a = diag.planted_unit(p, q, rng.randrange(1, p), rng.choice([1, q - 1]))
        if gcd(a, N) != 1:
            continue
        S = convergents(cf(a, N))[1][:-1]
        hit_p = any(s % p == 0 for s in S)
        hit_q = any(s % q == 0 for s in S)
        assert not hit_q                      # q never divides a proper convergent denominator (q > p)
        if hit_p:
            assert diag.vis_predicate(a, p, q)
            assert path_gcd_hit(a, N) in (p, q)
            planted_hits += t % 2 == 0
    assert planted_hits > 50


def test_gar_identity():
    rng = random.Random(8)
    for _ in range(200):
        p, q = sorted((_rp(rng, 20), _rp(rng, 20)))
        if p == q:
            continue
        x = rng.randrange(0, p * q)
        lhs, rhs = diag.gar_phase(x, p, q)
        assert lhs == rhs


def test_battery_stats_operational_and_deterministic():
    from src.classE.battery import MENU_NAMES, compute_stats
    N = 1000000016000000063  # 10^9+7 times 10^9+9 (public test value)
    a, _ = compute_stats(N)
    b, _ = compute_stats(N)
    assert a == b and set(a) == set(MENU_NAMES)


def test_operational_isolation():
    banned = {"src.common.instances", "src.classE.harness", "src.classE.diag", "src.classE.toy3", "sympy"}
    for fn in ("structure.py", "fd.py", "battery.py"):
        tree = ast.parse(open(os.path.join(ROOT, "src", "classE", fn)).read())
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                assert node.module not in banned, (fn, node.module)
            if isinstance(node, ast.Import):
                for al in node.names:
                    assert al.name not in banned, (fn, al.name)


def test_prop_split_genus_constraint_and_product_size():
    """Prop SPLIT: a s_k = (-1)^k r_{k+1} (mod N) gives ((-1)^k s_k r_{k+1} | l) = (a | l) for l in {p, q};
    and s_k r_{k+1} >= N / (a_{k+1} + 2)."""
    from src.classE.structure import arguments
    rng = random.Random(10)
    for _ in range(150):
        p, q = sorted((_rp(rng, 22), _rp(rng, 22)))
        if p == q:
            continue
        N = p * q
        for a in (arguments(N)["fermres"], rng.randrange(2, N)):
            if gcd(a, N) != 1:
                continue
            qs = cf(a, N)
            H, S = convergents(qs)
            rs = remainders(a, N) + [0]
            for k in range(len(qs)):
                v = (-1) ** k * S[k] * rs[k + 1]
                for l in (p, q):
                    if v % l:
                        assert gmpy2.legendre(v % l, l) == gmpy2.legendre(a % l, l)
                if k < len(qs):          # qs[k] = a_{k+1};  N = s_{k+1} r_{k+1} + s_k r_{k+2} <= (a_{k+1}+2) s_k r_{k+1}
                    assert (qs[k] + 2) * S[k] * rs[k + 1] >= N
