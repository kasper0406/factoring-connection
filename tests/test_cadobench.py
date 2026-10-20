"""Tests for src/cadobench (harness code; CADO itself is not needed except where skipped)."""
import math
import os
import random
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from cadobench import memwatch  # noqa: E402
from cadobench.lasbench import parse_las  # noqa: E402
from cadobench.numbers import BITS, NUMBERS, factors_ok  # noqa: E402
from cadobench.relstats import large_ideals, parse_rel  # noqa: E402
from cadobench.run_cado import group_of, parse_log  # noqa: E402
from cadobench.smoothcount import (QLattice, hom_eval, is_smooth, murphy_alpha, primes_upto,  # noqa: E402
                                   product_tree, remainders, rho, smooth_part)
from cadobench.validate import fast_hom  # noqa: E402


def _largest_pf(n):
    m, p, best = n, 2, 1
    while p * p <= m:
        while m % p == 0:
            m //= p
            best = p
        p += 1
    return max(best, m) if m > 1 else best


def test_numbers_digits_and_balance():
    for d, (N, p, q, seed, bits) in NUMBERS.items():
        assert len(str(N)) == d and N.bit_length() == BITS[d] == bits
        assert p * q == N and p < q and q / p < 8 / 3  # top-2-bit primes, bit lengths differ by <= 1
        assert factors_ok(d, [q, p]) and not factors_ok(d, [1, N])


def test_primes_and_trees():
    ps = primes_upto(1000).tolist()
    assert ps[:5] == [2, 3, 5, 7, 11] and len(ps) == 168
    xs = [random.Random(1).randrange(2, 10 ** 12) for _ in range(37)]
    assert product_tree(xs)[-1][0] == math.prod(xs)
    P = math.prod(ps)
    assert [int(r) for r in remainders(P, xs)] == [P % x for x in xs]


def test_smoothness_exact_against_brute_force():
    rng = random.Random(7)
    B = 1 << 10
    P = math.prod(primes_upto(B).tolist())
    xs = [rng.randrange(2, 1 << 30) for _ in range(400)] + [2 ** 29, 1021 * 1019 * 1013, 1031 * 3]
    got = is_smooth(xs, P)
    want = [_largest_pf(x) <= B for x in xs]
    assert got == want
    sp = smooth_part(xs, P)
    for x, s in zip(xs, sp):
        assert x % s == 0 and _largest_pf(x // s) > B if x != s else True


def test_rho_values():
    assert rho(0.5) == 1.0
    assert abs(rho(2) - (1 - math.log(2))) < 1e-6
    assert abs(rho(3) - 0.0486083882911) < 1e-6
    assert abs(rho(4) - 0.00491092564776) / 0.0049109 < 1e-4
    assert abs(rho(5) - 0.000354724700456) / 0.00035472 < 1e-3


def test_alpha_linear_and_cado_value():
    # linear form with coprime coefficients: alpha = sum_p log p / (p^2 - 1) ~ 0.5700 (p <= 2000)
    lin = murphy_alpha([-166188820558869, 3961719737])
    ref = sum(math.log(p) / (p * p - 1) for p in primes_upto(2000).tolist())
    assert abs(lin - ref) < 1e-3 and abs(lin - 0.57) < 0.01  # k-sum truncated at p^k <= 1e6
    # the c60 CADO polynomial: CADO's polyselect prints alpha -4.81
    f1 = [-302524656478155132, 498760510238416, -10081219960, 3104727, 1080]
    assert abs(murphy_alpha(f1) - (-4.81)) < 0.01


def test_norms_match_a_cado_relation():
    # c60 relation printed by las: 317419,39 with side-0 and side-1 factorizations
    f0 = [-166188820558869, 3961719737]
    f1 = [-302524656478155132, 498760510238416, -10081219960, 3104727, 1080]
    line = "317419,39:2,2,2,2,2,7,7,29,83,3871,a7b1:3,3,b,b,b,11,2f,265,287,4a9,bd9b,186a3"
    a, b, p0, p1 = parse_rel(line)
    assert abs(hom_eval(f0, a, b)) == math.prod(p0) == abs(fast_hom(f0, a, b))
    assert abs(hom_eval(f1, a, b)) == math.prod(p1) == abs(fast_hom(f1, a, b))
    li = large_ideals(a, b, p0, p1, 78682, 100000)
    assert li == {(1, 0x186a3, (a * pow(b, -1, 0x186a3)) % 0x186a3)}


def test_qlattice_roundtrip():
    lat = QLattice(q=101999, rho=81287, a0=20712, b0=-1, a1=-1561, b1=5, I=1024, J=512)
    for i, j in [(3, 7), (-500, 1), (511, 511), (1, 0)]:
        a, b = lat.ab(i, j)
        assert (a - lat.rho * b) % lat.q == 0
        assert lat.ij(a, b) == (i, j)


LAS_SNIPPET = """# Sieving side-1 q=101999; rho=81287; a0=20712; b0=-1; a1=-1561; b1=5; I=1024; J=512;
317419,39:2,2,2,2,2,7,7,29,83,3871,a7b1:3,3,b,b,b,11,2f,265,287,4a9,bd9b,186a3
317419,39:2,2,2,2,2,7,7,29,83,3871,a7b1:3,3,b,b,b,11,2f,265,287,4a9,bd9b,186a3
921019,70:6d,15d,239,16f9,f503:2,3,5,11,13,29,9f5,2363,4c4b,e611,186a3
# Average J=459 for 155 special-q's, max bucket fill -bkmult 1,1s:1.330
# Discarded 0 special-q's out of 155 pushed
# Categorized timings (total counted time 0.18):
# 2 sieving on side 0: 0.05
# 5 cofactoring on side 0: 0.01
# survivors before_sieve: 72810496
# survivors enter_cofactoring: 4390
# survivors smooth: 3489
# Total cpu time 1.31s, useful 0.18s [norm 0.00+0.0, sieving 0.1 (0.0+0.1 + 0.0), factor 0.0 (0.0+0.0 + 0.0), rest 0.0], wasted+waited 0.00s, rest 0.00s
# Total elapsed time 0.98s, per special-q 0.0063529, per relation 0.00028223
# PeakMemusage (MB) = 40
# Total 3489 reports [0.000375s/r, 22.5r/sq] in 0.985 elapsed s [132.9% CPU]
"""


def test_parse_las():
    r = parse_las(LAS_SNIPPET, keep_relations=True)
    assert r["relations"] == 3 and r["unique_ab"] == 2
    assert r["n_special_q"] == 155 and r["avg_J"] == 459 and r["discarded_q"] == 0
    assert r["funnel"] == {"before_sieve": 72810496, "enter_cofactoring": 4390, "smooth": 3489}
    assert r["categories"]["2 sieving on side 0"] == 0.05
    assert r["cpu_s"] == 1.31 and r["elapsed_s"] == 0.98 and r["las_peak_mem_mb"] == 40
    assert r["special_q"][0]["q"] == 101999 and r["special_q"][0]["J"] == 512


LOG_SNIPPET = """PID1 2026-10-20 15:40:58,699 Info:Polynomial Selection (size optimized): Starting
PID1 2026-10-20 15:41:00,208 Info:Polynomial Selection (size optimized): Total time: 0.48
PID1 2026-10-20 15:41:00,300 Debug:Command: Process with PID 4 finished successfully
PID1 2026-10-20 15:41:02,554 Info:Lattice Sieving: Starting
PID1 2026-10-20 15:41:12,554 Info:Lattice Sieving: Marking workunit x as ok
PID1 2026-10-20 15:41:13,554 Info:Linear Algebra: Starting
PID1 2026-10-20 15:41:20,554 Info:Square Root: Starting
PID1 2026-10-20 15:41:21,554 Info:Square Root: Factors: 3 5
"""


def test_parse_log(tmp_path):
    p = tmp_path / "c.log"
    p.write_text(LOG_SNIPPET)
    r = parse_log(str(p))
    g = r["wall_by_group_s"]
    assert abs(g["polyselect"] - (1.509 + 2.346)) < 1e-6  # until the next "Starting"
    assert abs(g["sieve"] - 11.0) < 1e-6 and abs(g["linalg"] - 7.0) < 1e-6 and abs(g["sqrt"] - 1.0) < 1e-6
    assert r["factors"] == [3, 5] and abs(r["wall_total_s"] - 22.855) < 1e-6
    assert group_of("Filtering - Merging") == "filter" and group_of("Quadratic Characters") == "linalg"


def test_memwatch_ok_and_kill(tmp_path):
    ok = memwatch.run([sys.executable, "-c", "import time; time.sleep(1.5)"], track=None)
    assert ok["returncode"] == 0 and ok["killed"] is None
    hog = memwatch.run([sys.executable, "-c", "x = bytearray(400 * 2**20); import time; time.sleep(20)"],
                       max_rss_gb=0.2, track=None)
    assert hog["killed"] and "tree RSS" in hog["killed"] and hog["wall_s"] < 15


CADO_LAS = os.path.expanduser("~/src/cado-nfs/build/kmini/sieve/las")
C60 = os.path.expanduser("~/cadobench-work/c60_smoke")


@pytest.mark.skipif(not (os.path.exists(CADO_LAS) and os.path.isdir(C60)), reason="needs the local CADO build")
def test_las_deterministic_yield():
    from cadobench.lasbench import Testbed, run_las
    tb = Testbed.from_workdir("c60_smoke")
    r1 = run_las(tb, 100000, 100300, threads=2, timings=False, keep_relations=True)
    r2 = run_las(tb, 100000, 100300, threads=2, timings=False, keep_relations=True)
    assert r1["relations"] > 0 and sorted(r1["relation_lines"]) == sorted(r2["relation_lines"])
