import math
import random

import numpy as np

from src.sinonobliv.core import (DensityCurve, ch_log_threshold, count_members_in_window, divisors_from_factored,
                                 family_pieces, family_value, is_smooth_int, johnson_log_threshold,
                                 lcm_upto_factored, smooth_mask_values, smooth_mask_window, value_of_factored)
from src.sinonobliv.msample import (_divisors_of, inheritance_scan, joint_tuples_indep, predicted_tuples)
from src.sinonobliv.theory import MNFS, si_cost, tau_threshold, ms_case_bound, prop_cd_log_hits


def test_smooth_masks_agree_with_trial_division():
    rng = random.Random(1)
    vals = [rng.randrange(1, 1 << 40) for _ in range(300)] + [2 ** 39, 3 ** 20, 1]
    m = smooth_mask_values(np.array(vals, dtype=np.uint64), 97)
    assert [bool(b) for b in m] == [is_smooth_int(v, 97) for v in vals]
    lo = (1 << 33) + 12345
    w = smooth_mask_window(lo, 5000, 97)
    assert [bool(b) for b in w] == [is_smooth_int(lo + i, 97) for i in range(5000)]


def test_family_identities():
    for t in range(2, 200):
        a, b = family_pieces("aur", np.array([t], dtype=object))
        assert int(a[0]) * int(b[0]) == 4 * t ** 4 + 1
        assert int(family_value("cp3", np.array([t], dtype=object))[0]) == t * (t + 1) * (t + 2)


def test_count_members_in_window_bruteforce():
    for fam, F in (("cp2", lambda t: t * (t + 1)), ("cp3", lambda t: t * (t + 1) * (t + 2)),
                   ("aur", lambda t: 4 * t ** 4 + 1)):
        for a, b in ((100, 5000), (7, 7), (10 ** 6, 10 ** 6 + 12345)):
            brute = sum(1 for t in range(2, 2000) if a <= F(t) <= b)
            assert count_members_in_window(fam, a, b) == brute


def test_inheritance_lemma_small_exhaustive():
    rng = random.Random(7)
    for _ in range(5):
        N = rng.randrange(1 << 20, 1 << 21) | 1
        res = inheritance_scan(N, [1, 2, 3, 5], 40)
        assert res["violations"] == 0
    # constructed inherited event: d | N + s, s small, gives g > Kh at offsets proportional to k
    N = 1000003 * 7919 - 3  # N + 3 = 1000003 * 7919
    res = inheritance_scan(N, [1, 2, 4], 20)
    assert res["violations"] == 0 and res["big_inherited"] > 0


def test_thresholds():
    assert abs(ch_log_threshold(5.0, 40.0, 1) - johnson_log_threshold(5.0, 40.0)) < 1e-12
    # thresholds decrease with m and tend to t
    vals = [ch_log_threshold(5.0, 40.0, m) for m in range(1, 30)]
    assert all(a > b for a, b in zip(vals, vals[1:])) and vals[-1] > 5.0


def test_predicted_tuples_exact_average_tiny():
    # exhaustive average over all window positions x1, x2 in [X, 2X) for tiny parameters
    X, h, B, m = 600, 6, 8, 2
    fac = lcm_upto_factored(B)
    divs = divisors_from_factored(fac)
    sig = math.log(h + 2)  # d >= h+2 > h+1: at most one multiple per window
    lim = math.exp(sig)
    D = [int(d) for d in divs if int(d) >= lim]
    cnt = {}
    for x in range(X, 2 * X):
        for d in D:
            c = (x + h) // d - (x - 1) // d
            cnt[(x, d)] = c
    tot = 0
    for d in D:
        s = sum(cnt[(x, d)] for x in range(X, 2 * X))
        tot += s * s  # independent x1, x2: (sum_x c)^2
    avg = tot / X ** 2
    pred = predicted_tuples(divs, h, m, sig)
    assert abs(avg / pred - 1) < 0.05


def test_joint_tuples_bruteforce_tiny():
    X, h, B, y, m = 1 << 14, 30, 12, 13, 2
    fac = lcm_upto_factored(B)
    P = value_of_factored(fac)
    sig = math.log(40)
    rng1 = np.random.default_rng(5)
    res = joint_tuples_indep(X, h, P, B, m, sig, y, 300, rng1)
    rng2 = np.random.default_rng(5)
    tot = 0
    for _ in range(300):
        xs = [int(v) for v in rng2.integers(X, 2 * X, size=m)]
        for d in [int(d) for d in divisors_from_factored(fac)]:
            if d < 40:
                continue
            c = 1
            for x in xs:
                c *= (x + h) // d - (x - 1) // d
            tot += c
    assert res["tuples"] == tot


def test_divisors_of():
    assert sorted(_divisors_of(360, [2, 3, 5])) == sorted(d for d in range(1, 361) if 360 % d == 0)


def test_cost_model_constants():
    c, k = si_cost(0.0)
    assert abs(c - (16 / 9) ** (1 / 3)) < 1e-4
    assert abs(tau_threshold() - 0.757) < 2e-3
    assert abs(MNFS - 1.9019) < 1e-4


def test_ms_case_split():
    # case 2 (sigma_m <= ell < j1) forces pi <= (m/(m-1) + 1/t)^{m+1} t <= 8 e^2 t when m <= t
    import random as _r
    rng = _r.Random(3)
    for _ in range(2000):
        t = rng.uniform(3, 50); m = rng.randint(2, int(t)); pi = t * rng.uniform(1, 500)
        r0 = ms_case_bound(t, pi, m, 0.0)
        ell = rng.uniform(r0["sigma"], max(pi, r0["sigma"]))
        if ell > pi:
            continue
        r = ms_case_bound(t, pi, m, ell)
        assert ell >= t - 1e-9
        if r["case"] == 2:
            assert pi <= r["case2_pi_cap"] + 1e-6 and r["case2_pi_cap"] <= 8 * math.e ** 2 * t + 1e-6


def test_prop_cd_sparse_families_vanish():
    x = 1024 * math.log(2)
    # one-parameter, degree-2 family: hopeless; m = d: placement term is O(F_{2/3})
    assert prop_cd_log_hits(x, 1, 2, 2)["total"] < -300
    assert abs(prop_cd_log_hits(x, 3, 3, 3)["placement"]) < 200


def test_hg_certificate_planted():
    from src.sinonobliv.uncert import classify_call
    P = value_of_factored(lcm_upto_factored(40))
    for k, off in ((12345, 5), (99991, 4000), (31337, 0)):
        x = (P // 7) * k - off
        res = classify_call(x, 1 << 12, P, 3, 3, 1 << 10)
        assert res["cert_violations"] == 0 and res["cert_roots"] >= 1
