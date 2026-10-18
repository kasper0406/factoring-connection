"""Tests for src/pm1bsgs: correctness lemmas, operational code vs exact harness
predictors, cost model vs instrumented counts, exact probabilities vs exhaustive
enumeration, and a no-leakage check."""
import random
from fractions import Fraction
from math import gcd

import pytest

from src.common.instances import balanced_semiprime, is_prime
from src.common.metrics import Metrics
from src.pm1bsgs import costmodel as cm
from src.pm1bsgs import harness as H
from src.pm1bsgs import operational as op


def _small_primes(lo, hi):
    return [r for r in range(lo, hi) if is_prime(r)]


# ---------------------------------------------------------------- lemmas
@pytest.mark.parametrize("B", [2, 5, 16, 60])
def test_lemma1_order_divides_R(B):
    """ord_p(a^E) | R = (p-1)/gcd(p-1,E) for every a, exhaustively for p < 400."""
    E = op.stage1_exponent(B)
    for p in _small_primes(3, 400):
        S = gcd(p - 1, E)
        R = (p - 1) // S
        for a in range(1, p):
            assert pow(a, E * R, p) == 1


def test_lemma2_N_congruent_q():
    rng = random.Random(1)
    ps = _small_primes(3, 3000)
    for _ in range(2000):
        p, q = rng.sample(ps, 2)
        N = p * q
        assert (N - q) % (p - 1) == 0          # N = q (mod p-1), hence mod every R | p-1


def test_lemma3_equivalence_with_order():
    """b^N = b^x (mod p)  <=>  ord_p(b) | N - x."""
    rng = random.Random(2)
    ps = _small_primes(3, 500)
    for _ in range(300):
        p, q = rng.sample(ps, 2)
        N = p * q
        b = rng.randrange(1, p)
        o = H.PrimeInfo.of(p).order(b)
        for x in range(0, 3 * p):
            assert (pow(b, N, p) == pow(b, x, p)) == ((N - x) % o == 0)


def test_interval_pure_hit_bruteforce():
    rng = random.Random(3)
    for _ in range(4000):
        dp, dq = rng.randrange(1, 40), rng.randrange(1, 40)
        lo = rng.randrange(0, 200)
        hi = lo + rng.randrange(0, 120)
        brute = any((t % dp == 0) != (t % dq == 0) for t in range(lo, hi + 1))
        assert H.interval_pure_hit(dp, dq, lo, hi) == brute


# ---------------------------------------------------------------- operational vs predictor, cost model
@pytest.mark.parametrize("n", [20, 28, 36])
def test_operational_matches_predictor_and_costmodel(n):
    for s in range(12):
        inst = balanced_semiprime(n, s)
        N = inst.N
        p, q = inst.labels["factors"]
        P, Q = H.PrimeInfo.of(p), H.PrimeInfo.of(q)
        bs = H.base_state(N, P, Q, s)
        for B1 in (1, 8, 64):
            for (mb, J, mode) in [(4, 4, "N"), (8, 16, "order"), (32, 32, "N"), (16, 64, "order")]:
                met = Metrics()
                r = op.pm1_bsgs(N, s, B1, mb, J, mode, met)
                lab = H.outcome_from_orders(bs.op, bs.opf, bs.oq, bs.oqf, N, "bsgs",
                                            {"B1": B1, "m": mb, "J": J, "mode": mode})
                assert (r["factor"] is not None) == (lab in ("stage1", "success2")), (s, B1, mb, J, mode, lab)
                if r["factor"]:
                    assert 1 < r["factor"] < N and N % r["factor"] == 0
                if lab in ("nohit",):
                    model = cm.base_cost() + cm.stage1_cost(B1) + cm.bsgs_stage2_cost(mb, J, mode, N)
                    assert model == met.fields["eq_ops"]
            for B2 in (B1 * 8, 2048):
                if B2 <= B1:
                    continue
                met = Metrics()
                r = op.pm1_prime(N, s, B1, B2, met)
                lab = H.outcome_from_orders(bs.op, bs.opf, bs.oq, bs.oqf, N, "prime", {"B1": B1, "B2": B2})
                assert (r["factor"] is not None) == (lab in ("stage1", "success2"))
                if lab == "nohit":
                    assert cm.base_cost() + cm.stage1_cost(B1) + cm.prime_stage2_cost(B1, B2) == met.fields["eq_ops"]
            st = H.stage1_outcome(bs, B1)
            for K in (8, 64, 512):
                met = Metrics()
                r = op.pm1_birthday(N, s, B1, K, met)
                if st[0] == "ok":
                    fp = H.birthday_first_pure(N, s, st[1], st[2], K)
                    pred = fp is not None
                    if not pred and "backtrack_gcd" not in met.counts:
                        assert cm.base_cost() + cm.stage1_cost(B1) + cm.birthday_stage2_cost(K, N, s) == met.fields["eq_ops"]
                else:
                    pred = st[0] == "stage1"
                assert (r["factor"] is not None) == pred


def test_naive_equals_multipoint():
    for s in range(6):
        inst = balanced_semiprime(28, s)
        a = op.pm1_bsgs(inst.N, s, 32, 16, 16, "N", naive=True)
        b = op.pm1_bsgs(inst.N, s, 32, 16, 16, "N", naive=False)
        assert a["factor"] == b["factor"]


# ---------------------------------------------------------------- exact probabilities vs exhaustive bases
@pytest.mark.parametrize("pq", [(103, 107), (211, 223), (157, 181)])
def test_exact_probs_exhaustive(pq):
    p, q = pq
    N = p * q
    P, Q = H.PrimeInfo.of(p), H.PrimeInfo.of(q)
    for params in ({"B1": 1, "m": 4, "J": 8, "mode": "N"}, {"B1": 4, "m": 4, "J": 4, "mode": "order"},
                   {"B1": 8, "m": 8, "J": 4, "mode": "N"}):
        counts = {}
        units = 0
        for a in range(2, N - 1):
            if gcd(a, N) != 1:
                continue
            units += 1
            bs = H.BaseState(a, P.order(a), Q.order(a), {}, {})
            bs.opf = H.factor_of_order(bs.op, P.fac)
            bs.oqf = H.factor_of_order(bs.oq, Q.fac)
            lab = H.outcome_from_orders(bs.op, bs.opf, bs.oq, bs.oqf, N, "bsgs", params)
            r = op.pm1_bsgs(N, 0, params["B1"], params["m"], params["J"], params["mode"], base=a)
            assert (r["factor"] is not None) == (lab in ("stage1", "success2"))
            counts[lab] = counts.get(lab, 0) + 1
        ex = H.exact_outcome_probs(P, Q, N, "bsgs", params)
        # exhaustive over a in [2, N-2] units = all units minus {1, N-1}
        tot = (p - 1) * (q - 1)
        corr = {}
        for a in (1, N - 1):
            lab = H.outcome_from_orders(P.order(a), H.factor_of_order(P.order(a), P.fac), Q.order(a),
                                        H.factor_of_order(Q.order(a), Q.fac), N, "bsgs", params)
            corr[lab] = corr.get(lab, 0) + 1
        for lab, fr in ex.items():
            assert fr * tot == counts.get(lab, 0) + corr.get(lab, 0)


def test_equal_order_probability():
    p, q = 211, 223
    P, Q = H.PrimeInfo.of(p), H.PrimeInfo.of(q)
    B1 = 2
    E = op.stage1_exponent(B1)
    Rp = (p - 1) // gcd(p - 1, E)
    Rq = (q - 1) // gcd(q - 1, E)
    # b uniform over E-th powers: enumerate b = a^E for a over units mod p and mod q
    from collections import Counter
    cp = Counter(P.order(pow(a, E, p)) for a in range(1, p))
    cq = Counter(Q.order(pow(a, E, q)) for a in range(1, q))
    brute = Fraction(sum(cp[d] * cq[d] for d in cp), (p - 1) * (q - 1))
    assert brute == H.prob_equal_orders(P, Q, B1)


# ---------------------------------------------------------------- leakage
def test_operational_has_no_factor_access():
    import inspect
    for mod in (op, cm):
        src = inspect.getsource(mod)
        for bad in ("sympy", "factorint", "harness", "labels", "is_prime", "n_order", ".factor("):
            assert bad not in src, (mod.__name__, bad)


def test_shift_control_matches_predictor():
    for s in range(10):
        inst = balanced_semiprime(32, s)
        N = inst.N
        p, q = inst.labels["factors"]
        P, Q = H.PrimeInfo.of(p), H.PrimeInfo.of(q)
        bs = H.base_state(N, P, Q, s)
        T = op.shift_target(N, s)
        for B1 in (1, 32):
            st = H.stage1_outcome(bs, B1)
            for mb, J in ((16, 16), (64, 128)):
                r = op.pm1_bsgs(N, s, B1, mb, J, "N", target=T)
                if st[0] == "ok":
                    lo, hi = H.window(mb, J, "N", T)
                    pred = H.interval_pure_hit(st[1], st[2], lo, hi)
                else:
                    pred = st[0] == "stage1"
                assert (r["factor"] is not None) == pred
