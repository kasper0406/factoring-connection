import json
import random
from math import gcd

from src.lemmaM import descent as dsc
from src.lemmaM.harness_e1 import run_instance, crt_eps
from src.lemmaM.unit_twist import cf_period, classify_unit
from src.nfsfruit.instances import composite


def test_gauss_short_is_in_lattice_and_short():
    rng = random.Random(1)
    for _ in range(200):
        N = rng.randrange(10**8, 10**9) | 1
        R = rng.randrange(1, N)
        u, v = dsc.gauss_short(R, N)
        assert (u - R * v) % N == 0
        assert u * u + v * v <= 2 * N  # Minkowski-type bound for the shortest vector


def test_tonelli():
    for p in [7, 13, 17, 97, 1009, 7681, 65537]:
        for n in range(1, 60):
            if pow(n, (p - 1) // 2, p) == 1:
                r = dsc.tonelli(n, p)
                assert r * r % p == n % p


def test_factor_cofactor():
    assert dsc.factor_cofactor(101 * 103, 200, 3) == {101: 1, 103: 1}
    assert dsc.factor_cofactor(101 * 1009, 200, 3) is None
    assert dsc.factor_cofactor(101 ** 4, 200, 3) is None


def test_descent_congruences_exact_and_flip():
    cfg = {"excess": 8, "n_mixD": 10, "max_trials": 20000, "kernel_samples": 15,
           "plain_n": 5, "36": {"B": 300, "B1": 90000, "J": 400, "qs_M": 200000}}
    for cls in ["blum", "mixed", "tri"]:
        for mode in ["pure", "mixed"]:
            rec = run_instance(cls, 36, 3, cfg, mode)
            if rec.get("accidental_factor"):
                continue
            res = rec["res"]
            assert res["nonsquare"] == 0
            assert res["flip_fail"] == 0
            assert res["flip_ok"] == res["D"][0]


def test_descent_output_depends_on_R_only():
    N, fs = composite("blum", 40, 5)
    sysm = dsc.DescentSystem(N, 400, 160000, 500)
    eps = crt_eps(fs, N)
    rng = random.Random(7)
    hits = 0
    for _ in range(300):
        r = rng.randrange(2, N - 1)
        R = r * r % N
        a = sysm.split_and_descend(R)
        b = sysm.split_and_descend((eps * r) ** 2 % N)
        assert a == b
        if a is not None:
            W, E = a
            # (r W)^2 equals the signed rational T
            T = 1
            for p, e in E.items():
                if p == -1:
                    T = -T if e % 2 else T
                else:
                    T = T * pow(p, e, N) % N
            assert (r * W) ** 2 % N == T % N
            hits += 1
    assert hits > 10


def test_unit_twist_identities():
    for cls in ["blum", "mixed", "both1"]:
        for s in range(3):
            N, fs = composite(cls, 28, s)
            r = cf_period(N, 10**6)
            assert r["ok_Q"] and r["ok_A"] and r["ok_unit"]
            c = classify_unit(N, r["a"], r["sign"], fs)
            assert c in ("nontrivial", "trivial", "norm-1")
            if c == "nontrivial":
                assert (r["a"] * r["a"] - 1) % N == 0
