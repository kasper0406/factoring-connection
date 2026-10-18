import random
import pathlib

from src.nfsfruit.gf2 import left_kernel, coords, rank_of
from src.nfsfruit.cubic import Order
from src.nfsfruit.smooth import FactorBase
from src.nfsfruit.harness import run_instance, phi_class
from src.nfsfruit import families as fam
from src.nfsfruit.instances import composite


def test_gf2_kernel_random():
    rng = random.Random(1)
    for _ in range(30):
        M, C = rng.randint(1, 40), rng.randint(1, 30)
        rows = [rng.getrandbits(C) for _ in range(M)]
        rk, kern = left_kernel(rows)
        assert rk + len(kern) == M
        assert rank_of(rows) == rk
        for S in kern:
            acc, i, T = 0, 0, S
            while T:
                if T & 1:
                    acc ^= rows[i]
                T >>= 1
                i += 1
            assert acc == 0
        if kern:
            c = rng.getrandbits(len(kern)) or 1
            S = 0
            for j in range(len(kern)):
                if (c >> j) & 1:
                    S ^= kern[j]
            assert coords(S, kern) == c


def test_smooth():
    fb = FactorBase(50)
    assert fb.factor(2 ** 5 * 3 * 47) == {2: 5, 3: 1, 47: 1}
    assert fb.factor(53 * 2) is None
    assert fb.factor(-12) == {2: 2, 3: 1}


def test_cubic_sqrt_roundtrip():
    rng = random.Random(2)
    O = Order([7, -3, 2, 1])            # x^3 + 2x^2 - 3x + 7
    ls = O.inert_primes(1 << 20, 2)
    for _ in range(20):
        v = [rng.randint(-10 ** 30, 10 ** 30) for _ in range(3)]
        g = O.mul(v, v)
        w = O.sqrt(g, ls[0])
        assert w is not None and (w == v or w == [-x for x in v])
    # a non-square: alpha itself (norm -7 is not a square)
    assert O.sqrt([0, 1, 0], ls[0]) is None


def test_phi_class_semiprime():
    p, q = 7, 11
    N = p * q
    # X/Y = 1 mod 7, -1 mod 11 -> class 1
    r = 43  # 43 = 1 mod 7, 43 = -1 mod 11
    assert phi_class(r, 1, N, [p, q]) == 1
    assert phi_class(1, 1, N, [p, q]) == 0
    assert phi_class(N - 1, 1, N, [p, q]) == 0


def test_controls_and_lemma_h_small():
    cfg = dict(B=300, excess=30, max_trials=400000)
    r = run_instance("trivial", "blum", 40, 123, cfg)
    assert r["status"] == "ok" and r["phi_rank"] == 0
    r = run_instance("dixon", "mixed", 40, 123, cfg)
    assert r["status"] == "ok" and r["lin_bad"] == 0 and r["phi_rank"] == 1
    r = run_instance("lin_rand", "tri", 40, 123, dict(B=1500, A=600, Bb=200, excess=30))
    assert r["status"] == "ok" and r["lin_bad"] == 0 and r["phi_rank"] == 2


def test_nfs3_congruences_exact():
    cfg = dict(B=1500, A=1000, Bb=400, nchar=40, excess=20)
    r = run_instance("nfs3_base", "both1", 40, 5, cfg, subset_trials=5, lin_checks=3)
    assert r["status"] == "ok" and r["nonsquare_basis"] == 0 and r["lin_bad"] == 0


def test_no_factor_leakage_in_operational_modules():
    here = pathlib.Path(__file__).resolve().parents[1] / "src" / "nfsfruit"
    for name in ("families.py", "gf2.py", "smooth.py", "cubic.py"):
        src = (here / name).read_text()
        assert "instances" not in src and "harness" not in src
        assert "factorint" not in src and ".factor(" not in src.replace("fb.factor(", "").replace(
            "nmod_poly(self.f, l).factor(", "").replace("nmod_poly(self.f, q).factor(", "").replace(
            "fmpz_poly(f).factor(", "")
