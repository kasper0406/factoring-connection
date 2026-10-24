"""Tests for agent `geometry` (phase 26)."""
import ast
import os

from src.common.instances import balanced_semiprime
from src.geometry import modcurve as mc
from src.geometry import curves as cv
from src.geometry import harness as H
from src.geometry import decoders as D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_decoders_import_no_factoring():
    src = open(os.path.join(ROOT, "src", "geometry", "decoders.py")).read()
    tree = ast.parse(src)
    mods = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods |= {a.name for a in node.names}
        elif isinstance(node, ast.ImportFrom):
            mods.add(node.module or "")
    assert mods <= {"fractions", "itertools", "math", "typing"}, mods
    for bad in ("factorint", "isprime", "is_prime", "harness", "modcurve", "sympy", "gmpy2", "flint"):
        assert bad not in src


def test_orbifold_formulas_match_p1_action():
    for N in list(range(2, 120)) + [143, 221, 385, 1001]:
        f = mc.factor_small(N)
        ob = mc.orbifold_from_p1(N)
        assert ob == {"mu": mc.psi(f), "nu2": mc.nu2(f), "nu3": mc.nu3(f), "nu_inf": mc.nu_inf(f),
                      "genus": mc.genus(f)}


def test_trace_formula_level_one_and_11a1():
    for l in (2, 3, 5, 7, 11, 13, 97):
        assert mc.trace_formula({}, l) == 0
    for l in (2, 3, 5, 7, 13, 17, 19, 23):
        assert mc.trace_formula({11: 1}, l) == cv.ap(cv.CREMONA["11a1"], l)


def test_manin_symbols_small():
    for N, l in ((11, 2), (11, 3), (22, 3), (15, 2), (35, 3), (26, 5)):
        f = mc.factor_small(N)
        ms = mc.manin_trace(N, l)
        assert ms["dim"] == 2 * mc.genus(f) + 2 ** len(f) - 1
        assert ms["trace"] == 2 * mc.trace_formula(f, l) + (l + 1) * (2 ** len(f) - 1)


def test_trace_formula_square_index():
    for l in (2, 3, 5):
        assert mc.trace_formula_n({}, l * l) == 0
        assert mc.trace_formula_n({11: 1}, l * l) == cv.ap(cv.CREMONA["11a1"], l) ** 2 - l
        assert mc.x0_count_sq({11: 1}, l) == cv.count_points_fl2(cv.CREMONA["11a1"], l)


def test_reductions_split():
    for seed in range(5):
        inst = balanced_semiprime(48, 1000 + seed)
        N = inst.N
        p, q = inst.labels["factors"]
        for res in (D.from_genus(N, H.genus_oracle(p, q)), D.from_index(N, H.index_oracle(p, q)),
                    D.from_pillai(N, H.pillai_oracle(p, q)), D.from_sigma(N, H.sigma_oracle(p, q)),
                    D.from_necklace(N, lambda k, r: H.necklace_oracle(p, q, k, r))[0],
                    D.from_catmap(N, lambda t, r: H.catmap_oracle(p, q, t, r))[0],
                    D.from_x0_count_sq(N, 2, H.x0_count_sq_oracle(p, q, 2))[0]):
            assert res == (p, q)


def test_x0_fl_counts_give_quadratic_symbols_only():
    done = 0
    seed = 2000
    while done < 5:
        inst = balanced_semiprime(48, seed)
        seed += 1
        if inst.N % 4 != 1:
            continue
        p, q = inst.labels["factors"]
        counts = {l: H.x0_count_oracle(p, q, l) for l in (2, 3, 5, 7, 11, 13)}
        eps, ncons, unk = D.decode_eps_from_x0_counts(inst.N, counts)
        assert eps == D.kron(-4, p)
        done += 1


def test_trace_formula_brute_vs_product_branch():
    for N in (5 * 7, 7 * 11, 11 * 13, 13 * 17, 7 * 11 * 13, 17 * 19):
        f = mc.factor_small(N)
        for n in (2, 3, 4, 9):
            if N % n == 0 or min(f) ** 2 <= 4 * n:
                continue
            assert mc.trace_formula_n(f, n, brute=True) == mc.trace_formula_n(f, n, brute=False)
