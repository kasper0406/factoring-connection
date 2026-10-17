from math import comb

from src.circuit_sections import build_power, run_path, section_once
from src.common import instances as I
from src.seed_checks import kummer_support, factorize_reference
from src.trackA.diagnostic import section_coeffs_explicit, section_support, support_explicit
from src.trackA.sections import descent_split, section_values, section_values_ext
from src.trackA.selectors import ExactScalarSelector, descend
from src.trackA.theory_checks import berlekamp_massey, min_lc_masked
from src.trackA.zero_tests import scalar_all_zero, sections_of


def F(N):
    return [0] + [comb(N, k) % N for k in range(1, N)] + [0]


def test_circuit_sections_exact():
    for N in (143, 561, 1001):
        Fc = F(N)
        for kar in (False, True):
            C, o = build_power(N, subtract_endpoints=True)
            r = 0
            for k, b in enumerate([1, 0, 1, 1, 0]):
                C, o = section_once(C, o, b, kar)
                r += b << k
                want = [Fc[i] for i in range(r, N + 1, 2 ** (k + 1))]
                while len(want) > 1 and want[-1] == 0:
                    want.pop()
                assert C.expand(o) == want


def test_circuit_growth_is_superlinear():
    rows = run_path(40009 * 40013, [0] * 6)
    n = [x["nodes"] for x in rows]
    assert all(n[i + 1] > 3 * n[i] for i in range(2, 6))


def test_section_values_match_explicit():
    for N in (143, 899, 1001):
        for M in (1, 4, 32, 2048):
            for c in (0, 3, N - 2):
                v = section_values(N, M, c)
                for r in range(min(M, 64)):
                    S = section_coeffs_explicit(N, r, M)
                    assert v[r] == sum(a * pow(c, j, N) for j, a in S.items()) % N
        T = section_values_ext(N, 8, 3, 2)
        for r in range(8):
            S = section_coeffs_explicit(N, r, 8)
            w = [0] * 3
            for j, x in S.items():
                w[j % 3] = (w[j % 3] + x * pow(2, j // 3, N)) % N
            assert T[r] == w
        e, o = descent_split(N, 16, 5)
        assert e + o == section_values(N, 32, 25)
    import pytest
    from src.common.ring import NonUnit
    with pytest.raises(NonUnit) as ei:
        descent_split(1001, 4, 7)  # 2s = 14 shares 7 with 1001
    assert ei.value.divisor == 7


def test_diagnostic_support():
    for N in range(15, 700, 2):
        f = factorize_reference(N)
        if len(f) >= 2 and all(x == 1 for x in f.values()):
            fs = sorted(f)
            sup = kummer_support(N)
            assert support_explicit(N, fs) == sorted(sup)
            for M in (4, 64):
                for r in range(M):
                    assert section_support(N, fs, r, M) == sorted(i for i in sup if i % M == r)


def test_props_C1_C2_small():
    for seed in range(3):
        inst = I.balanced_semiprime(20, seed)
        p, q = inst.factors
        sup = support_explicit(inst.N, [p, q])
        assert len(sup) == 2 * q - p - 1
        M = 1
        while M <= inst.N:
            cls = {}
            for i in sup:
                cls[i % M] = cls.get(i % M, 0) + 1
            if M < p:
                assert len(cls) == M
            if M > p:
                assert max(cls.values()) <= 3
            M *= 2


def test_scalar_false_zero_none_for_balanced():
    for N, pr in ((143, [11, 13]), (323, [17, 19]), (899, [29, 31])):
        co = [comb(N, i) % N for i in range(N + 1)]
        M = 1
        while M <= 2 * N:
            for r, sec in sections_of(N, co, M).items():
                assert not scalar_all_zero(sec, pr)
            M *= 2
    co = [comb(561, i) % 561 for i in range(562)]
    assert scalar_all_zero(sections_of(561, co, 1)[0], [3, 11, 17])


def test_exact_selector_reaches_support_16bit():
    inst = I.balanced_semiprime(16, 2)
    r, M, logs, disc = descend(inst.N, ExactScalarSelector(), seed=0)
    assert M > inst.N and 0 < r < inst.N and 1 < __import__("math").gcd(r, inst.N) < inst.N


def test_gap_lemma_small():
    assert berlekamp_massey([0, 0, 0, 1], 2) == 4
    assert min_lc_masked([int(i % 4 == 1) for i in range(12)], 2, 4) == 4
