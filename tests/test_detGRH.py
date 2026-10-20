"""Tests for src/detGRH (Lemma MR, Prop FL, cocycle identity, Z route / Miller)."""
import itertools
from math import gcd

import pytest

from src.detGRH import core, diag


def _brute_sections(N):
    """Independent backtracking count of homomorphic sections of squaring on (Z/N)^{*2}."""
    units = [x for x in range(1, N) if gcd(x, N) == 1]
    sq = sorted({x * x % N for x in units})
    roots = {s: [z for z in units if z * z % N == s] for s in sq}
    gens, span = [], {1}
    for s in sq:
        if s not in span:
            gens.append(s)
            new = set(span)
            while True:
                nxt = {a * s % N for a in new} | new
                if nxt == new:
                    break
                new = nxt
            span = new
    cnt = 0
    for imgs in itertools.product(*[roots[g] for g in gens]):
        m, ok, changed = {1: 1}, True, True
        while changed and ok:
            changed = False
            for x, y in list(m.items()):
                for g, im in zip(gens, imgs):
                    k, v = x * g % N, y * im % N
                    if k in m:
                        if m[k] != v:
                            ok = False
                            break
                    else:
                        m[k] = v
                        changed = True
                if not ok:
                    break
        if ok and all(m[s] * m[s] % N == s for s in sq):
            cnt += 1
    return cnt


def test_lemma_MR_counts_match_independent_brute_force():
    for N in range(3, 160, 2):
        pred = 1 if diag.all_3mod4(N) else 0
        assert diag.count_sections(N) == pred == _brute_sections(N), N


def test_principal_root_multiplicative_blum():
    for N in (21, 33, 57, 77, 3 * 7 * 11):
        assert diag.principal_root_multiplicative(N)


def test_miller_with_phi_factors():
    for N in (77, 65, 5 * 7 * 11 * 13, 41 * 61):
        fs = [f for f in range(2, N) if N % f == 0 and all(f % d for d in range(2, f))]
        phi = 1
        for f in fs:
            phi *= f - 1
        g, a = core.miller(N, phi, 50)
        assert g is not None and 1 < g < N and N % g == 0


def test_f2_basis_express_and_deps():
    rows = [[1, 0, 1, 0], [0, 1, 1, 0], [1, 1, 0, 0], [0, 0, 0, 2], [1, 0, 1, 1]]
    b = core.F2Basis(rows, 4)
    assert len(b.deps) == 2          # rows 0+1+2 and row 3
    for target in range(16):
        combo = b.express(target)
        if combo is not None:
            s = b.combo_sum(combo)
            assert core.bits_of(s) == target


@pytest.mark.parametrize("cls,bits,seed", [("blum", 28, 1), ("both1", 28, 2), ("tri", 32, 3)])
def test_pipeline_identities(cls, bits, seed):
    N, fs = diag.make_instance(cls, bits, seed, min_factor=core.fb_bound(1 << bits))
    fb = core.FactorBase.make(N)
    rel = core.collect(N, fb, K=600, X=fb.B, T=3000)
    # every kernel relation evaluates to 1
    for v in rel.kernel:
        assert fb.ev(v) == 1
    rk, D = core.lattice_index(rel.kernel, fb.n)
    assert rk == fb.n
    Lam = diag.full_lattice(fb, fs)
    phi = 1
    for p in fs:
        phi *= p - 1
    assert diag.lattice_det(Lam) == phi                 # G_B = (Z/N)^*
    assert D % phi == 0                                 # Prop Z: index is a multiple of |G_B|
    assert diag.dim_W(Lam) == len(fs)                   # 2-rank of G_B is omega (Prop FL)
    g, a = core.miller(N, D, core.grh_bound(N))
    assert g is not None and N % g == 0 and 1 < g < N
    c = diag.two_rank_quotient(Lam, rel.kernel)
    _, nfr, _ = core.f2_fruitful(fb, rel.kernel)
    if c < len(fs) - 1:
        assert nfr > 0                                  # Prop FL codimension criterion
    # doubled control: F2 route dead, Z route alive
    _, nfr2, _ = core.f2_fruitful(fb, [[2 * x for x in v] for v in rel.kernel])
    assert nfr2 == 0
    g2, _ = core.miller(N, D << fb.n, core.grh_bound(N))
    assert g2 is not None and 1 < g2 < N


def test_cocycle_identity_canonical_root():
    N, fs = diag.make_instance("mixed", 28, 4, min_factor=core.fb_bound(1 << 28))
    fb = core.FactorBase.make(N)
    rel = core.collect(N, fb, K=600, X=fb.B, T=3000)
    basis = core.F2Basis(rel.kernel, fb.n)
    ys = [core.canonical_root(fb, basis, e) for (_, _, e) in rel.squares]
    seen = 0
    for i in range(len(ys)):
        for j in range(i + 1, min(len(ys), i + 40)):
            if ys[i] is None or ys[j] is None:
                continue
            R = rel.squares[i][1] * rel.squares[j][1] % N
            e = fb.vector(R)
            if e is None:
                continue
            y = core.canonical_root(fb, basis, e)
            if y is None:
                continue
            seen += 1
            assert (y * y - R) % N == 0
            pr = ys[i] * ys[j] % N
            d = y * pow(pr, -1, N) % N
            assert d * d % N == 1
            if d not in (1, N - 1):
                g = gcd(y - pr, N)
                assert 1 < g < N
    assert seen > 0
