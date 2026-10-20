"""Tests for agent escapeR2 (src/escapeR2)."""
from fractions import Fraction as Fr

from src.escapeR2 import resplit, sizemodel as SM, transfer as TR


# ---------------------------------------------------------------- transfer (Lemma TI)
def test_group_closure_orders():
    G, blocks = TR.direct_product([3, 4])
    assert len(G) == 6 * 24
    assert len(TR.symmetric(4)) == 24
    assert {k: len(v) for k, v in TR.census_groups().items()} == {"S4": 24, "D4": 8, "A4": 12, "C2^3": 8,
                                                                    "S3xC2": 12}


def test_double_coset_basis_partitions():
    G = TR.symmetric(4)
    H = TR.stabiliser(G, 0)
    for D in TR.cyclic_subgroups(G):
        basis = TR.image_basis(G, H, D)
        cos, _ = TR.left_cosets(G, D)
        union = 0
        for v in basis:
            assert union & v == 0
            union |= v
        assert union == (1 << len(cos)) - 1


def test_v4_transfer_exists():
    G, Hs = TR.v4()
    D = TR.closure([TR.identity(4)], 4)
    rep = TR.transfer_report(G, Hs, D)
    assert rep["tau"] == 1 and rep["trivial_meet"] == 0
    # pairwise there is no transfer: Q(sqrt a) ∩ Q(sqrt b) = Q
    assert TR.transfer_report(G, Hs[:2], D)["tau2"] == 0
    assert TR.transfer_report(G, Hs[:2], D)["tau"] == 0


def test_disjoint_fields_no_transfer():
    G, blocks = TR.direct_product([3, 3, 3])
    Hs = [TR.stabiliser(G, b[0]) for b in blocks]
    for D in TR.cyclic_subgroups(G):
        rep = TR.transfer_report(G, Hs, D)
        assert rep["trivial_meet"] == 1 and rep["tau"] == 0


def test_tau_bounded_by_meet_field():
    G = TR.census_groups()["D4"]
    subs = TR.subgroups_two_generated(G)
    for D in TR.cyclic_subgroups(G):
        for i in range(len(subs)):
            for j in range(len(subs)):
                for k in range(len(subs)):
                    if len({i, j, k}) == 3:
                        rep = TR.transfer_report(G, [subs[i], subs[j], subs[k]], D)
                        assert rep["tau"] <= rep["bound"]


# ---------------------------------------------------------------- size model
def test_lp_basic():
    assert SM.lp_feasible([[Fr(1)]], [Fr(1)])
    assert SM.lp_feasible([[Fr(-1)]], [Fr(-2)])           # -x <= -2, i.e. x >= 2
    assert not SM.lp_feasible([[Fr(1)], [Fr(-1)]], [Fr(1), Fr(-2)])  # x <= 1 and x >= 2


def test_one_level_floor_is_two_fifths():
    assert SM.chain_feasible(Fr(2, 5), 0, None)
    assert not SM.chain_feasible(Fr(2, 5) - Fr(1, 10 ** 6), 0, None)


def test_level_thresholds():
    for r in range(4):
        a = SM.closed_no_budget(r)
        assert SM.chain_feasible(a, r, None)
        assert not SM.chain_feasible(a - Fr(1, 10 ** 6), r, None)


def test_budget_one_gives_no_gain():
    for r in range(5):
        assert not SM.chain_feasible(Fr(2, 5) - Fr(1, 10 ** 6), r, Fr(1))


def test_budget_three_halves_reaches_three_eighths():
    assert SM.chain_feasible(Fr(3, 8), 1, Fr(3, 2))
    assert not SM.chain_feasible(Fr(3, 8), 1, Fr(3, 2) - Fr(1, 100))


# ---------------------------------------------------------------- re-split obstruction
def test_resplit_coset_obstruction_small():
    d = resplit.build(1)
    S = resplit.products(d["good"], 2, d["Yp"] // 2, d["Yp"])
    assert len(S) > 0
    nonsq = (-1,) + (1,) * (len(d["ps"]) - 1)
    assert resplit.pair_count(S, d["ps"], nonsq) == 0
    assert resplit.pair_count(S, d["ps"], (1,) * len(d["ps"])) > 0
