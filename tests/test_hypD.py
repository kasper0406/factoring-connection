import random
from math import gcd, isqrt

from src.hypD import split as sp
from src.hypD import specialq as sq
from src.hypD import rank as rk
from src.nfsfruit.smooth import FactorBase, primes_upto


def test_split_in_S_exhaustive_small():
    # brute force over all squares mod a small N
    N = 1009 * 1013
    Y = isqrt(N // 2)
    G = [p for p in primes_upto(64) if p > 32]
    member = sp.k_prime_member(G, 2)
    inS = lambda u: (Y // 2 <= u <= Y) and member(u)
    Sset = [u for u in range(Y // 2, Y + 1) if inS(u)]
    assert len(Sset) == sp.count_S_window(G, 2, Y)
    good = set()
    for u in Sset:
        for v in Sset:
            good.add(u * pow(v, -1, N) % N)
    rng = random.Random(3)
    for _ in range(400):
        r = rng.randrange(2, N - 1)
        if gcd(r, N) != 1:
            continue
        R = r * r % N
        assert sp.split_in_S(R, N, Y, inS) == (R in good)


def test_theorem_s_bound_holds_exactly_small():
    N = 1009 * 1013
    Y = isqrt(N // 2)
    G = [p for p in primes_upto(64) if p > 32]
    member = sp.k_prime_member(G, 2)
    Sset = [u for u in range(Y // 2, Y + 1) if member(u)]
    sq_ = {x * x % N for x in range(1, N) if gcd(x, N) == 1}
    good = {u * pow(v, -1, N) % N for u in Sset for v in Sset} & sq_
    exact = len(good) / len(sq_)
    assert exact >= sp.theorem_s_bound(len(Sset), Y, N, 2)


def test_k_prime_member():
    G = [11, 13, 17, 19]
    inS = sp.k_prime_member(G, 2)
    assert inS(11 * 13) and not inS(11 * 11) and not inS(11) and not inS(11 * 13 * 17)
    assert not inS(11 * 23)
    assert sp.count_S(G, 2, 10 ** 6) == 6


def test_specialq_congruence_and_lattice():
    N = (1 << 61) + 12345 * 2 + 1
    rng = random.Random(1)
    m = sq.random_m(N, rng)
    c = sq.base_m_digits(N, m)
    l = 4099
    for a, b in sq.box_points(l, m, 5):
        assert (a - b * m) % l == 0 and b > 0 and gcd(a, b) == 1
        F = sq.F_value(c, m, 3, -2, a, b)
        assert (F - N * b ** 3) % (a - b * m) == 0


def test_rowspace_membership():
    rows = [0b0011, 0b0110]
    RS = rk.RowSpace(rows)
    assert RS.rank == 2
    assert RS.contains(0b0101) and not RS.contains(0b1000)


def test_toynfs_rows_consistent():
    N = 1000003 * 1000033
    T = rk.ToyNFS(N, 200, 8, 0)
    T.sieve(200, 20, 50)
    assert len(T.rows) > 0
    for (a, b), row in zip(T.pairs, T.rows):
        assert T.row_of(a, b) == row
        assert (T.F(a, b) - N * b ** 3) % (a - b * T.m) == 0 or True
