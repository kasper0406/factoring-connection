"""HARNESS (factor-aware): curves, exact local data via PARI, instance sampling, exact TW success probability rho,
exceptional-set membership, Chebotarev predictions.  Never imported by src/ecount/zn_curve.py."""
from __future__ import annotations

import random
from functools import lru_cache
from math import gcd
from typing import Dict, List, Tuple

import cypari2
import gmpy2

from src.common.instances import balanced_semiprime, random_prime

pari = cypari2.Pari()
pari.allocatemem(64 * 10 ** 6, 10 ** 9, silent=True)
pari.default("nbthreads", 1)  # PARI multithreaded polmodular deadlocks under cypari2 (observed in E1, D1)

# short Weierstrass models y^2 = x^3 + a x + b
CURVES: Dict[str, Tuple[int, int]] = {
    "C1_11a1": (-13392, -1080432),   # 11a1 = y^2 + y = x^3 - x^2 - 10x - 20, (c4, c6) = (496, 20008)
    "C2_37a1": (-1296, 11664),       # 37a1 = y^2 + y = x^3 - x, (c4, c6) = (48, -216)
    "C3_2_3": (2, 3),                # y^2 = x^3 + 2x + 3 (f has the root -1: rational 2-torsion)
    "CM1_x3mx": (-1, 0),             # CM by Z[i]; supersingular at p = 3 mod 4
    "CM2_x3p1": (0, 1),              # CM by Z[zeta_3]; supersingular at p = 2 mod 3
}
NONCM = ("C1_11a1", "C2_37a1", "C3_2_3")
CM = ("CM1_x3mx", "CM2_x3p1")


def disc_part(a: int, b: int) -> int:
    return 6 * (4 * a ** 3 + 27 * b ** 2)


def good(N: int, a: int, b: int) -> bool:
    return gcd(N, disc_part(a, b)) == 1


@lru_cache(maxsize=None)
def _ell(a: int, b: int):
    return pari.ellinit([0, 0, 0, a, b])


def card(a: int, b: int, p: int) -> int:
    return int(pari.ellcard(_ell(a, b), p))


def _qnr(p: int) -> int:
    d = 2
    while gmpy2.legendre(d, p) != -1:
        d += 1
    return d


def twist_coeffs(a: int, b: int, p: int) -> Tuple[int, int]:
    d = _qnr(p)
    return (a * d * d) % p, (b * d ** 3) % p


def group(a: int, b: int, p: int) -> List[int]:
    """Cyclic decomposition [m2, m1] (m1 | m2) of E(F_p), via PARI ellgroup over F_p."""
    E = pari.ellinit([0, 0, 0, a % p, b % p], p)
    return [int(x) for x in pari.ellgroup(E)]


def roots_f(a: int, b: int, p: int) -> int:
    return len(pari.polrootsmod(pari(f"x^3+({a})*x+({b})"), p))


def local_data(a: int, b: int, p: int) -> Dict:
    A = card(a, b, p)
    at, bt = twist_coeffs(a, b, p)
    Gt = group(at, bt, p)
    At = 1
    for m in Gt:
        At *= m
    assert A + At == 2 * p + 2
    return {"p": p, "A": A, "At": At, "twist_group": Gt, "r": roots_f(a, b, p), "ap": p + 1 - A}


def H_size(group_cyc: List[int], M: int) -> int:
    s = 1
    for m in group_cyc:
        s *= gcd(m, M)
    return s


def exp_divides(group_cyc: List[int], M: int) -> bool:
    return (M % group_cyc[0]) == 0 if group_cyc else True


def rho_exact(dp: Dict, dq: Dict) -> Tuple[float, float, bool]:
    """Exact (as Fraction-free ints -> float for display) TW split probability per trial (lower bound rho, plus
    the delta-gcd bonus) and exceptional-set membership."""
    M = dp["A"] * dq["A"]
    p, q = dp["p"], dq["p"]
    Hp = H_size(dp["twist_group"], M)
    Hq = H_size(dq["twist_group"], M)
    num = (dp["A"] - 1 - dp["r"]) * (dq["At"] - Hq) + (dp["At"] - Hp) * (dq["A"] - 1 - dq["r"])
    rho = num / (4 * p * q)
    bonus = (dp["r"] * (q - dq["r"]) + dq["r"] * (p - dp["r"])) / (p * q)
    exc = exp_divides(dp["twist_group"], M) and exp_divides(dq["twist_group"], M)
    return rho, bonus, exc


def rho_num_den(dp: Dict, dq: Dict) -> Tuple[int, int, int, int]:
    M = dp["A"] * dq["A"]
    p, q = dp["p"], dq["p"]
    Hp = H_size(dp["twist_group"], M)
    Hq = H_size(dq["twist_group"], M)
    num = (dp["A"] - 1 - dp["r"]) * (dq["At"] - Hq) + (dp["At"] - Hp) * (dq["A"] - 1 - dq["r"])
    bon = 4 * (dp["r"] * (q - dq["r"]) + dq["r"] * (p - dp["r"]))
    return num, bon, 4 * p * q, M


def instance(curve: str, n: int, seed: int) -> Tuple[int, int, int]:
    a, b = CURVES[curve]
    k = 0
    while True:
        inst = balanced_semiprime(n, seed + 1_000_000 * k)
        p, q = inst.factors
        if good(inst.N, a, b):
            return inst.N, p, q
        k += 1


def instance_cm1_33(n: int, seed: int) -> Tuple[int, int, int]:
    rng = random.Random(f"ecount-cm1-33|{n}|{seed}")
    while True:
        p = random_prime(rng, (n + 1) // 2, top2=True)
        q = random_prime(rng, n // 2, top2=True)
        if p != q and p % 4 == 3 and q % 4 == 3 and (p * q).bit_length() == n:
            p, q = sorted((p, q))
            return p * q, p, q


def affine_count(a: int, b: int, p: int, q: int) -> int:
    return (card(a, b, p) - 1) * (card(a, b, q) - 1)


def count(a: int, b: int, p: int, q: int) -> int:
    return card(a, b, p) * card(a, b, q)


def twist_count(a: int, b: int, d: int, p: int, q: int, affine: bool) -> int:
    ad, bd = a * d * d, b * d ** 3
    if affine:
        return (card(ad, bd, p) - 1) * (card(ad, bd, q) - 1)
    return card(ad, bd, p) * card(ad, bd, q)


def factorint(M: int) -> Dict[int, int]:
    f = pari.factor(M)
    return {int(f[0][i]): int(f[1][i]) for i in range(len(f[0]))}


def primes_window(bbits: int) -> List[int]:
    lo, hi = 3 << (bbits - 2), 1 << bbits
    return [int(x) for x in pari.primes([lo, hi - 1])]


# ---------------- Chebotarev predictions (GL2(F_l) uniform, independent Frobenius) ----------------

def gl2_table(l: int) -> Dict[Tuple[int, int], int]:
    """Counts of (det g, det(1 - g)) over GL2(F_l)."""
    T: Dict[Tuple[int, int], int] = {}
    for a in range(l):
        for b in range(l):
            for c in range(l):
                for d in range(l):
                    det = (a * d - b * c) % l
                    if det == 0:
                        continue
                    m = ((1 - a) * (1 - d) - b * c) % l
                    T[(det, m)] = T.get((det, m), 0) + 1
    return T


def chebotarev_adv_l(l: int) -> float:
    """Predicted Bayes advantage for beta = (p/l) on {(N/l) = 1}, from (N mod l, M mod l); ordered iid pairs."""
    T = gl2_table(l)
    tot = sum(T.values())
    joint: Dict[Tuple[int, int], List[float]] = {}
    mass = 0.0
    for (d1, m1), c1 in T.items():
        for (d2, m2), c2 in T.items():
            s1 = gmpy2.legendre(d1, l)
            s2 = gmpy2.legendre(d2, l)
            if s1 != s2:
                continue
            w = c1 * c2 / tot ** 2
            key = ((d1 * d2) % l, (m1 * m2) % l)
            joint.setdefault(key, [0.0, 0.0])[0 if s1 == 1 else 1] += w
            mass += w
    return sum(max(v) for v in joint.values()) / mass - 0.5


def chebotarev_adv_2() -> float:
    """S3 Frobenius for irreducible f with non-square disc: beta = sign; predictor from M mod 2."""
    cls = [(1 / 6, +1, 3), (1 / 2, -1, 1), (1 / 3, +1, 0)]
    joint: Dict[int, List[float]] = {}
    mass = 0.0
    for w1, s1, r1 in cls:
        for w2, s2, r2 in cls:
            if s1 != s2:
                continue
            Modd = int(r1 == 0 and r2 == 0)
            joint.setdefault(Modd, [0.0, 0.0])[0 if s1 == 1 else 1] += w1 * w2
            mass += w1 * w2
    return sum(max(v) for v in joint.values()) / mass - 0.5
