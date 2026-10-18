"""Part 1: the exact hidden quantity behind s mod 128.  DIAGNOSTIC (knows p and q).

Setting: p != q odd primes with p = q (mod 4), eps = (-1/p), p = eps + 4u, q = eps + 4v, N = pq, s = p + q.
Exact identity (REPORT.md 1.1):  s = eps(N+1) - 16 eps uv,  so s mod 128 <=> uv mod 8 given (N, eps).

2-adic characters.  (Z/32)^* = {+-1} x <5>, 5 of order 8.  Write
    p = eps_p' * 5^a (mod 32),  q = . * 5^b (mod 32),  N = 5^c (mod 32) when N = 1 (mod 4),
with a, b, c in Z/8, c = a + b.  (The sign is (-1/p) = eps because (-1/5) = 1.)
  omega  = order-8 even character mod 32, omega(-1) = 1, omega(5) = zeta_8:   omega(p) = zeta_8^a
  psi    = omega^2 = order-4 even character mod 16, psi(5) = i:              psi(p)   = i^(a mod 4)
  psi^2  = (2/.), the quadratic character mod 8.
Only the products omega(N) = omega(p) omega(q), psi(N), (2/N) are public.

Theorem T1 (proved in notes/tower_results.md section 1; checked exhaustively here):
  (i)   u = eps * a(2a - 1)  (mod 8)                                     [p = q (mod 4)]
  (ii)  uv = (5 - 2c) * a(c - a)  (mod 8)                                 [p = q (mod 4)]
  (iii) if N = 1 (mod 8) (c even), a(c - a) mod 8 depends on a only mod 4, so uv mod 8 -- hence s mod 128 --
        is a function of (N mod 128, eps, psi(p)); omega(p) itself is NOT needed.
  (iv)  N = 1 (mod 8): given (N mod 128, eps, QR_2) the bit h = bit 2 of uv mod 8 (equivalently bit 6 of s beyond
        the known baseline) is undetermined iff psi(p) = psi(q) (iff p = q (mod 16)), and then
            N = 1 (mod 16), (2/p) = +1:  h = [psi(p) = -1]                  = [p = +-7 (mod 16)]
            N = 9 (mod 16), (2/p) = -1:  h = [psi(p) = i] XOR [N = 25 (mod 32)]   ([psi(p) = i] = [p = 5, 11 (mod 16)])
        In the other cells h is a public function:  N = 1 (16), (2/p) = -1: h = [N = 17 (mod 32)];
                                                     N = 9 (16), (2/p) = +1: h = 0.
  (v)   N = 5 (mod 8) (c odd, out of the task's scope): a(c - a) mod 8 changes by 4 when a -> a + 4, so there the
        order-8 value omega(p) (p mod 32) IS needed.
Lemma T2 (proved): for N = 1 (mod 8) and k >= 6, s mod 2^k is a function of (N mod 2^k, p mod 2^(k-3)).

Run:  .venv/bin/python -m src.tower.formula --bound 2000 --out experiments/raw/tower128_formula.json
"""
from __future__ import annotations

import argparse
import json
import time
from collections import defaultdict
from typing import Dict, List, Tuple

from src.repair.verify_r4 import legendre, odd_primes_below

# discrete log base 5 modulo 32 on the subgroup {1 mod 4}
_POW5 = [pow(5, j, 32) for j in range(8)]
_DLOG5 = {x: j for j, x in enumerate(_POW5)}


def index32(x: int) -> Tuple[int, int]:
    """Return (sign, a) with x = sign * 5^a (mod 32), sign = (-1/x) in {+1, -1}, a in Z/8.  x odd."""
    r = x % 32
    if r % 4 == 1:
        return 1, _DLOG5[r]
    return -1, _DLOG5[(-r) % 32]


def omega_exp(x: int) -> int:
    """omega(x) = zeta_8^omega_exp(x): the even order-8 Dirichlet character mod 32 with omega(5) = zeta_8."""
    return index32(x)[1]


def psi_exp(x: int) -> int:
    """psi(x) = i^psi_exp(x): the even order-4 character mod 16 with psi(5) = i (psi = omega^2)."""
    return index32(x)[1] % 4


def hidden_bit(p: int, q: int) -> int:
    """h = bit 2 of (uv mod 8).  With uv mod 4 public (REPORT 1.1), s mod 128 <=> h given (N mod 128, eps, QR_2)."""
    eps = 1 if p % 4 == 1 else -1
    u, v = (p - eps) // 4, (q - eps) // 4
    return ((u * v) % 8) >> 2


def uv_mod8_formula(p: int, q: int) -> int:
    """T1(ii): uv = (5 - 2c) a (c - a) mod 8, using only the 2-adic indices a (of p) and c (of N)."""
    a = omega_exp(p)
    c = omega_exp(p * q)
    return ((5 - 2 * c) * a * (c - a)) % 8


def h_from_psi(N: int, eps: int, psi_p: int) -> int:
    """T1(iv) as a function of public (N mod 32, eps) and the hidden psi(p) exponent (N = 1 mod 8).  eps unused:
    the rule does not depend on eps (kept in the signature to make the key explicit)."""
    two_sq = psi_p % 2 == 0  # psi(p)^2 = (2/p)
    if N % 16 == 1:
        return int(psi_p == 2) if two_sq else int(N % 32 == 17)
    assert N % 16 == 9
    return int(psi_p == 1) ^ int(N % 32 == 25) if not two_sq else 0


def run(bound: int = 2000) -> Dict[str, object]:
    t0 = time.time()
    P = odd_primes_below(bound)
    fails: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    counts: Dict[str, int] = defaultdict(int)
    groups_T2: Dict[int, Dict[tuple, set]] = {k: defaultdict(set) for k in range(6, 12)}
    groups_pub: Dict[tuple, set] = defaultdict(set)
    groups_pub_psi: Dict[tuple, set] = defaultdict(set)
    groups_pub_psi_amb: Dict[tuple, set] = defaultdict(set)
    amb_cells = defaultdict(int)
    groups_n5_psi: Dict[tuple, set] = defaultdict(set)
    groups_n5_omega: Dict[tuple, set] = defaultdict(set)
    for i, p in enumerate(P):
        for q in P[i + 1:]:
            if p % 4 != q % 4:
                continue
            N, s = p * q, p + q
            eps = 1 if p % 4 == 1 else -1
            u, v = (p - eps) // 4, (q - eps) // 4
            counts["p_eq_q_mod4"] += 1
            # T1(i), T1(ii): all pairs p = q (mod 4)
            a, b, c = omega_exp(p), omega_exp(q), omega_exp(N)
            if index32(p)[0] != eps or index32(q)[0] != eps or index32(N)[0] != 1 or (a + b - c) % 8:
                fails["index_bookkeeping"].append((p, q))
            if u % 8 != (eps * a * (2 * a - 1)) % 8 or v % 8 != (eps * b * (2 * b - 1)) % 8:
                fails["T1_i"].append((p, q))
            if (u * v) % 8 != uv_mod8_formula(p, q):
                fails["T1_ii"].append((p, q))
            if s % 128 != (eps * (N + 1) - 16 * eps * uv_mod8_formula(p, q)) % 128:
                fails["s128_from_formula"].append((p, q))
            if N % 8 != 1:
                counts["N_5_mod_8"] += 1
                groups_n5_psi[(N % 128, eps, psi_exp(p))].add(s % 128)
                groups_n5_omega[(N % 128, eps, omega_exp(p))].add(s % 128)
                continue
            counts["N_1_mod_8"] += 1
            two_sq = legendre(2, p) == 1
            if legendre(2, q) != legendre(2, p):
                fails["2p_eq_2q"].append((p, q))
            h = hidden_bit(p, q)
            # T1(iv): rule for h
            if h != h_from_psi(N, eps, psi_exp(p)) or h != h_from_psi(N, eps, psi_exp(q)):
                fails["T1_iv_rule"].append((p, q))
            # public key and the determination tests
            key = (N % 128, eps, int(two_sq))
            groups_pub[key].add(s % 128)
            groups_pub_psi[key + (psi_exp(p),)].add(s % 128)  # T1(iii): (N mod 128, eps, psi(p)) determines s mod 128
            same_psi = psi_exp(p) == psi_exp(q)
            if same_psi != (p % 16 == q % 16):
                fails["psi_equal_iff_p_eq_q_mod16"].append((p, q))
            amb_cells[(N % 16, int(two_sq), same_psi)] += 1
            # T2: s mod 2^k function of (N mod 2^k, p mod 2^(k-3))
            for k in groups_T2:
                groups_T2[k][(N % (1 << k), p % (1 << (k - 3)))].add(s % (1 << k))
                groups_T2[k][(N % (1 << k), q % (1 << (k - 3)))].add(s % (1 << k))
    def amb(groups):
        return sum(1 for vals in groups.values() if len(vals) > 1)
    # predicted ambiguity: undetermined iff psi(p) = psi(q).  Check against the actual table of keys.
    pub_amb_keys = {k for k, vals in groups_pub.items() if len(vals) > 1}
    return {
        "bound": bound,
        "counts": dict(counts),
        "failures": {k: {"n": len(v), "first": v[:5]} for k, v in fails.items()},
        "checks": ["index_bookkeeping", "T1_i", "T1_ii", "s128_from_formula", "2p_eq_2q", "T1_iv_rule",
                   "psi_equal_iff_p_eq_q_mod16"],
        "T1_iii_keys": len(groups_pub_psi), "T1_iii_ambiguous_keys": amb(groups_pub_psi),
        "pub_keys": len(groups_pub), "pub_ambiguous_keys": len(pub_amb_keys),
        "T1_v_N5mod8_by_psi": {"keys": len(groups_n5_psi), "ambiguous_keys": amb(groups_n5_psi)},
        "T1_v_N5mod8_by_omega": {"keys": len(groups_n5_omega), "ambiguous_keys": amb(groups_n5_omega)},
        "T2": {k: {"keys": len(g), "ambiguous_keys": amb(g)} for k, g in groups_T2.items()},
        "cells_N16_twoSq_samePsi": {f"N%16={k[0]},QR2={k[1]},psi_p==psi_q={k[2]}": n
                                    for k, n in sorted(amb_cells.items())},
        "seconds": round(time.time() - t0, 2),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bound", type=int, default=2000)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    res = run(args.bound)
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(res, fh, sort_keys=True, indent=1)
            fh.write("\n")
    print(json.dumps(res, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
