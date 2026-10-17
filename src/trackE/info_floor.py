"""HARNESS: exact information-floor census for per-modulus quantities.

Model: for a public modulus m with gcd(m, N) = 1, the hidden pair is
(u, v) = (p mod m, q mod m) in (Z/m)^* x (Z/m)^* with u v ≡ c := N mod m.
Lemma F1 (notes/trackE_analysis.md §3): a function g(u, v) is determined by N mod m
iff it is constant on the fiber {(u, c/u)}; conditional on c and on u being uniform
(true if p mod m, q mod m are independent uniform units), the new information
carried by g is the entropy H(g(U, c/U)).

This module enumerates fibers exactly (small m) and writes one JSONL record per
(m, quantity) with min/mean/max over c of H(g | N mod m), whether g is symmetric,
and whether g is a function of s = u + v given c.  It also checks (empirically) the
uniformity of p mod l given N mod l over random balanced semiprimes.

Run: .venv/bin/python -m src.trackE.info_floor
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict
from math import gcd
from typing import Callable, Dict, List

from src.common.metrics import append_jsonl


def _units(m: int) -> List[int]:
    return [u for u in range(1, m) if gcd(u, m) == 1]


def _ord(u: int, m: int) -> int:
    k, x = 1, u % m
    while x != 1:
        x = x * u % m; k += 1
    return k


def _leg(a: int, l: int) -> int:
    a %= l
    return 0 if a == 0 else (1 if pow(a, (l - 1) // 2, l) == 1 else -1)


def _jac(a: int, n: int) -> int:
    """Jacobi symbol (a/n), n odd positive."""
    a %= n; t = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a %= n
    return t if n == 1 else 0


def quantities(m: int) -> Dict[str, Callable[[int, int], object]]:
    is_prime = m > 2 and all(m % i for i in range(2, int(m ** 0.5) + 1))
    Q: Dict[str, Callable[[int, int], object]] = {
        "N mod m": lambda u, v: u * v % m,
        "p mod m": lambda u, v: u,
        "s=p+q mod m": lambda u, v: (u + v) % m,
        "(p-q)^2 mod m": lambda u, v: (u - v) ** 2 % m,
        "phi(N) mod m": lambda u, v: (u - 1) * (v - 1) % m,
        "p^2+q^2 mod m": lambda u, v: (u * u + v * v) % m,
        "{p,q} mod m (unordered)": lambda u, v: tuple(sorted((u, v))),
        "[m | p-1]": lambda u, v: u == 1,
        "[m | phi(N)]": lambda u, v: u == 1 or v == 1,
        "[p = q mod m]": lambda u, v: u == v,
        "ord_m(p)": lambda u, v: _ord(u, m),
        "lcm(ord_m p, ord_m q)": lambda u, v: math.lcm(_ord(u, m), _ord(v, m)),
        "min(p mod m, q mod m)": lambda u, v: min(u, v),
    }
    if is_prime:
        Q["(N/m) Legendre"] = lambda u, v: _leg(u * v, m)
        Q["(p/m) Legendre"] = lambda u, v: _leg(u, m)
        Q["(p/m)+(q/m)"] = lambda u, v: _leg(u, m) + _leg(v, m)
        for e in (3, 4):
            if (m - 1) % e == 0:
                Q[f"p^((m-1)/{e}) (e={e} power class)"] = (lambda e: lambda u, v: pow(u, (m - 1) // e, m))(e)
    if m % 4 == 0:   # Jacobi symbols (a/N) for a | m/4 are functions of N mod 4a
        a = m // 4
        if a % 2 == 1 and a > 1:
            Q[f"Jacobi ({a}/N)"] = lambda u, v: _jac(a, u * v % m)
            Q[f"Jacobi ({a}/p)"] = lambda u, v: _jac(a, u)
    return Q


def _H(counter: Counter) -> float:
    tot = sum(counter.values())
    return -sum(c / tot * math.log2(c / tot) for c in counter.values())


def census(m: int) -> List[Dict]:
    U = _units(m)
    out = []
    for name, g in quantities(m).items():
        Hs, Hs_given_s, sym = [], [], True
        for c in U:
            cinv = pow(c, -1, m)
            fib = [(u, c * pow(u, -1, m) % m) for u in U]
            cnt = Counter(g(u, v) for u, v in fib)
            Hs.append(_H(cnt))
            bys = defaultdict(Counter)
            for u, v in fib:
                bys[(u + v) % m][g(u, v)] += 1
                if g(u, v) != g(v, u):
                    sym = False
            tot = len(fib)
            Hs_given_s.append(sum(sum(cn.values()) / tot * _H(cn) for cn in bys.values()))
        out.append({
            "m": m, "quantity": name,
            "H_min": round(min(Hs), 4), "H_mean": round(sum(Hs) / len(Hs), 4), "H_max": round(max(Hs), 4),
            "determined_by_N_mod_m": max(Hs) < 1e-12,
            "symmetric": sym,
            "function_of_s_given_c": max(Hs_given_s) < 1e-12,
            "log2_units": round(math.log2(len(U)), 4),
        })
    return out


def uniformity_check(bits: int = 64, count: int = 4000, primes=(3, 5, 7, 11, 13)) -> List[Dict]:
    """Empirical: distribution of p mod l given N mod l over random balanced semiprimes."""
    from src.common.instances import balanced_semiprime
    res = []
    tab = {l: defaultdict(Counter) for l in primes}
    for seed in range(count):
        inst = balanced_semiprime(bits, 10_000 + seed)
        p = inst.labels["p"]
        for l in primes:
            if inst.N % l:
                tab[l][inst.N % l][p % l] += 1
    for l in primes:
        for c, cnt in sorted(tab[l].items()):
            tot = sum(cnt.values()); k = l - 1
            chi2 = sum((cnt.get(u, 0) - tot / k) ** 2 / (tot / k) for u in range(1, l))
            res.append({"l": l, "c": c, "n": tot, "chi2": round(chi2, 3), "dof": k - 1,
                        "H_emp": round(_H(cnt), 4), "H_unif": round(math.log2(k), 4)})
    return res


def main():
    out = "experiments/raw/trackE_info_floor.jsonl"
    import os
    if os.path.exists(out):
        os.remove(out)
    for m in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 8, 16, 32, 64, 12, 20, 28, 9, 25, 27):
        for r in census(m):
            append_jsonl(out, r)
    for r in uniformity_check():
        append_jsonl(out, dict(r, kind="uniformity"))
    print("wrote", out)


if __name__ == "__main__":
    main()
