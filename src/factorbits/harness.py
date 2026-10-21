"""HARNESS oracles for factorbits (they factor or brute-force discrete logs; never used by equiv.py)."""
from __future__ import annotations

from typing import Dict

import sympy

from .equiv import Vec, value


def fact(x: int) -> Vec:
    return {int(p): int(e) for p, e in sympy.factorint(x).items()} if x > 1 else {}


def add_oracle(a: Vec, b: Vec) -> Vec:
    return fact(value(a) + value(b))


def succ_oracle(a: Vec) -> Vec:
    return fact(value(a) + 1)


def support_oracle(a: Vec, b: Vec) -> set:
    return set(fact(value(a) + value(b)))


def make_expmod_oracle(m: int):
    def expmod(a: Vec, b: Vec) -> Vec:
        return {p: e % m for p, e in fact(value(a) + value(b)).items() if e % m}
    return expmod


def split_oracle(x: int) -> int:
    f = sorted(fact(x))
    return f[0] if f[0] != x else None


def make_dlog_table(p: int, g: int) -> Dict[int, int]:
    tab, y = {}, 1
    for k in range(p - 1):
        tab[y] = k
        y = y * g % p
    assert len(tab) == p - 1, "g must be a generator"
    return tab


def make_zech_oracle(p: int, g: int, tab: Dict[int, int]):
    def zech(k: int) -> int:
        return tab[(1 + pow(g, k, p)) % p]
    return zech
