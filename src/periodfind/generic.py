"""HARNESS: generic-ring observables (decision trees with zero tests) and Theorem GR-F checks.

A tree node tests whether T(y) = y^m - c is 0 mod N and branches; leaves carry complex outputs.
`run` evaluates the tree mod N, mod p and mod q, and records whether any test on the N-run splits N
(zero at exactly one of p, q). Theorem GR-F: mass of F = Psi(a^t) on frequencies j with
(r/g) not dividing j is <= 8 * eps_split, with g = gcd(r_p, r_q).
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class Node:
    m: int = 0
    c: int = 0
    yes: Optional["Node"] = None
    no: Optional["Node"] = None
    leaf: Optional[complex] = None


def random_tree(depth: int, rng: random.Random, N: int, exps: list) -> Node:
    if depth == 0 or rng.random() < 0.2:
        return Node(leaf=complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) / math.sqrt(2))
    m = rng.choice(exps)
    c = 1 if rng.random() < 0.7 else rng.randrange(1, N)
    return Node(m=m, c=c, yes=random_tree(depth - 1, rng, N, exps),
                no=random_tree(depth - 1, rng, N, exps))


def run(tree: Node, y: int, mod: int):
    """Return (output, split_seen) of the run mod `mod` (split_seen only meaningful for mod N)."""
    node = tree
    while node.leaf is None:
        v = (pow(y, node.m, mod) - node.c) % mod
        node = node.yes if v == 0 else node.no
    return node.leaf


def run_N(tree: Node, y: int, N: int, p: int, q: int):
    node, split = tree, False
    while node.leaf is None:
        v = (pow(y, node.m, N) - node.c) % N
        zp, zq = v % p == 0, v % q == 0
        if zp != zq:
            split = True
        node = node.yes if v == 0 else node.no
    return node.leaf, split


def check_tree(tree: Node, orb: np.ndarray, N: int, p: int, q: int, r: int, g: int):
    vals = np.empty(r, dtype=np.complex128)
    nsplit = 0
    for t in range(r):
        out, sp = run_N(tree, int(orb[t]), N, p, q)
        vals[t] = out
        nsplit += sp
    eps = nsplit / r
    Fh = np.fft.fft(vals) / r
    P = np.abs(Fh) ** 2
    j = np.arange(r)
    visible = (j % (r // g)) != 0
    vis_mass = float(P[visible].sum())
    return {"eps_split": eps, "visible_mass": vis_mass, "bound": 8 * eps,
            "ok": vis_mass <= 8 * eps + 1e-9, "nonconstant": float(np.var(vals)) > 1e-12}
