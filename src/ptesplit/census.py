"""Exhaustive PTE fibre census (DIAGNOSTIC; no N, no factoring).

Multisets of size M from {0..R} (= [-B, B] after translation, R = 2B).  For s prescribed power sums
(k = M - s free coefficients) we group multisets by (p_1..p_s) using two independent 64-bit
polynomial hashes of the power-sum vector (false-merge probability ~ T^2 / 2^128, disclosed).
"""
from __future__ import annotations

import math
from itertools import combinations_with_replacement

import numpy as np

_C1 = np.array([0x9E3779B97F4A7C15, 0xC2B2AE3D27D4EB4F, 0x165667B19E3779F9, 0xD6E8FEB86659FD93,
                0xFF51AFD7ED558CCD, 0xC4CEB9FE1A85EC53, 0x2545F4914F6CDD1D, 0x94D049BB133111EB,
                0xBF58476D1CE4E5B9, 0x7FB5D329728EA185], dtype=np.uint64)
_C2 = np.array([0xA0761D6478BD642F, 0xE7037ED1A0B428DB, 0x8EBC6AF09C88C6E3, 0x589965CC75374CC3,
                0x1D8E4E27C47D124F, 0x2127599BF4325C37, 0x880355F21E6D1965, 0xDB4F0B9175AE2165,
                0x6C8E9CF570932BD5, 0x9C4E2D1F3B5A7E61], dtype=np.uint64)


def multisets(M: int, R: int) -> np.ndarray:
    """All nondecreasing M-tuples over {0..R}, shape (T, M), dtype uint8. Built column by column."""
    assert R < 256
    # level-wise construction: arrays of tuples ending at value v
    cur = np.arange(R + 1, dtype=np.uint8).reshape(-1, 1)
    for _ in range(M - 1):
        last = cur[:, -1].astype(np.int64)
        reps = (R + 1 - last)  # next value in [last, R]
        idx = np.repeat(np.arange(cur.shape[0]), reps)
        starts = np.repeat(last, reps)
        offs = np.arange(idx.shape[0]) - np.repeat(np.cumsum(reps) - reps, reps)
        nxt = (starts + offs).astype(np.uint8)
        cur = np.concatenate([cur[idx], nxt.reshape(-1, 1)], axis=1)
    return cur


def count_multisets(M: int, R: int) -> int:
    return math.comb(R + M, M)


def power_sum(arr: np.ndarray, i: int) -> np.ndarray:
    out = np.zeros(arr.shape[0], dtype=np.int64)
    for j in range(arr.shape[1]):
        c = arr[:, j].astype(np.int64)
        out += c ** i
    return out


def fibre_hashes(arr: np.ndarray, s_list: list[int]) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """Cumulative double hashes of (p_1..p_s) for each s in s_list."""
    smax = max(s_list)
    h1 = np.zeros(arr.shape[0], dtype=np.uint64)
    h2 = np.zeros(arr.shape[0], dtype=np.uint64)
    out = {}
    with np.errstate(over="ignore"):
        for i in range(1, smax + 1):
            p = power_sum(arr, i).astype(np.uint64)
            h1 = h1 * np.uint64(0x100000001B3) + p * _C1[(i - 1) % 10]
            h2 = (h2 ^ (p * _C2[(i - 1) % 10])) * np.uint64(0x9FB21C651E98DF25)
            if i in s_list:
                out[i] = (h1.copy(), h2.copy())
    return out


def group_sizes(h1: np.ndarray, h2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (order, sizes-per-element-group): sorted order and group boundaries."""
    order = np.lexsort((h2, h1))
    a, b = h1[order], h2[order]
    brk = np.ones(len(a), dtype=bool)
    brk[1:] = (a[1:] != a[:-1]) | (b[1:] != b[:-1])
    starts = np.flatnonzero(brk)
    sizes = np.diff(np.append(starts, len(a)))
    return order, np.stack([starts, sizes], axis=1)


def box_model_pairs(M: int, R: int, s: int) -> float:
    """Independent-uniform-box prediction of unordered colliding pairs: T^2 / (2 prod_i (M R^i + 1))."""
    T = count_multisets(M, R)
    logbox = sum(math.log(M * R ** i + 1) for i in range(1, s + 1))
    return 0.5 * T * T * math.exp(-logbox)


def injection_bound(R: int, k: int) -> int:
    """|fibre| <= C(R + k, k): multisets of size k over {0..R}."""
    return math.comb(R + k, k)


def ideal_constants(arr: np.ndarray, members: np.ndarray) -> list[int]:
    """For a k = 1 fibre (members share p_1..p_{M-1}): the constants e_M differences, i.e.
    P_A - P_B = (-1)^M (prod A - prod B).  Returns |prod A_i - prod A_0| for i >= 1 (exact ints)."""
    prods = [math.prod(int(v) for v in arr[i]) for i in members]
    return [abs(p - prods[0]) for p in prods[1:]]


def census_cell(M: int, R: int, ks=(1, 2, 3), keep_examples: int = 3) -> list[dict]:
    arr = multisets(M, R)
    T = arr.shape[0]
    s_list = [M - k for k in ks if M - k >= 1]
    hs = fibre_hashes(arr, s_list)
    rows = []
    for k in ks:
        s = M - k
        if s < 1:
            continue
        h1, h2 = hs[s]
        order, gs = group_sizes(h1, h2)
        sizes = gs[:, 1]
        nontriv = sizes[sizes >= 2]
        pairs = int((nontriv * (nontriv - 1) // 2).sum())
        maxf = int(sizes.max())
        row = dict(M=M, R=R, k=k, s=s, T=T, n_fibres=int(len(sizes)), n_nontrivial=int(len(nontriv)),
                   pairs=pairs, max_fibre=maxf, box_pairs=box_model_pairs(M, R, s),
                   injection_bound=injection_bound(R, k),
                   size_hist={int(a): int(b) for a, b in zip(*np.unique(sizes[sizes >= 2], return_counts=True))})
        if k == 1:
            # Kleiman check: every constant difference divisible by (M-1)!, plus min nonzero constant
            fact = math.factorial(M - 1)
            viol = 0
            cmin = None
            exs = []
            big = np.flatnonzero(sizes >= 2)
            for g in big:
                st, sz = gs[g]
                mem = order[st:st + sz]
                consts = ideal_constants(arr, mem)
                for c in consts:
                    if c == 0 or c % fact:
                        viol += 1
                    cmin = c if cmin is None or (c and c < cmin) else cmin
                if len(exs) < keep_examples:
                    exs.append([[int(v) for v in arr[i]] for i in mem])
            row.update(kleiman_violations=viol, min_constant=cmin, examples=exs,
                       bp_proxy=2 * math.sqrt(R + 1) + 2)
        rows.append(row)
    return rows


def brute_fibres(M: int, R: int, s: int) -> dict:
    """Pure-Python reference (small cells only)."""
    d: dict = {}
    for a in combinations_with_replacement(range(R + 1), M):
        key = tuple(sum(v ** i for v in a) for i in range(1, s + 1))
        d.setdefault(key, []).append(a)
    return d
