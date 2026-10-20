"""Straight-line programs over Z: generation, evaluation mod p (vectorised), polynomial reduction.

An SLP is a tuple (leaves, ops): leaves is a list of 'x' or int constants; ops is a list of
(op, i, j) with op in '+-*' referring to earlier nodes (leaves first). The output is the last node.
Constant-free in the slpcirc sense when the integer leaves are small fixed constants.
"""
from __future__ import annotations
import random
from math import lcm
import numpy as np


def random_slp(size: int, rng: random.Random, consts=(1, 2, 3)):
    leaves = ['x'] + list(consts)
    ops = []
    n = len(leaves)
    for k in range(size):
        op = rng.choice('+-*')
        # bias towards recent nodes so the output depends on most of the program
        i = rng.randrange(n + k)
        j = rng.randrange(n + k)
        if k > 0 and rng.random() < 0.5:
            i = n + k - 1
        ops.append((op, i, j))
    return (leaves, ops)


def composition_slp(cs):
    """(x^2 + c_k) o ... o (x^2 + c_1) for cs = (c_1, ..., c_k)."""
    leaves = ['x'] + sorted(set(cs))
    idx = {c: 1 + sorted(set(cs)).index(c) for c in set(cs)}
    ops = []
    cur = 0
    n = len(leaves)
    for c in cs:
        ops.append(('*', cur, cur)); cur = n + len(ops) - 1
        ops.append(('+', cur, idx[c])); cur = n + len(ops) - 1
    return (leaves, ops)


def power_slp(M: int, plus_x: bool = False):
    """x^M by left-to-right square-and-multiply; optionally x^M + x."""
    leaves = ['x']
    ops = []
    cur = 0
    for bit in bin(M)[3:]:
        ops.append(('*', cur, cur)); cur = len(leaves) + len(ops) - 1
        if bit == '1':
            ops.append(('*', cur, 0)); cur = len(leaves) + len(ops) - 1
    if plus_x:
        ops.append(('+', cur, 0))
    return (leaves, ops)


def formal_degree(slp) -> int:
    leaves, ops = slp
    d = [1 if l == 'x' else 0 for l in leaves]
    for op, i, j in ops:
        d.append(d[i] + d[j] if op == '*' else max(d[i], d[j]))
    return d[-1]


def eval_mod(slp, xs: np.ndarray, p: int) -> np.ndarray:
    """Evaluate on an int64 array xs (values in [0, p)) mod p. Requires p < 2^31."""
    leaves, ops = slp
    vals = [xs % p if l == 'x' else np.full_like(xs, l % p) for l in leaves]
    for op, i, j in ops:
        a, b = vals[i], vals[j]
        if op == '+':
            vals.append((a + b) % p)
        elif op == '-':
            vals.append((a - b) % p)
        else:
            vals.append((a * b) % p)
    return vals[-1]


def eval_int(slp, x: int, m: int) -> int:
    leaves, ops = slp
    vals = [x % m if l == 'x' else l % m for l in leaves]
    for op, i, j in ops:
        a, b = vals[i], vals[j]
        vals.append((a + b) % m if op == '+' else (a - b) % m if op == '-' else (a * b) % m)
    return vals[-1]


def _cyc_mul(a, b, n, p):
    """Multiply polynomials in F_p[y]/(y^n - 1), coefficient arrays of length <= n."""
    c = np.convolve(a, b) % p
    if len(c) > n:
        out = np.zeros(n, dtype=np.int64)
        for s in range(0, len(c), n):
            seg = c[s:s + n]
            out[:len(seg)] += seg
        return out % p
    return c


def _cyc_add(a, b, p, sign=1):
    n = max(len(a), len(b))
    out = np.zeros(n, dtype=np.int64)
    out[:len(a)] += a
    out[:len(b)] += sign * b
    return out % p


def reduced_poly(slp, p: int):
    """Coefficients of the SLP output as a function on F_p^*, i.e. mod (y^{p-1} - 1, p).
    O(len^2) per multiplication; intended for p < 2^14 or small formal degree."""
    n = p - 1
    leaves, ops = slp
    vals = []
    for l in leaves:
        if l == 'x':
            vals.append(np.array([0, 1], dtype=np.int64))
        else:
            vals.append(np.array([l % p], dtype=np.int64))
    for op, i, j in ops:
        a, b = vals[i], vals[j]
        if op == '*':
            vals.append(_cyc_mul(a, b, n, p))
        else:
            vals.append(_cyc_add(a, b, p, 1 if op == '+' else -1))
    out = vals[-1]
    if len(out) > n:
        full = np.zeros(n, dtype=np.int64)
        for s in range(0, len(out), n):
            seg = out[s:s + n]
            full[:len(seg)] += seg
        out = full % p
    return out


def weil_width(coeffs, p: int) -> int:
    """Minimal a + b over lifts of the nonzero exponents of a function on F_p^* (exponents mod p-1,
    constant term ignored) to a Laurent polynomial with exponents in [-a, b]. This is the minimal
    length of an arc of Z/(p-1) containing the support and 0. Returns 0 for a constant function."""
    n = p - 1
    E = [e for e in range(1, len(coeffs)) if coeffs[e] % p]
    if not E:
        return 0
    E.sort()
    best = min(E[-1], n - E[0])
    for i in range(len(E) - 1):
        best = min(best, E[i] + (n - E[i + 1]))
    return best


def ph_phi(ys: np.ndarray, w: int, p: int, e: int, m: int):
    """Pohlig-Hellman as a straight-line program (two inputs y, w; division only by the constant 2).

    z = y^m, zeta = w^m. If p - 1 = 2^e * u with u | m (u odd) and w is a non-residue, zeta generates
    mu_{2^e} and the program outputs L(y) in [0, 2^e) with y^m = zeta^L. Each bit:
      s_i = z_i^(2^(e-1-i)) in {+1, -1};  b_i = (1 - s_i)/2;
      z_{i+1} = z_i * ((1 + s_i) + (1 - s_i) * zeta_i^{-1}) / 2,  zeta_{i+1} = zeta_i^2,
    with zeta_0^{-1} = zeta^(2^e - 1). Output sum_i b_i 2^i. Returns (L array, op count).
    All exponents/constants are N-independent integers; only ring ops and halving are used."""
    inv2 = (p + 1) // 2
    ops = 0
    z = np.array([pow(int(y), m, p) for y in ys], dtype=np.int64)
    ops += 2 * m.bit_length()
    zeta = pow(w, m, p); ops += 2 * m.bit_length()
    zinv = pow(zeta, 2 ** e - 1, p); ops += 2 * e
    L = np.zeros_like(z)
    for i in range(e):
        s = z.copy()
        for _ in range(e - 1 - i):
            s = s * s % p
        ops += e - 1 - i
        b = (1 - s) % p * inv2 % p          # in {0, 1}
        L = (L + b * (2 ** i)) % p
        sel = ((1 + s) % p + (1 - s) % p * zinv % p) % p * inv2 % p
        z = z * sel % p
        zinv = zinv * zinv % p
        ops += 10
    return L, ops


def chain_slp(size: int, rng: random.Random, consts=(1, 2, 3), pmul: float = 0.7, psq: float = 0.6):
    """Post-hoc generator (E1b): each op takes the previous node as left operand. '*' is chosen with
    probability pmul and is a squaring with probability psq, otherwise the right operand is a uniformly
    random earlier node. Degrees grow like 2^{Theta(size)}; the output depends on the whole chain."""
    leaves = ['x'] + list(consts)
    n = len(leaves)
    ops = []
    cur = 0
    for k in range(size):
        if rng.random() < pmul:
            op = '*'
            j = cur if rng.random() < psq else rng.randrange(n + k)
        else:
            op = rng.choice('+-')
            j = rng.randrange(n + k)
            if j == cur:
                j = 0
        ops.append((op, cur, j))
        cur = n + k
    return (leaves, ops)
