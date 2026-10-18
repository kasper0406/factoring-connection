"""Deterministic single-polynomial quadratic sieve (OPERATIONAL part) plus
harness statistics (escape 4: what a derandomization proof would need).

`qs_factor(N, cfg)` sees only N.  It is fully deterministic: the polynomial is
Q(x) = (s + x)^2 - N with s = ceil(sqrt N), x in [-X0, X0); factor base = -1, 2
and odd primes l <= y with (N/l) = 1; sieve with INTEGER weights
(bit lengths) as a candidate filter, then exact trial division; exact F_2
Gaussian elimination (Python-int bitsets) for a left-kernel basis; every kernel
basis vector gives X = prod(s+x_i), Y = prod l^{e_l/2} (mod N) and a gcd.
No randomness, no floating point.  A prime l | N met while building the factor
base is output and classified ACCIDENTAL (trial division) by the harness.

Facts used by the harness (classical, Blum/Rabin): if p == q == 3 (mod 4) and
X^2 == Y^2 (mod N) with gcd(XY, N) = 1, then gcd(X - Y, N) is proper iff the
Jacobi symbol (XY/N) = -1.  The sign map S -> X_S/Y_S (mod p, mod q) is a
homomorphism on the kernel, so all of the 2^k kernel vectors are trivial iff all
k basis vectors are.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional

import gmpy2
import numpy as np

from .avgbound import primes_upto, sqrt_mod_prime


@dataclass(frozen=True)
class QSConfig:
    y: int            # factor-base bound
    X0: int           # initial half-width of the sieve interval
    extra: int = 24   # relations beyond |FB| before stopping
    X0_cap: int = 1 << 22
    block: int = 4096


def default_config(n: int) -> QSConfig:
    """Parameters as a function of n only (frozen in the pre-registration)."""
    table = {40: (1 << 9, 1 << 13), 48: (1 << 10, 1 << 15), 56: (1 << 11, 1 << 16),
             64: (1 << 12, 1 << 18), 72: (1 << 13, 1 << 19), 80: (1 << 14, 1 << 20)}
    if n in table:
        y, x0 = table[n]
    else:
        y = 1 << max(6, n // 5 - 1)
        x0 = 1 << max(10, (n * 5) // 18)
    return QSConfig(y=y, X0=x0)


def _factor_base(N: int, y: int):
    fb = [-1, 2]
    roots = {}
    hit = None
    for l in primes_upto(y):
        l = int(l)
        if l == 2:
            continue
        if N % l == 0:
            if 1 < l < N:
                hit = l
            continue
        if gmpy2.legendre(N, l) == 1:
            fb.append(l)
            roots[l] = sqrt_mod_prime(N % l, l)
    return fb, roots, hit


def _sieve(N: int, s: int, fb, roots, X0: int, block: int, slack: int):
    """Return candidate x values in [-X0, X0) (integer-weight sieve)."""
    L = 2 * X0
    w = np.zeros(L, dtype=np.int32)
    xs0 = -X0
    for l in fb[2:]:
        bl = l.bit_length()
        r = roots[l]
        for t in {r, (l - r) % l}:
            # s + x == t (mod l)  =>  x == t - s
            start = (t - s - xs0) % l
            w[start::l] += bl
    # factor 2: Q(x) = (s+x)^2 - N is even iff s + x == N (mod 2)
    start = (N - s - xs0) % 2
    w[start::2] += 1
    cands = []
    for b0 in range(0, L, block):
        b1 = min(L, b0 + block)
        xa, xb = xs0 + b0, xs0 + b1 - 1
        qmax = max(abs((s + xa) ** 2 - N), abs((s + xb) ** 2 - N), 1)
        thr = qmax.bit_length() - slack
        idx = np.nonzero(w[b0:b1] >= thr)[0]
        cands.extend(int(xs0 + b0 + i) for i in idx)
    return cands


def _trial(Qv: int, fb) -> Optional[List[int]]:
    e = [0] * len(fb)
    if Qv < 0:
        e[0] = 1
        Qv = -Qv
    if Qv == 0:
        return None
    for i in range(1, len(fb)):
        l = fb[i]
        while Qv % l == 0:
            Qv //= l
            e[i] += 1
        if Qv == 1:
            break
    return e if Qv == 1 else None


def _kernel(rows: List[int], nrows: int) -> List[int]:
    """Left kernel basis of the F_2 matrix whose rows are bitmasks `rows`.
    Returns combinations as bitmasks over row indices."""
    piv: Dict[int, tuple] = {}
    kern = []
    for i, v in enumerate(rows):
        comb = 1 << i
        while v:
            hb = v.bit_length() - 1
            if hb in piv:
                pv, pc = piv[hb]
                v ^= pv
                comb ^= pc
            else:
                piv[hb] = (v, comb)
                break
        if v == 0:
            kern.append(comb)
    return kern


def qs_factor(N: int, cfg: QSConfig) -> dict:
    """OPERATIONAL: deterministic QS on N.  Returns a record with the factor (or
    None) and the raw statistics (kernel dimension, per-basis-vector outcomes)."""
    out = {"N": N, "factor": None, "how": None}
    if N % 2 == 0:
        out.update(factor=2, how="even")
        return out
    s = math.isqrt(N - 1) + 1
    if s * s == N:
        out.update(factor=s, how="square")
        return out
    fb, roots, hit = _factor_base(N, cfg.y)
    if hit is not None:
        out.update(factor=hit, how="trial")
        return out
    need = len(fb) + cfg.extra
    slack = 2 * cfg.y.bit_length()
    X0 = cfg.X0
    rels = {}
    while True:
        # candidates in order of |x| (smallest |Q| first); stop at `need`
        for x in sorted(_sieve(N, s, fb, roots, X0, cfg.block, slack), key=lambda t: (abs(t), t)):
            if x in rels:
                continue
            Qv = (s + x) ** 2 - N
            e = _trial(Qv, fb)
            if e is not None:
                rels[x] = e
                if len(rels) >= need:
                    break
        if len(rels) >= need or 2 * X0 > cfg.X0_cap:
            break
        X0 *= 2
    xs = sorted(rels)
    out.update(fb_size=len(fb), n_rels=len(xs), X0_final=X0)
    rows = []
    for x in xs:
        m = 0
        for i, ei in enumerate(rels[x]):
            if ei & 1:
                m |= 1 << i
        rows.append(m)
    kern = _kernel(rows, len(rows))
    out["kdim"] = len(kern)
    res = []
    for comb in kern:
        X, tot = 1, [0] * len(fb)
        for i, x in enumerate(xs):
            if comb >> i & 1:
                X = X * (s + x) % N
                for j, ej in enumerate(rels[x]):
                    tot[j] += ej
        assert all(t % 2 == 0 for t in tot)
        Y = 1
        for j in range(1, len(fb)):
            if tot[j]:
                Y = Y * pow(fb[j], tot[j] // 2, N) % N
        d = math.gcd((X - Y) % N, N)
        jac = int(gmpy2.jacobi(X * Y % N, N)) if math.gcd(X * Y, N) == 1 else 0
        res.append({"nontriv": 1 < d < N, "jacobi": jac})
        if 1 < d < N and out["factor"] is None:
            out.update(factor=d, how="dependency")
    out["deps"] = res
    return out
