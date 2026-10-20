"""Johnson-radius CRT decoders against the lattice-query crux (agent pilist, Task 2 / 3).

Contents
- `log_moment(P, sigma)`: log E[gcd(f(t), P)^sigma] in the random-divisibility model for a linear form
  (Pr[q^j | f(t)] = q^-j, independent over q). Exact, finite sums.
- `rankin_log_bound(logV, P, logS, lam)`: log of the model upper bound
      Gamma <= min_s V S^-s E[d^(s+lam)]
  on the expected relative yield of one decoder call that returns every t of a region of size V with
  d(t) = gcd(f(t), P) >= S (Prop JR). `lam` is the per-unit-log cofactor gain exponent
  log(C u log u)/log B (smoothenum Prop LQ, Dickman step).
- `weighted_log_bound`: the same for weighted (soft-decision, Guruswami-Sahai-Sudan) agreement.
- `count_model(P, logS, dens)`: exact expected fraction of points with d >= S, from local densities.
- `hg_row_decoder(c, P, X, m, t)`: Howgrave-Graham / Boneh list decoder for f(x) = x + c modulo an unknown
  divisor of P (OPERATIONAL in the sense of AGENT_RULES: it sees only public data P, c, X).
- `sieve_row(...)`: exact smooth-part data for one NFS row (HARNESS: used to score outputs).

Exact integer arithmetic everywhere except the log-sums of the bounds (ANALYSIS, rule 1).
"""
from __future__ import annotations

from math import exp, gcd, isqrt, log, sqrt
from typing import Dict, List, Sequence, Tuple

import flint
import numpy as np

# --------------------------------------------------------------------------------- model bounds


def primes_upto(z: int) -> List[int]:
    if z < 2:
        return []
    s = bytearray([1]) * (z + 1)
    s[0] = s[1] = 0
    for i in range(2, isqrt(z) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(z + 1) if s[i]]


def primorial_exps(z: int, power: int = 1) -> Dict[int, int]:
    return {q: power for q in primes_upto(z)}


def log_P(P: Dict[int, int]) -> float:
    return sum(e * log(q) for q, e in P.items())


def log_moment(P: Dict[int, int], sigma: float) -> float:
    """log E[d^sigma], d = gcd(f(t), P), f linear, Pr[v_q(f) >= j] = q^-j independent over q."""
    tot = 0.0
    for q, e in P.items():
        # Pr[min(v, e) = j] = q^-j (1 - 1/q) for j < e, q^-e for j = e
        s = 0.0
        for j in range(e):
            s += q ** (-j) * (1 - 1 / q) * q ** (j * sigma)
        s += q ** (-e) * q ** (e * sigma)
        tot += log(s)
    return tot


def rankin_log_bound(logV: float, P: Dict[int, int], logS: float, lam: float,
                     s_grid: Sequence[float] = None) -> Tuple[float, float]:
    """(log bound, argmin s) for Gamma <= V S^-s E[d^(s+lam)] (any s >= 0; s = 0 gives the trivial V E[d^lam])."""
    if s_grid is None:
        s_grid = [i / 400 for i in range(0, 801)]
    best, arg = None, None
    for s in s_grid:
        val = logV - s * logS + log_moment(P, s + lam)
        if best is None or val < best:
            best, arg = val, s
    return best, arg


def boneh_threshold(logV: float, P: Dict[int, int], D: int = 1, eps: float = 0.0) -> float:
    """Smallest log S inside the proven range log V <= (log S)^2 / (D log P) (Boneh 2002 Thm 3.1/4.1;
    Coppersmith / Howgrave-Graham)."""
    return sqrt((1 + eps) * D * logV * log_P(P))


def jr_log_bound(logV: float, P: Dict[int, int], lam: float, D: int = 1) -> Tuple[float, float]:
    """Worst case over the proven decoding range: the bound is decreasing in S, so take S at the threshold."""
    return rankin_log_bound(logV, P, boneh_threshold(logV, P, D), lam)


def weighted_log_bound(logV: float, P: Dict[int, int], w: Dict[int, float], lam: float,
                       s_max: float = None, iters: int = 200) -> Tuple[float, float]:
    """Soft-decision CRT decoding (Guruswami-Sahai-Sudan 2000): outputs t with
        W(t) = sum_{q^j | gcd} w_q log q >= T = sqrt(log V * sum_q e_q w_q^2 log q)
    (one weight per prime, applied to every power level). Model bound (Rankin, any s >= 0)
        Gamma <= V e^{-sT} prod_q E[e^{s A_q} q^{lam j}],  A_q = j w_q log q.
    The exponent is convex in s (log-MGF plus a linear term), so a ternary search on [0, s_max] finds the
    minimum. Returns (log bound, argmin s)."""
    T = sqrt(logV * sum(e * w.get(q, 0.0) ** 2 * log(q) for q, e in P.items()))
    if s_max is None:
        wmin = min([v for v in w.values() if v > 0] or [1.0])
        s_max = 60.0 / wmin

    def F(s):
        tot = logV - s * T
        for q, e in P.items():
            wq = w.get(q, 0.0)
            lq = log(q)
            # log-sum-exp for stability
            terms = [-j * lq + log(1 - 1 / q) + j * (s * wq * lq + lam * lq) for j in range(e)]
            terms.append(-e * lq + e * (s * wq * lq + lam * lq))
            mx = max(terms)
            tot += mx + log(sum(exp(t - mx) for t in terms))
        return tot
    lo, hi = 0.0, s_max
    for _ in range(iters):
        m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
        if F(m1) <= F(m2):
            hi = m2
        else:
            lo = m1
    s = (lo + hi) / 2
    return min(F(s), F(0.0)), s


def lam_of(u: float, logB: float, C: float = 1.0) -> float:
    """lambda = log(C u log u) / log B (smoothenum Prop LQ step 1)."""
    return log(C * u * log(u)) / logB


# --------------------------------------------------------------------------------- exact count model


def count_model_fraction(P: Dict[int, int], logS: float, dens: Dict[int, List[float]] = None,
                         res: float = 1e-4) -> Tuple[float, float]:
    """Bracket [lo, hi] for Pr[d >= S] under independence over q, with Pr[min(v_q, e_q) = j] = dens[q][j]
    (default: linear form, q^-j (1 - 1/q), q^-e). Dynamic programme on a grid of `res` nats: each log q
    is rounded down (for hi) and up (for lo), so lo <= true value <= hi."""
    items = []
    for q, e in sorted(P.items()):
        pr = dens[q] if (dens is not None and q in dens) else [q ** (-j) * (1 - 1 / q) for j in range(e)] + [q ** (-e)]
        items.append((log(q), pr))
    tot = sum(lq * (len(pr) - 1) for lq, pr in items)
    nb = int(tot / res) + len(items) + 2
    import math
    res_out = {}
    for mode in ("down", "up"):
        dist = np.zeros(nb)
        dist[0] = 1.0
        for lq, pr in items:
            step = int(math.floor(lq / res)) if mode == "down" else int(math.floor(lq / res)) + 1
            nd = np.zeros(nb)
            for j, pj in enumerate(pr):
                if pj == 0.0:
                    continue
                sh = j * step
                if sh < nb:
                    nd[sh:] += pj * dist[: nb - sh]
            dist = nd
        if mode == "down":
            # floor-rounded sum * res <= true sum: [floor-sum >= ceil(thr/res)] implies [true >= thr]
            res_out["lo"] = float(dist[int(math.ceil(logS / res)):].sum())
        else:
            # ceil-rounded sum * res >= true sum: [true >= thr] implies [ceil-sum >= ceil(thr/res)]
            res_out["hi"] = float(dist[int(math.ceil(logS / res)):].sum())
    lo, hi = res_out["lo"], res_out["hi"]
    return lo, hi


# --------------------------------------------------------------------------------- HG / Boneh decoder


def hg_row_decoder(c: int, P: int, X: int, m: int, t: int, nvec: int = 2):
    """All integer x, |x| <= X, that are roots of the first `nvec` reduced vectors of the Howgrave-Graham
    lattice for f(x) = x + c modulo an unknown divisor of P (basis P^{m-i} f^i, i <= m; x^j f^m, 1 <= j < t).
    Returns (roots, l1norms) where l1norms[k] = ||v_k(xX)||_1: every x with |x| <= X and
    gcd(f(x), P)^m > l1norms[k] is a root of v_k (HG lemma), so the call is certified for
    S_cert = l1norms[0]^(1/m)."""
    n = m + t
    polys = []
    f = flint.fmpz_poly([c, 1])
    for i in range(m + 1):
        polys.append(flint.fmpz_poly([P ** (m - i)]) * f ** i)
    for j in range(1, t):
        polys.append(flint.fmpz_poly([0] * j + [1]) * f ** m)
    rows = []
    for g in polys:
        co = [int(a) for a in g.coeffs()] + [0] * (n - g.degree() - 1)
        rows.append([co[i] * X ** i for i in range(n)])
    red = flint.fmpz_mat(rows).lll().tolist()
    roots, norms = set(), []
    for k in range(min(nvec, n)):
        v = [int(a) for a in red[k]]
        norms.append(sum(abs(a) for a in v))
        co = []
        for i in range(n):
            assert v[i] % (X ** i) == 0
            co.append(v[i] // X ** i)
        p = flint.fmpz_poly(co)
        if p.degree() <= 0:
            continue
        for r, _mult in p.roots():
            r = int(r)
            if abs(r) <= X:
                roots.add(r)
    return sorted(roots), norms


# --------------------------------------------------------------------------------- harness: one NFS row


def row_gcd_with_P(Y1: int, Y0: int, b: int, X: int, P: Dict[int, int]) -> Tuple[np.ndarray, np.ndarray]:
    """For a in [-X, X]: exact log gcd(Y1 a + Y0 b, P) (float array) and the exact integer gcd for the
    points whose log gcd is positive is recomputed by the caller when needed. Uses strided sieving."""
    L = 2 * X + 1
    lg = np.zeros(L)
    for q, e in P.items():
        qj = 1
        for j in range(1, e + 1):
            qj *= q
            # solve Y1 a + Y0 b == 0 mod q^j
            if Y1 % q == 0:
                break  # leading coefficient not a unit: handled exactly by the caller (skipped here)
            r = (-Y0 * b * pow(Y1, -1, qj)) % qj
            start = (r + X) % qj  # index of the smallest a >= -X with a == r
            lg[start::qj] += log(q)
    return lg, np.arange(-X, X + 1)


def trial_smooth(n: int, primes: Sequence[int]) -> bool:
    n = abs(n)
    if n == 0:
        return False
    for q in primes:
        if q * q > n:
            break
        while n % q == 0:
            n //= q
    return n == 1 or n <= primes[-1]


def sieve_row_smooth(Y1: int, Y0: int, b: int, X: int, primes: Sequence[int], slack: float = 1.0):
    """HARNESS. Exact list of a in [-X, X] with gcd(a, b) = 1 and Y1 a + Y0 b B-smooth (B = primes[-1]).
    Log-sieve with all prime powers below the value size, then exact trial-division confirmation."""
    L = 2 * X + 1
    lg = np.zeros(L)
    vmax = abs(Y1) * X + abs(Y0) * b
    for q in primes:
        if Y1 % q == 0:
            continue
        lq = log(q)
        qj = q
        inv = None
        while qj <= vmax:
            inv = pow(Y1, -1, qj)
            r = (-Y0 * b * inv) % qj
            lg[(r + X) % qj::qj] += lq
            qj *= q
    a = np.arange(-X, X + 1)
    vals = np.abs(Y1 * a.astype(object) + Y0 * b)
    out = []
    thr = np.log(np.maximum(np.abs(Y1 * a.astype(float) + Y0 * float(b)), 2.0)) - slack
    cand = np.nonzero(lg >= thr)[0]
    for i in cand:
        ai = int(a[i])
        if gcd(ai, b) != 1:
            continue
        if trial_smooth(int(vals[i]), primes):
            out.append(ai)
    return out
