"""MODEL (predictions only; floating point allowed here as a fit/model).

Local-factor model for Z_N = #{t <= T : Q_t(N) y-smooth}, Q_t = (m+t)^2 - N.

For a prime l <= y and e >= 0 put g_{l,u}(e) = rho(u - e log l / log y) / rho(u)
(the smoothness chance of a size-y^u integer whose l-part is l^e, relative to
rho(u)).  Under the uniform law on Z_l, Pr(v_l = e) = (1 - 1/l) l^{-e}; under
the law of x^2 - N with x uniform in Z_l it is P_N(e), computed exactly by
Hensel counting.  The local factor is

    lambda_{l,u}(N) = sum_e P_N(e) g(e) / sum_e (1 - 1/l) l^{-e} g(e),

and the model for the conditional mean is E[Z | N] ~ sum_t rho(u_t) prod_{l<=y} lambda_{l,u_t}(N),
u_t = log Q_t / log y.  Only the N-dependent product L(N) is used in the
pre-registered tests (normalised by its sample mean).
"""
from __future__ import annotations

import math
from functools import lru_cache
from typing import Dict, List, Sequence, Tuple

from .sieve import primes_upto

_H = 1e-3


@lru_cache(maxsize=1)
def _rho_table(umax: float = 30.0) -> List[float]:
    n = int(umax / _H) + 2
    k = int(round(1.0 / _H))
    r = [1.0] * n
    for i in range(k + 1, n):
        u = i * _H
        # rho'(u) = -rho(u-1)/u, trapezoid
        u0 = (i - 1) * _H
        r[i] = r[i - 1] - 0.5 * _H * (r[i - 1 - k] / u0 + r[i - k] / u)
    return r


def rho(u: float) -> float:
    if u <= 1.0:
        return 1.0
    t = _rho_table()
    x = u / _H
    i = int(x)
    if i + 1 >= len(t):
        return 0.0
    f = x - i
    return t[i] * (1 - f) + t[i + 1] * f


def vdist_square_minus_N(l: int, N: int, E: int) -> List[float]:
    """Exact Pr(v_l(x^2 - N) = e) for e = 0..E-1 (last entry = Pr(v >= E-1)... see code),
    x uniform mod l^E.  Returns probs for e=0..E (index E collects v >= E)."""
    mod = l ** E
    cnt = [0] * (E + 1)
    if l ** E > 200000:
        raise ValueError("modulus too large for enumeration")
    Nm = N % mod
    for x in range(mod):
        a = (x * x - Nm) % mod
        if a == 0:
            cnt[E] += 1
            continue
        e = 0
        while a % l == 0:
            a //= l
            e += 1
        cnt[e] += 1
    return [c / mod for c in cnt]


def _vdist(l: int, N: int, E: int) -> List[float]:
    """Pr(v_l(x^2-N) = e), e = 0..E (E = tail >= E). Closed form for odd l not dividing N."""
    if l != 2 and N % l != 0:
        chi = 1 if pow(N % l, (l - 1) // 2, l) == 1 else -1
        if chi == -1:
            return [1.0] + [0.0] * E
        # Pr(v >= e) = 2 l^{-e} for e >= 1 (Hensel)
        ge = [1.0] + [2.0 * l ** (-e) for e in range(1, E + 1)]
        return [ge[e] - ge[e + 1] for e in range(E)] + [ge[E]]
    # l = 2 or l | N: enumerate modulo l^E with l^E <= 2^17
    EE = E
    while l ** EE > 131072:
        EE -= 1
    d = vdist_square_minus_N(l, N, EE)
    return d + [0.0] * (E - EE)


def local_factor(l: int, N: int, u: float, logy: float, E: int = 12) -> float:
    d = _vdist(l, N, E)
    r0 = rho(u)
    if r0 <= 0:
        return 1.0
    g = [rho(u - e * math.log(l) / logy) / r0 for e in range(E + 1)]
    num = sum(d[e] * g[e] for e in range(E + 1))
    den = sum((1 - 1 / l) * l ** (-e) * g[e] for e in range(E)) + l ** (-E) * g[E]
    return num / den


def model_LN(N: int, y: int, us: Sequence[float], weights: Sequence[float]) -> float:
    """sum_j w_j prod_l lambda_{l,u_j}(N) / sum_j w_j  (w_j = rho(u_j) * block length)."""
    logy = math.log(y)
    ps = primes_upto(y)
    tot = 0.0
    for u, w in zip(us, weights):
        lp = 0.0
        for l in ps:
            lp += math.log(local_factor(l, N, u, logy))
        tot += w * math.exp(lp)
    return tot / sum(weights)


def u_blocks(n: int, y: int, T: int, nblocks: int = 16) -> Tuple[List[float], List[float]]:
    """Representative u_t and weights rho(u_t)*len over geometric t-blocks (model only)."""
    m = 2.0 ** (n / 2.0) * 1.2247  # ~ sqrt of a typical n-bit N (1.5 * 2^{n-1})
    edges = sorted({max(1, int(round(T ** (j / nblocks)))) for j in range(nblocks + 1)} | {1, T})
    us, ws = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        tm = math.sqrt(a * b)
        Qm = 2 * m * tm + tm * tm
        u = math.log(Qm) / math.log(y)
        us.append(u)
        ws.append(rho(u) * (b - a))
    return us, ws


def expected_Z(n: int, y: int, T: int) -> float:
    us, ws = u_blocks(n, y, T, 64)
    return sum(ws)


def pair_ratio_model(l: int, u: float, logy: float, E: int) -> Dict[int, float]:
    """Unconditional local pair factor R_l(j) = E[w(v) w(v')] / E[w]^2 over N, x uniform
    mod l^E, v = x^2 - N, v' = (x+k)^2 - N, for k with v_l(k) = j (j = 0..E-1).
    w(v) = g(v_l(v)), g(e) = rho(u - e log l/log y)/rho(u); v_l capped at E."""
    mod = l ** E
    r0 = rho(u)
    g = [rho(u - e * math.log(l) / logy) / r0 for e in range(E + 1)]

    def val(a: int) -> int:
        if a == 0:
            return E
        e = 0
        while a % l == 0:
            a //= l
            e += 1
        return e
    vt = [val(a) for a in range(mod)]
    w = [g[vt[a]] for a in range(mod)]
    Ew = sum(w) / mod
    out = {}
    for j in range(E):
        k = l ** j  # representative: only v_l(k) matters (k = l^j * unit, unit absorbed by x -> x*unit? see note)
        s = 0.0
        for N in range(mod):
            for x in range(mod):
                a = (x * x - N) % mod
                b = ((x + k) * (x + k) - N) % mod
                s += w[a] * w[b]
        out[j] = s / (mod * mod) / (Ew * Ew)
    return out


def pair_ratio_avg(l: int, u: float, logy: float, E: int) -> float:
    """Average of R_l over k with l | k (v_l(k) = j has weight (1-1/l) l^{-(j-1)})."""
    R = pair_ratio_model(l, u, logy, E)
    tot, wt = 0.0, 0.0
    for j in range(1, E):
        wj = (1 - 1 / l) * l ** (-(j - 1))
        tot += wj * R[j]
        wt += wj
    # tail j >= E: use R[E-1]
    tail = l ** (-(E - 1))
    tot += tail * R[E - 1]
    wt += tail
    return tot / wt
