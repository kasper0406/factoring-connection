"""DIAGNOSTIC: exact leak probabilities of primality tests from the factorisation.

Uses p, q and factorisations of group orders (sympy.factorint on p-1, p+1, ...).
Never imported by operational code.  Theorem numbers refer to
notes/primality_analysis.md.

Model.  The base is a uniform element x of a group G(Z/N) = prod_i G_i with each
G_i cyclic of order n_i (F_p^* for Fermat/MR/SS, the norm-one torus T_p of order
p - (D/p) for Lucas).  A *check* (E, zeta) compares x^E with zeta in {+1, -1};
its gcd observable gcd(x^E - zeta, N) equals the product of the primes where
x_i^E = zeta.

Order-pattern lemma (Thm 1).  The indicator [x_i^E = zeta] depends only on
d_i = ord(x_i): [d | E] for zeta = 1, [d | 2E and d not| E] for zeta = -1.  The d_i
are independent with P(d_i = d) = phi(d)/n_i.  Hence a family of checks reveals a
proper divisor iff the indicator vectors of the primes are not all equal, and
    P(leak) = 1 - sum_v prod_i P_i(v).
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd
from typing import Dict, Iterable, List, Sequence, Tuple

import sympy

INF = 10 ** 9  # tau = "never"


def v2(x: int) -> int:
    return (abs(x) & -abs(x)).bit_length() - 1


def odd(x: int) -> int:
    return abs(x) >> v2(x)


_DIV_CACHE: Dict[int, List[Tuple[int, int]]] = {}


def divisors_with_phi(n: int) -> List[Tuple[int, int]]:
    if n not in _DIV_CACHE:
        ds = sympy.divisors(n)
        _DIV_CACHE[n] = [(int(d), int(sympy.totient(d))) for d in ds]
    return _DIV_CACHE[n]


# ------------------------------------------------------------ generic engine
def indicator(d: int, E: int, zeta: int) -> int:
    if zeta == 1:
        return int(E % d == 0)
    return int((2 * E) % d == 0 and E % d != 0)


def indicator_dist(n: int, checks: Sequence[Tuple[int, int]], exclude_identity: bool = False
                   ) -> Dict[Tuple[int, ...], Fraction]:
    """Law of the indicator vector ([x^E = zeta])_checks for x uniform in C_n
    (or in C_n minus the identity)."""
    tot = n - (1 if exclude_identity else 0)
    out: Dict[Tuple[int, ...], int] = {}
    for d, ph in divisors_with_phi(n):
        cnt = ph - (1 if (exclude_identity and d == 1) else 0)
        if cnt == 0:
            continue
        v = tuple(indicator(d, E, z) for E, z in checks)
        out[v] = out.get(v, 0) + cnt
    return {v: Fraction(c, tot) for v, c in out.items()}


def leak_prob(dists: Sequence[Dict[Tuple[int, ...], Fraction]]) -> Fraction:
    """P(indicator vectors not all equal) for independent components."""
    keys = set(dists[0])
    for D in dists[1:]:
        keys &= set(D)
    same = Fraction(0)
    for v in keys:
        pr = Fraction(1)
        for D in dists:
            pr *= D[v]
        same += pr
    return 1 - same


def tau_dist(n: int, e: int, exclude_identity: bool = False) -> Dict[int, Fraction]:
    """tau(x) = least i with x^(m 2^i) = 1 (e = 2^t m), INF if odd(ord x) does not divide m.
    Unbounded in i (tau can exceed t)."""
    m = odd(e)
    tot = n - (1 if exclude_identity else 0)
    out: Dict[int, int] = {}
    for d, ph in divisors_with_phi(n):
        cnt = ph - (1 if (exclude_identity and d == 1) else 0)
        if cnt == 0:
            continue
        tau = v2(d) if m % odd(d) == 0 else INF
        out[tau] = out.get(tau, 0) + cnt
    return {k: Fraction(c, tot) for k, c in out.items()}


def ladder_events(taus: Sequence[Dict[int, Fraction]], t: int) -> Dict[str, Fraction]:
    """Exact probabilities of the MR-ladder events for independent tau laws.

    sqrt   : max tau <= t and taus not all equal  (classical nontrivial sqrt of 1)
    single : some but not all tau <= t             (Fermat-type split: x^e = 1 mod some primes only)
    prefix : sqrt + single = some gcd(b_i - 1, N) proper, i <= t
    pm     : some gcd(b_i -+ 1, N) proper, i <= t
    liar   : all taus equal and <= t               (strong liar)
    liar_level[j] : all taus equal to j
    fermat_pass : max tau <= t
    """
    import itertools
    r = len(taus)
    res = {"sqrt": Fraction(0), "single": Fraction(0), "pm": Fraction(0), "liar": Fraction(0),
           "fermat_pass": Fraction(0)}
    levels: Dict[int, Fraction] = {}
    items = [list(T.items()) for T in taus]
    for combo in itertools.product(*items):
        vals = [c[0] for c in combo]
        pr = Fraction(1)
        for c in combo:
            pr *= c[1]
        le = [v <= t for v in vals]
        if all(le):
            res["fermat_pass"] += pr
            if len(set(vals)) == 1:
                res["liar"] += pr
                levels[vals[0]] = levels.get(vals[0], 0) + pr
            else:
                res["sqrt"] += pr
        elif any(le):
            res["single"] += pr
        # pm: exists i in [0,t]: indicator sets for (b_i = 1) or (b_i = -1) proper
        leak = False
        for i in range(t + 1):
            s1 = [v <= i for v in vals]
            sm = [v == i + 1 for v in vals]
            if (any(s1) and not all(s1)) or (any(sm) and not all(sm)):
                leak = True
                break
        if leak:
            res["pm"] += pr
    res["prefix"] = res["sqrt"] + res["single"]
    res["liar_level"] = levels
    return res


# ----------------------------------------------------------- per-test exact
def mr_exact(primes: Sequence[int], exclude_identity: bool = False) -> Dict[str, Fraction]:
    """Miller-Rabin on squarefree N = prod primes, base uniform in (Z/N)^*."""
    N = 1
    for p in primes:
        N *= p
    e = N - 1
    t = v2(e)
    return ladder_events([tau_dist(p - 1, e) for p in primes], t)


def fermat_exact(primes: Sequence[int], chain: bool = False) -> Fraction:
    N = 1
    for p in primes:
        N *= p
    e = N - 1
    if chain:
        checks, pref = [], 0
        for bit in bin(e)[2:]:
            pref = 2 * pref + (bit == "1")
            checks.append((pref, 1))
    else:
        checks = [(e, 1)]
    return leak_prob([indicator_dist(p - 1, checks) for p in primes])


def ss_exact(primes: Sequence[int]) -> Fraction:
    N = 1
    for p in primes:
        N *= p
    E = (N - 1) // 2
    return leak_prob([indicator_dist(p - 1, [(E, 1), (E, -1)]) for p in primes])


def lucas_exact(primes: Sequence[int], D: int, exclude_identity: bool = True) -> Dict[str, Fraction]:
    """Strong Lucas ladder on N+1 in the norm-one tori T_p (order p - (D/p)).
    exclude_identity=True is the law induced by uniform P with Q = (P^2 - D)/4 a unit."""
    N = 1
    for p in primes:
        N *= p
    e = N + 1
    t = v2(e)
    taus = [tau_dist(p - int(sympy.legendre_symbol(D % p, p)), e, exclude_identity) for p in primes]
    return ladder_events(taus, t)


def frobenius_exact(p: int, q: int, Delta: int) -> Fraction:
    """Structural leak of the quadratic Frobenius condition, alpha uniform in R^*."""
    def kappa(pp, qq):
        if int(sympy.legendre_symbol(Delta % pp, pp)) == -1:
            return Fraction(gcd(qq - 1, pp * pp - 1), pp * pp - 1)
        return Fraction(gcd(qq * qq - 1, pp - 1), (pp - 1) ** 2)
    kp, kq = kappa(p, q), kappa(q, p)
    return kp * (1 - kq) + kq * (1 - kp)


# ------------------------------------------------------------ closed forms
def closed_ladder(n_p: int, n_q: int, e: int) -> Dict[str, Fraction]:
    """Thm 3 closed forms for two cyclic components of orders n_p, n_q, exponent e = 2^t m.

    D_i = gcd(e, n_i), u_i = min(v2 n_i, t), u = min(u_p, u_q), kappa_i = D_i / n_i:
      single = kappa_p (1-kappa_q) + kappa_q (1-kappa_p)
      sqrt   = kappa_p kappa_q [1 - (4^u + 2) / (3 * 2^(u_p+u_q))]
      liar   = odd(D_p) odd(D_q) (4^u + 2) / (3 n_p n_q)
    """
    t = v2(e)
    Dp, Dq = gcd(e, n_p), gcd(e, n_q)
    up, uq = min(v2(n_p), t), min(v2(n_q), t)
    u = min(up, uq)
    kp, kq = Fraction(Dp, n_p), Fraction(Dq, n_q)
    single = kp * (1 - kq) + kq * (1 - kp)
    sqrt = kp * kq * (1 - Fraction(4 ** u + 2, 3 * 2 ** (up + uq)))
    liar = Fraction(odd(Dp) * odd(Dq) * (4 ** u + 2), 3 * n_p * n_q)
    return {"single": single, "sqrt": sqrt, "prefix": single + sqrt, "liar": liar,
            "fermat_pass": kp * kq}


def closed_mr_semiprime(p: int, q: int) -> Dict[str, Fraction]:
    """Cor 3a: N = pq, g = gcd(p-1, q-1) = gcd(N-1, p-1) = gcd(N-1, q-1), s_< = min(v2(p-1), v2(q-1)):
      fermat split = g/(p-1) + g/(q-1) - 2 g^2/phi
      sqrt         = (2/3)(1 - 4^(-s_<)) g^2/phi
      liar         = h^2 (4^(s_<) + 2) / (3 phi),  h = odd(g)   (Monier)
    """
    g = gcd(p - 1, q - 1)
    phi = (p - 1) * (q - 1)
    s = min(v2(p - 1), v2(q - 1))
    h = odd(g)
    single = Fraction(g, p - 1) + Fraction(g, q - 1) - Fraction(2 * g * g, phi)
    sqrt = Fraction(2, 3) * (1 - Fraction(1, 4 ** s)) * Fraction(g * g, phi)
    liar = Fraction(h * h * (4 ** s + 2), 3 * phi)
    return {"single": single, "sqrt": sqrt, "prefix": single + sqrt, "liar": liar,
            "sqrt_minus1": Fraction(h * h * (4 ** s - 4), 3 * phi) if s >= 2 else Fraction(0)}


def closed_ss_semiprime(p: int, q: int) -> Fraction:
    """Cor 3c: e = (N-1)/2.  d_i = gcd(N-1, p_i-1) = g, d'_i = gcd((N-1)/2, p_i-1):
    split = 1 - [d'_p d'_q + (g-d'_p)(g-d'_q)]/phi - (1 - g/(p-1))(1 - g/(q-1))."""
    N = p * q
    g = gcd(p - 1, q - 1)
    dp, dq = gcd((N - 1) // 2, p - 1), gcd((N - 1) // 2, q - 1)
    phi = (p - 1) * (q - 1)
    return 1 - Fraction(dp * dq + (g - dp) * (g - dq), phi) - (1 - Fraction(g, p - 1)) * (1 - Fraction(g, q - 1))


def torus_exponent(p: int, q: int, dp: int, dq: int) -> int:
    """Thm 5: gamma^N = gamma^(dq*p + dp*q - dp*dq) on T_N = T_p x T_q, |T_i| = p_i - d_i."""
    return dq * p + dp * q - dp * dq
