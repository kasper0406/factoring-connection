"""OPERATIONAL decoders: they receive N, the semiprime promise (omega = 2) and oracle answers, nothing else.

Each decoder turns one geometric quantity into p + q (hence the factors) or into a quadratic symbol. They show that
the quantity is as hard as factoring (or as hard as the symbol); they are not factoring algorithms.
No factoring or primality-with-factor routines are imported here (checked by an AST test).
"""
from fractions import Fraction
from itertools import product
from math import gcd, isqrt
from typing import Callable, Dict, List, Optional, Set, Tuple


def _jacobi(a: int, n: int) -> int:
    a %= n
    res = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                res = -res
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            res = -res
        a %= n
    return res if n == 1 else 0


def kron(D: int, n: int) -> int:
    """Kronecker symbol (D/n), n odd positive."""
    return _jacobi(D, n)


def split_from_s(N: int, s: int) -> Optional[Tuple[int, int]]:
    """p + q = s and pq = N  =>  p, q (rule 3: also accept any nontrivial gcd)."""
    disc = s * s - 4 * N
    if disc < 0:
        return None
    r = isqrt(disc)
    if r * r != disc or (s - r) % 2:
        return None
    p, q = (s - r) // 2, (s + r) // 2
    if p > 1 and p * q == N:
        return (p, q)
    return None


# -- closed-form oracles -----------------------------------------------------------------------------------------

def from_index(N: int, mu: int):
    return split_from_s(N, mu - N - 1)


def from_sigma(N: int, sig: int):
    return split_from_s(N, sig - N - 1)


def from_pillai(N: int, P: int):
    # (2p-1)(2q-1) = 4N - 2s + 1
    if (4 * N + 1 - P) % 2:
        return None
    return split_from_s(N, (4 * N + 1 - P) // 2)


def from_genus(N: int, g: int):
    """N = pq odd squarefree: 12 g = 12 + psi - 3 nu2 - 4 nu3 - 24, nu2, nu3 in {0, 4}."""
    for n2, n3 in ((0, 0), (4, 0), (0, 4), (4, 4)):
        mu = 12 * (g + 1) + 3 * n2 + 4 * n3
        res = from_index(N, mu)
        if res:
            return res
    return None


# -- Moebius-sum oracles over small moduli -------------------------------------------------------------------------

def _small_primes(lo: int, count: int, avoid: int) -> List[int]:
    out = []
    r = lo
    while len(out) < count:
        if all(r % d for d in range(2, isqrt(r) + 1)) and avoid % r:
            out.append(r)
        r += 1
    return out


def _crt_merge(cands: Set[Tuple[int, int]], new: Set[int], m: int) -> Set[Tuple[int, int]]:
    """cands = {(x, M)} residues mod M; merge with residues mod m."""
    out = set()
    for (x, M) in cands:
        g = gcd(M, m)
        L = M // g * m
        for y in new:
            if (x - y) % g:
                continue
            # solve z = x mod M, z = y mod m
            t = ((y - x) // g * pow(M // g, -1, m // g)) % (m // g) if m // g > 1 else 0
            out.add(((x + M * t) % L, L))
    return out


def _finish(N: int, cands: Set[Tuple[int, int]]):
    for (x, M) in cands:
        res = split_from_s(N, x)
        if res:
            return res
    return None


def from_necklace(N: int, oracle: Callable[[int, int], int], ks=(2, 3, 5, 7), r_lo: int = 11):
    """oracle(k, r) = sum_{d|N} phi(d) k^{N/d} mod r. Decode s = p + q mod r(r-1) for small primes r; CRT."""
    target = 4 * (isqrt(N) + 2) + 1
    cands = {(0, 1)}
    queries = 0
    for r in _small_primes(r_lo, 200, N * 7 * 5 * 3 * 2):
        kk = [k for k in ks if k % r]
        vals = {k: oracle(k, r) for k in kk}
        queries += len(kk)
        m1 = r - 1
        new = set()
        for a in range(1, r):
            qa = N * pow(a, -1, r) % r
            for b in range(1, m1):
                if gcd(b, m1) != 1:
                    continue
                qb = N * pow(b, -1, m1) % m1
                ok = True
                for k in kk:
                    pred = (pow(k, N % m1, r) + (a - 1) * pow(k, qb, r) + (qa - 1) * pow(k, b, r)
                            + (a - 1) * (qa - 1) * k) % r
                    if pred != vals[k]:
                        ok = False
                        break
                if ok:
                    sr = (a + qa) % r
                    sm = (b + qb) % m1
                    # CRT to mod r*m1
                    z = (sr + r * ((sm - sr) * pow(r, -1, m1) % m1)) % (r * m1)
                    new.add(z)
        cands = _crt_merge(cands, new, r * m1)
        M = next(iter(cands))[1] if cands else 1
        if not cands:
            return None, queries
        if M > target:
            return _finish(N, cands), queries
    return None, queries


def from_catmap(N: int, oracle: Callable[[int, int], int], ts=(3, 4, 5, 6), r_lo: int = 11):
    """oracle(t, r) = primitive N-periodic points of a trace-t toral automorphism, mod r."""
    target = 4 * (isqrt(N) + 2) + 1
    cands = {(0, 1)}
    queries = 0

    def V(t, n, r):
        a, b = 2 % r, t % r
        for bit in bin(n)[2:]:
            if bit == "1":
                a, b = (a * b - t) % r, (b * b - 2) % r
            else:
                a, b = (a * a - 2) % r, (a * b - t) % r
        return a

    for r in _small_primes(r_lo, 200, N * 7 * 5 * 3 * 2):
        tt = [t for t in ts if (t * t - 4) % r]
        if not tt:
            continue
        vals = {t: oracle(t, r) for t in tt}
        queries += len(tt)
        L = (r * r - 1) // 2
        # V_n(t) mod r depends on n mod (r - (t^2-4 / r)), which divides 2L
        FN = {t: (V(t, N, r) - 2) % r for t in tt}
        F1 = {t: (t - 2) % r for t in tt}
        per = {t: r - kron(t * t - 4, r) for t in tt}
        new = set()
        for x in range(1, 2 * L):
            if gcd(x, 2 * L) != 1:
                continue
            y = N * pow(x, -1, 2 * L) % (2 * L)
            ok = True
            for t in tt:
                m = per[t]
                pred = (FN[t] - (V(t, x % m, r) - 2) - (V(t, y % m, r) - 2) + F1[t]) % r
                if pred != vals[t]:
                    ok = False
                    break
            if ok:
                new.add((x + y) % (2 * L))
        cands = _crt_merge(cands, new, 2 * L)
        if not cands:
            return None, queries
        M = next(iter(cands))[1]
        if M > target:
            return _finish(N, cands), queries
    return None, queries


# -- point counts of X_0(N) over F_l: quadratic symbols only ----------------------------------------------------

def _hw(D: int) -> Fraction:
    h = 0
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or (c == a and b < 0) or gcd(gcd(a, abs(b)), c) != 1:
                continue
            h += 1
        a += 1
    if D == -4:
        return Fraction(h, 2)
    if D == -3:
        return Fraction(h, 3)
    return Fraction(h)


def fundamental(D: int) -> int:
    """Fundamental discriminant of the order of discriminant D < 0."""
    d = D
    f = 1
    k = 2
    while k * k <= -d:
        while d % (k * k) == 0 and (d // (k * k)) % 4 in (0, 1):
            d //= k * k
        k += 1
    return d


def x0_count_terms(l: int) -> List[Tuple[int, Fraction]]:
    """[(D0, weight)]: #X_0(N)(F_l) = l + 1 - sigma(l) + 2^omega + (1/2) sum weight * prod_{r|N}(1 + chi_D0(r))."""
    out: Dict[int, Fraction] = {}
    tmax = isqrt(4 * l - 1)
    for t in range(-tmax, tmax + 1):
        D = t * t - 4 * l
        f = 1
        while f * f <= -D:
            if D % (f * f) == 0 and (D // (f * f)) % 4 in (0, 1):
                D0 = fundamental(D)
                out[D0] = out.get(D0, Fraction(0)) + _hw(D // (f * f))
            f += 1
    return sorted(out.items())


def decode_eps_from_x0_counts(N: int, counts: Dict[int, int]):
    """counts[l] = #X_0(N)(F_l). Returns (eps or None, number of consistent patterns, unknown D0 list).

    For N = pq: prod_{r|N}(1 + chi(r)) = 0 if chi(N) = -1, else 2(1 + chi(p)); chi(N) is a Jacobi symbol."""
    eqs = []
    unknown: Set[int] = set()
    for l, cnt in counts.items():
        terms = x0_count_terms(l)
        sig = l + 1
        # cnt = l + 1 - a_l,  a_l = sigma(l) - 4 - (1/2) sum w * prod  ->  sum w * prod = 2 (a_l - sigma(l) + 4)
        a_l = l + 1 - cnt
        rhs = 2 * (a_l - sig + 4)  # = - sum w * prod
        lin = []
        for D0, w in terms:
            if kron(D0, N) == 1:
                lin.append((D0, w))
                unknown.add(D0)
        eqs.append((lin, -rhs))
    unk = sorted(unknown)
    consistent = []
    for pat in product((-1, 1), repeat=len(unk)):
        chi = dict(zip(unk, pat))
        if all(sum(w * 2 * (1 + chi[D0]) for D0, w in lin) == val for lin, val in eqs):
            consistent.append(chi)
    eps_vals = {c.get(-4) for c in consistent}
    eps = eps_vals.pop() if len(eps_vals) == 1 and None not in eps_vals else None
    return eps, len(consistent), unk


# -- amendment A1: point counts over F_{l^2} ------------------------------------------------------------------------

def _terms_n(n: int) -> List[Tuple[int, Fraction]]:
    out: Dict[int, Fraction] = {}
    tmax = isqrt(4 * n - 1)
    for t in range(-tmax, tmax + 1):
        D = t * t - 4 * n
        f = 1
        while f * f <= -D:
            if D % (f * f) == 0 and (D // (f * f)) % 4 in (0, 1):
                D0 = fundamental(D)
                out[D0] = out.get(D0, Fraction(0)) + _hw(D // (f * f))
            f += 1
    return sorted(out.items())


def from_x0_count_sq(N: int, l: int, count: int):
    """count = #X_0(N)(F_{l^2}), N = pq odd with p, q > 2l. Returns (p, q) or None, and the number of candidates.

    count = (l-1) psi/12 + l^2 + 1 - sigma(l^2) + 2 lambda(l^2) + (1/2) sum_D0 w_D0 prod_D0 - l (1 + nu2/4 + nu3/3),
    prod_D0 = prod_{r|N}(1 + chi_D0(r)) = 0 if (D0/N) = -1 else 2 (1 + chi_D0(p)); nu2 = prod_{-4}, nu3 = prod_{-3}."""
    n = l * l
    terms = _terms_n(n)
    sig = 1 + l + l * l
    lam = l + 2
    D0s = sorted({D0 for D0, _ in terms} | {-4, -3})
    unk = [D0 for D0 in D0s if kron(D0, N) == 1]
    tried = 0
    for pat in product((-1, 1), repeat=len(unk)):
        chi = dict(zip(unk, pat))
        prod_ = lambda D0: 0 if kron(D0, N) != 1 else 2 * (1 + chi[D0])
        rest = (Fraction(l * l + 1 - sig + 2 * lam) + Fraction(1, 2) * sum(w * prod_(D0) for D0, w in terms)
                - l * (1 + Fraction(prod_(-4), 4) + Fraction(prod_(-3), 3)))
        psi12 = (count - rest) * 12 / (l - 1)
        tried += 1
        if psi12.denominator != 1:
            continue
        res = from_index(N, int(psi12))
        if res:
            return res, tried
    return None, tried
