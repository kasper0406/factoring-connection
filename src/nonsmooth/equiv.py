"""Equivalence map: factoring-equivalent invariants of N, with exact reductions.

Each ``reduce_*`` is OPERATIONAL: input N and the oracle answer (the invariant's value), output a
proper factor of N or None. It never sees p or q. Each ``oracle_*`` is DIAGNOSTIC: it computes
the invariant from the factorisation (harness only).

All reductions here are for squarefree semiprimes N = pq, the case that matters for the
equivalence map; general-N reductions are cited in notes/nonsmooth_results.md.
"""
from __future__ import annotations

import random
from fractions import Fraction
from math import gcd, isqrt

import gmpy2

from src.nonsmooth import forms


# ---------------------------------------------------------------- helpers (operational)

def split_from_s(N, s):
    """Given s = p + q, return p (the smaller root of x^2 - s x + N) if consistent, else None."""
    if s is None:
        return None
    disc = s * s - 4 * N
    if disc < 0:
        return None
    r = isqrt(disc)
    if r * r != disc or (s - r) % 2:
        return None
    p = (s - r) // 2
    if 1 < p < N and N % p == 0:
        return p
    return None


def jacobi(a, n):
    return int(gmpy2.jacobi(a, n))


# ---------------------------------------------------------------- arithmetic functions

def oracle_phi(p, q):
    return (p - 1) * (q - 1)


def reduce_phi(N, phi):
    """phi(N) = N - s + 1."""
    return split_from_s(N, N + 1 - phi)


def oracle_sigma(p, q):
    return (p + 1) * (q + 1)


def reduce_sigma(N, sig):
    """sigma(N) = N + s + 1 (Bach-Miller-Shallit 1986 for general N, randomized)."""
    return split_from_s(N, sig - N - 1)


def oracle_psi(p, q):
    """Dedekind psi = #P^1(Z/N) = [SL2(Z) : Gamma_0(N)]."""
    return (p + 1) * (q + 1)


reduce_psi = reduce_sigma


def oracle_r4(p, q):
    """r_4(N) = #{x in Z^4 : |x|^2 = N} = 8 sigma(N) for odd N (Jacobi)."""
    return 8 * (p + 1) * (q + 1)


def reduce_r4(N, r4):
    if N % 2 == 0 or r4 % 8:
        return None
    return reduce_sigma(N, r4 // 8)


def oracle_pillai(p, q):
    """Pillai's function P(N) = sum_{k=1}^{N} gcd(k, N) = (2p-1)(2q-1).

    Every summand is computable in polylog time by Euclid; the sum has N terms.
    """
    return (2 * p - 1) * (2 * q - 1)


def reduce_pillai(N, P):
    # P = 4N - 2s + 1
    num = 4 * N + 1 - P
    if num % 2:
        return None
    return split_from_s(N, num // 2)


# ---------------------------------------------------------------- modular-curve data

def oracle_genus_X0(p, q):
    """Genus of X_0(N) for N = pq, p, q > 3 distinct primes."""
    N = p * q
    psi = (p + 1) * (q + 1)
    nu2 = (1 + jacobi(-1, p)) * (1 + jacobi(-1, q))
    nu3 = (1 + jacobi(-3, p)) * (1 + jacobi(-3, q))
    nuinf = 4
    g = 1 + Fraction(psi, 12) - Fraction(nu2, 4) - Fraction(nu3, 3) - Fraction(nuinf, 2)
    assert g.denominator == 1
    return int(g)


def reduce_genus_X0(N, g):
    """dim S_2(Gamma_0(N)) = genus(X_0(N)) gives s up to <= 4 public candidates.

    nu2, nu3 take values in {0, 4} when the Jacobi symbol (-1/N) resp. (-3/N) is +1 and 0 when it
    is -1, so at most 2 x 2 candidates; each is tested.
    """
    if gcd(N, 6) != 1:
        return None
    opts2 = (0, 4) if jacobi(-1, N) == 1 else (0,)
    opts3 = (0, 4) if jacobi(-3, N) == 1 else (0,)
    for nu2 in opts2:
        for nu3 in opts3:
            # 12 g = 12 + psi - 3 nu2 - 4 nu3 - 24
            psi = 12 * g - 12 + 3 * nu2 + 4 * nu3 + 24
            p = split_from_s(N, psi - N - 1)
            if p:
                return p
    return None


# ---------------------------------------------------------------- trace formula side

def oracle_hurwitz_sum(p, q):
    """S(N) = sum_{t^2 < 4N} H(4N - t^2), computed by brute force (O(N^{3/2}); tiny N only)."""
    N = p * q
    S = Fraction(0)
    t = 0
    while t * t < 4 * N:
        h = forms.hurwitz_H(4 * N - t * t)
        S += h if t == 0 else 2 * h
        t += 1
    return S


def kronecker_hurwitz_rhs(n_divisors, n):
    """2 sigma(n) - sum_{d | n} min(d, n/d), from an explicit divisor list (diagnostic)."""
    return 2 * sum(n_divisors) - sum(min(d, n // d) for d in n_divisors)


def reduce_hurwitz_sum(N, S):
    """Kronecker-Hurwitz: S(N) = 2 sigma(N) - sum_{d|N} min(d, N/d) = 2N + 2q for N = pq, p < q."""
    S = Fraction(S)
    if S.denominator != 1:
        return None
    S = int(S)
    if S % 2:
        return None
    q = S // 2 - N
    if 1 < q < N and N % q == 0:
        return N // q
    return None


# ---------------------------------------------------------------- modular form coefficients

def tau_table(M):
    """tau(1..M) from Delta = q prod (1-q^n)^24, exact (O(M^2) python ints; diagnostic)."""
    # prod (1 - q^n) via Euler's pentagonal theorem
    e = [0] * (M + 1)
    k = 0
    while True:
        done = True
        for kk in (k, -k) if k else (0,):
            g = kk * (3 * kk - 1) // 2
            if g <= M:
                e[g] += -1 if kk % 2 else 1
                done = False
        if done and k > 0:
            break
        k += 1
    # raise to the 24th power by repeated squaring of truncated series
    def mul(a, b):
        out = [0] * (M + 1)
        for i, ai in enumerate(a):
            if ai:
                lim = M - i
                for j in range(lim + 1):
                    bj = b[j]
                    if bj:
                        out[i + j] += ai * bj
        return out
    e2 = mul(e, e)
    e4 = mul(e2, e2)
    e8 = mul(e4, e4)
    e16 = mul(e8, e8)
    e24 = mul(e16, e8)
    tau = [0] * (M + 1)
    for n in range(1, M + 1):
        tau[n] = e24[n - 1]
    return tau


def oracle_tau_pair(p, q, tau):
    """(tau(N), tau(N^2)) for N = pq, from tau(p), tau(q) (Hecke multiplicativity)."""
    tp, tq = tau[p], tau[q]
    tp2 = tp * tp - p ** 11
    tq2 = tq * tq - q ** 11
    return tp * tq, tp2 * tq2


def reduce_tau_pair(N, A, B):
    """Bach-Charles-type: from A = tau(N), B = tau(N^2).

    C = A^2 + N^11 - B = tau(p)^2 q^11 + tau(q)^2 p^11, and the two summands are the roots of
    z^2 - C z + A^2 N^11. gcd(tau(p)^2 q^11, N) = q unless p | tau(p), and symmetrically, so this
    fails only if BOTH p, q are non-ordinary for Delta (p | tau(p): 2, 3, 5, 7, 2411, ...).
    """
    C = A * A + N ** 11 - B
    disc = C * C - 4 * A * A * N ** 11
    if disc < 0:
        return None
    r = isqrt(disc)
    if r * r != disc:
        return None
    for root in ((C + r) // 2, (C - r) // 2):
        g = gcd(root, N)
        if 1 < g < N:
            return g
    return None


# ---------------------------------------------------------------- Frobenius of a Z/N-algebra

def oracle_quadratic_frobenius(p, q, d):
    """The componentwise Frobenius Phi of A = (Z/N)[y]/(y^2 - d) sends y to eps*y,
    eps = CRT((d/p), (d/q)). Returned as eps in Z/N (diagnostic)."""
    N = p * q
    ep, eq = jacobi(d, p), jacobi(d, q)
    # CRT
    return (ep * q * pow(q, -1, p) + eq * p * pow(p, -1, q)) % N


def reduce_quadratic_frobenius(N, frob_oracle, max_d=10_000):
    """Prop FROB: an oracle for the Frobenius automorphism of quadratic etale Z/N-algebras factors N.

    Pick the least d > 1 with (d/N) = -1 (deterministic search; least such d is small in
    practice, O(log^2 N) under ERH). Then eps^2 = 1 and eps = +1 mod one prime, -1 mod the other.
    Returns (factor, d, number_of_oracle_calls).
    """
    for d in range(2, max_d):
        g = gcd(d, N)
        if 1 < g < N:
            return g, d, 0
        if jacobi(d, N) == -1:
            eps = frob_oracle(d) % N
            g = gcd(eps - 1, N)
            if 1 < g < N:
                return g, d, 1
            return None, d, 1
    return None, None, 0


# ---------------------------------------------------------------- class numbers / class groups

def oracle_class_number_4N(p, q):
    return forms.class_number(-4 * p * q)


def reduce_class_number(N, h, rng, tries=64):
    """Shanks / Schnorr-Lenstra: h(-4N) and class-group arithmetic factor N (randomized).

    Take a random form f of discriminant -4N (from a small prime l with (-4N/l) = 1),
    g = f^{h_odd}; square until the identity; the last non-identity element is ambiguous,
    (a, 0, c), (a, a, c) or (a, b, a), and gcd(a, N) or gcd(c, N) is checked.
    Returns (factor, tries_used).
    """
    D = -4 * N
    h_odd = h
    while h_odd % 2 == 0:
        h_odd //= 2
    ident = forms.identity(D)
    primes = [l for l in range(3, 20000) if gmpy2.is_prime(l)]
    for t in range(1, tries + 1):
        # random form: product of a few small split prime forms with random exponents
        f = ident
        for _ in range(4):
            l = rng.choice(primes)
            g0 = gcd(l, N)
            if 1 < g0 < N:
                return g0, t
            if jacobi(D % l, l) != 1:
                continue
            # b^2 = D mod 4l, b even since D = 0 mod 4
            # find b with b^2 = D (mod 4l), b = 0 mod 2
            r = _sqrt_mod_prime(D % l, l)
            b = r if r % 2 == 0 else r + l   # make b even; then b^2 = D mod 4 too
            b %= 2 * l
            c = (b * b - D) // (4 * l)
            pf = forms.reduce_form((l, b, c))
            f = forms.compose(f, forms.power(pf, rng.randrange(1, 1 << 20)))
        g = forms.power(f, h_odd)
        prev = None
        steps = 0
        while g != ident and steps < 64:
            prev = g
            g = forms.compose(g, g)
            steps += 1
        if prev is None or g != ident:
            continue
        a, b, c = prev
        for x in (a, c, b):
            gg = gcd(x, N)
            if 1 < gg < N:
                return gg, t
    return None, tries


def _sqrt_mod_prime(a, l):
    a %= l
    if a == 0:
        return 0
    if l == 2:
        return a
    # Tonelli-Shanks
    s, e = l - 1, 0
    while s % 2 == 0:
        s //= 2
        e += 1
    z = 2
    while pow(z, (l - 1) // 2, l) != l - 1:
        z += 1
    m, c, t, r = e, pow(z, s, l), pow(a, s, l), pow(a, (s + 1) // 2, l)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % l
            i += 1
        b = pow(c, 1 << (m - i - 1), l)
        m, c, t, r = i, b * b % l, t * b * b % l, r * b % l
    return r


# ---------------------------------------------------------------- invariants hard at primes

def class_number_half_sum(N):
    """For prime p = 3 mod 4, p > 3: h(-p) = (1/(2 - (2/p))) sum_{a < p/2} (a/p) (Dirichlet).

    This interval sum of Jacobi symbols is the reason 'fast interval Jacobi sums' would solve the
    class-number problem at PRIME discriminants (not known to be easier than L[1/2]). O(N) here.
    """
    s = sum(jacobi(a, N) for a in range(1, (N + 1) // 2))
    denom = 2 - jacobi(2, N)
    return Fraction(s, denom)


# ---------------------------------------------------------------- division polynomials (diagnostic)

def division_poly_x(ell, a, b):
    """psi_ell(x) in Z[x] for odd ell on y^2 = x^3 + a x + b (sympy; small ell only).

    Used to check (disc psi_ell / p) = (p / ell): Frobenius permutes the roots of psi_ell with sign
    Legendre(det) = (p/ell), since the sign character of GL2(F_ell) on (F_ell^2 - 0)/{+-1} is the
    Legendre symbol of the determinant (checked for ell = 5, 7, 11, 13 in tests).
    """
    import sympy as sp
    x, y = sp.symbols("x y")
    F = x ** 3 + a * x + b
    psi = {0: sp.Integer(0), 1: sp.Integer(1), 2: 2 * y,
           3: 3 * x ** 4 + 6 * a * x ** 2 + 12 * b * x - a * a,
           4: 4 * y * (x ** 6 + 5 * a * x ** 4 + 20 * b * x ** 3 - 5 * a * a * x ** 2
                       - 4 * a * b * x - 8 * b * b - a ** 3)}

    def red(e):
        e = sp.expand(e)
        P = sp.Poly(e, y)
        out = 0
        for (k,), c in P.terms():
            out += c * F ** (k // 2) * y ** (k % 2)
        return sp.expand(out)

    def get(n):
        if n in psi:
            return psi[n]
        m = n // 2
        if n % 2:
            v = get(m + 2) * get(m) ** 3 - get(m - 1) * get(m + 1) ** 3
        else:
            v = sp.cancel(get(m) * (get(m + 2) * get(m - 1) ** 2
                                     - get(m - 2) * get(m + 1) ** 2) / (2 * y))
        psi[n] = red(v)
        return psi[n]

    g = get(ell)
    return sp.Poly(g, x)


def sign_character_check(ell):
    """True iff the GL2(F_ell) permutation sign on (F_ell^2 - 0)/{+-1} equals Legendre(det)
    on the diagonal generators diag(g, 1) (hence everywhere for ell >= 5, SL2 being perfect)."""
    from itertools import product
    pts, idx = [], {}
    for a_, b_ in product(range(ell), repeat=2):
        if (a_, b_) == (0, 0):
            continue
        key = min((a_, b_), ((-a_) % ell, (-b_) % ell))
        if key not in idx:
            idx[key] = len(pts)
            pts.append(key)
    for g in range(1, ell):
        perm = []
        for (u, v) in pts:
            w = (g * u % ell, v)
            perm.append(idx[min(w, ((-w[0]) % ell, (-w[1]) % ell))])
        seen, s = set(), 1
        for i in range(len(perm)):
            if i in seen:
                continue
            j, L = i, 0
            while j not in seen:
                seen.add(j)
                j = perm[j]
                L += 1
            if L % 2 == 0:
                s = -s
        if s != jacobi(g, ell):
            return False
    return True


# ---------------------------------------------------------------- Prop FROB(b): fixed-space rank

def _polymulmod(u, v, f, m):
    """(u * v) mod (f, m); polys as coefficient lists (low degree first), f monic."""
    n = len(f) - 1
    prod_ = [0] * (len(u) + len(v) - 1)
    for i, ui in enumerate(u):
        if ui:
            for j, vj in enumerate(v):
                prod_[i + j] = (prod_[i + j] + ui * vj) % m
    for k in range(len(prod_) - 1, n - 1, -1):
        c = prod_[k]
        if c:
            for i in range(n + 1):
                prod_[k - n + i] = (prod_[k - n + i] - c * f[i]) % m
    return (prod_[:n] + [0] * n)[:n]


def _polypowmod(base, e, f, m):
    n = len(f) - 1
    r = [1] + [0] * (n - 1)
    while e:
        if e & 1:
            r = _polymulmod(r, base, f, m)
        base = _polymulmod(base, base, f, m)
        e >>= 1
    return r


def oracle_frobenius_matrix(p, q, f):
    """DIAGNOSTIC: matrix of the componentwise Frobenius Phi_B on B = (Z/N)[x]/(f), f monic integer
    coefficient list (low first), in the basis 1, x, ..., x^{n-1}: column j = Phi(x^j) = CRT of
    x^{jp} mod (f, p) and x^{jq} mod (f, q)."""
    N = p * q
    n = len(f) - 1
    cols = []
    xp = _polypowmod([0, 1] + [0] * (n - 2), p, [c % p for c in f], p) if n > 1 else [0]
    xq = _polypowmod([0, 1] + [0] * (n - 2), q, [c % q for c in f], q) if n > 1 else [0]
    ip, iq = pow(q, -1, p), pow(p, -1, q)
    cp = [1] + [0] * (n - 1)
    cq = [1] + [0] * (n - 1)
    for j in range(n):
        cols.append([(cp[i] * q * ip + cq[i] * p * iq) % N for i in range(n)])
        cp = _polymulmod(cp, xp, [c % p for c in f], p)
        cq = _polymulmod(cq, xq, [c % q for c in f], q)
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def reduce_frobenius_fixed_space(N, M):
    """OPERATIONAL (given the oracle matrix M of Phi_B): Gaussian elimination on M - I over Z/N.

    If every pivot is a unit the kernel is free; by Prop FROB(b) it is not when
    (disc f / N) = -1, so a nonzero nonunit entry appears and its gcd with N is returned.
    Returns (factor or None, number of pivots).
    """
    n = len(M)
    A = [[(M[i][j] - (1 if i == j else 0)) % N for j in range(n)] for i in range(n)]
    rank = 0
    row = 0
    for col in range(n):
        piv = None
        for i in range(row, n):
            v = A[i][col]
            if v:
                g = gcd(v, N)
                if 1 < g < N:
                    return g, rank
                piv = i
                break
        if piv is None:
            continue
        A[row], A[piv] = A[piv], A[row]
        inv = pow(A[row][col], -1, N)
        A[row] = [x * inv % N for x in A[row]]
        for i in range(n):
            if i != row and A[i][col]:
                c = A[i][col]
                A[i] = [(A[i][j] - c * A[row][j]) % N for j in range(n)]
        row += 1
        rank += 1
    # remaining entries (all columns processed): any nonzero nonunit left?
    for i in range(n):
        for j in range(n):
            g = gcd(A[i][j], N)
            if 1 < g < N:
                return g, rank
    return None, rank
