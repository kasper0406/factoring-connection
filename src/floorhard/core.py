"""floorhard: reductions FROM exact polynomial floor sums (hardness side of Theorem L_j).

All oracles here are BRUTE FORCE (diagnostic/harness): they test identities and reductions, never
costs. Model cost = number of oracle calls.

Conventions. A floor-sum query is (f, D, h1, h2) with f an integer polynomial (coefficient list,
f[i] is the coefficient of h^i) and D >= 1: F = sum_{h1 <= h < h2} floor(f(h) / D).
Window identity (W): #{h in [h1,h2): f(h) mod D in [u, v)} = F(f - u) - F(f - v), 0 <= u < v <= D.
"""
from __future__ import annotations

import random
from math import gcd, isqrt

import gmpy2
import numpy as np

from src.nonsmooth import forms
from src.nonsmooth.equiv import reduce_class_number


# ----------------------------------------------------------------------------- basic oracles

def poly_eval(f, h):
    r = 0
    for c in reversed(f):
        r = r * h + c
    return r


def floor_sum_brute(f, D, h1, h2):
    """sum_{h1<=h<h2} floor(f(h)/D), exact big-int brute force (diagnostic)."""
    return sum(poly_eval(f, h) // D for h in range(h1, h2))


def floor_sum_np(f, D, h1, h2):
    """Same as floor_sum_brute for small integer polynomials (all |values| < 2^62), numpy int64."""
    h = np.arange(h1, h2, dtype=np.int64)
    v = np.zeros_like(h)
    for c in reversed(f):
        v = v * h + int(c)
    return int((v // int(D)).sum())


def shift(f, s):
    g = list(f)
    g[0] -= s
    return g


def window_count(oracle, f, D, h1, h2, u, v):
    """Identity (W): two oracle calls."""
    return oracle(shift(f, u), D, h1, h2) - oracle(shift(f, v), D, h1, h2)


def jacobi(a, n):
    return int(gmpy2.jacobi(a % n, n))


# ----------------------------------------------------------------------------- (A) full period

def squarefree_odd(m):
    if m % 2 == 0:
        return False
    d = 3
    while d * d <= m:
        if m % (d * d) == 0:
            return False
        d += 2
    return True


def divisors(m):
    ds = []
    for d in range(1, isqrt(m) + 1):
        if m % d == 0:
            ds.append(d)
            if d * d != m:
                ds.append(m // d)
    return sorted(ds)


def hstar(d):
    """h*(-d) = h(-d)*2/w: 1/3 for d = 3 (returned as Fraction-free pair), else h(-d)."""
    from fractions import Fraction
    if d == 3:
        return Fraction(1, 3)
    return Fraction(forms.class_number(-d))


def full_period_rhs(m, a):
    """Prop A right side: a(m-1)(2m-1)/6 - (m-1)/2 + sum_{d|m, d=3 mod 4} (a/d) h*(-d)."""
    from fractions import Fraction
    s = Fraction(a * (m - 1) * (2 * m - 1), 6) - Fraction(m - 1, 2)
    for d in divisors(m):
        if d > 1 and d % 4 == 3:
            s += jacobi(a, d) * hstar(d)
    return s


def full_period_brute(m, a):
    """sum_{h<m} floor(a h^2/m), exact via numpy int64 (needs a*m^2 < 2^63)."""
    assert a * m * m < (1 << 62)
    h = np.arange(m, dtype=np.int64)
    return int(((a * h * h) // m).sum())


def class_part(m, a, oracle=full_period_brute):
    """Oracle value minus the elementary part: sum_{d|m, d=3(4)} (a/d) h*(-d) (exact rational)."""
    from fractions import Fraction
    return Fraction(oracle(m, a)) - Fraction(a * (m - 1) * (2 * m - 1), 6) + Fraction(m - 1, 2)


def decode_class_number_semiprime(m, rng, oracle=full_period_brute, max_calls=40, amax=None):
    """m = pq odd squarefree, 3 does not divide m. Returns (M, h(-M), calls) with M = m (m = 3 mod 4)
    or M = 3m (m = 1 mod 4). Uses multipliers a with (a/M) = 1 (and (a/3) = 1 if M = 3m), so that
    the class part is h(-M) + (const) +- X with an unknown sign; two distinct values average to the
    constant part. Operational on m: uses only m, Jacobi symbols and oracle calls.
    """
    from fractions import Fraction
    M = m if m % 4 == 3 else 3 * m
    seen = set()
    calls = 0
    while calls < max_calls:
        a = rng.randrange(2, amax) if amax else rng.randrange(1, M)   # theorem: uniform mod M
        if gcd(a, M) != 1 or jacobi(a, M) != 1:
            continue
        if M != m and jacobi(a, 3) != 1:
            continue
        calls += 1
        v = class_part(M, a, oracle)
        if M != m:
            v -= Fraction(1, 3)          # (a/3) h*(-3) with (a/3) = 1
        seen.add(v)
        if len(seen) == 2:
            v1, v2 = sorted(seen)
            hM = (v1 + v2) / 2
            assert hM.denominator == 1
            return M, int(hM), calls
    return M, None, calls


def factor_from_class_number(m, M, hM, rng):
    """Ambiguous forms (Shanks / Schnorr-Lenstra) on disc -4M with h(-4M) = h(-M)(2 - (-M/2))."""
    chi2 = 1 if (-M) % 8 == 1 else -1
    h4 = hM * (2 - chi2)
    for _ in range(8):
        g, _t = reduce_class_number(M, h4, rng)
        if g is None:
            continue
        g = gcd(g, m)
        if 1 < g < m:
            return g
    return None


# ----------------------------------------------------------------------------- (B) #SubsetSum

def subset_sum_count(w, T):
    cnt = {0: 1}
    for x in w:
        new = dict(cnt)
        for s, c in cnt.items():
            new[s + x] = new.get(s + x, 0) + c
        cnt = new
    return cnt.get(T, 0)


def subset_sum_to_qc(w, T, deterministic=False):
    """Parsimonious reduction (Prop B). Returns dict with a, b, primes, theta, windows.

    #{x in windows : x^2 = a (mod b)} = #{delta in {0,1}^k : sum delta_i w_i = T}.
    S = sum w, M = 4S + 4, primes p_i > max(kM, 4M) distinct odd, b = prod p_i,
    theta_i = (b/p_i) u_i with u_i = round(p_i w_i / M); a = (sum theta_i)^2 mod b.
    Signs eps = 2 delta - 1; root x = sum eps_i theta_i mod b; x/b = sum eps_i w_i/M + eta mod 1,
    |eta| < 1/(2M). Target T' = 2T - S. Window: x/b in [T'/M - 1/(2M), T'/M + 1/(2M)) mod 1.
    """
    k = len(w)
    S = sum(w)
    assert k >= 1 and min(w) >= 1
    if not 0 <= T <= S:
        return None                            # count is 0 (referee fix: aliasing mod M otherwise)
    M = 4 * S + 4
    P = max(k * M, 4 * M) + 1
    primes = []                                # moduli q_i (primes, or prime powers if deterministic)
    if deterministic:
        # codex floorhard_2 repair: first k odd primes, each raised to a power > 8kM (poly-time, no
        # prime search); u_i adjusted by <= 1 to be a unit, |eta| < 3/(16M).
        ell = 3
        while len(primes) < k:
            q = ell
            while q <= 8 * k * M:
                q *= ell
            primes.append(q)
            ell = int(gmpy2.next_prime(ell))   # small primes (<= k-th prime): trial division is poly
    else:
        q = P
        while len(primes) < k:
            q = int(gmpy2.next_prime(q))
            primes.append(q)
    b = 1
    for p in primes:
        b *= p
    theta = []
    for wi, p in zip(w, primes):
        u = (2 * p * wi + M) // (2 * M)        # round(p w / M), in [1, p/4]
        if deterministic:
            base = next(l for l in range(3, p + 1, 2) if p % l == 0)
            if u % base == 0:
                u += 1
        assert 1 <= u < p
        theta.append((b // p) * u)
    a = pow(sum(theta), 2, b)
    Tp = 2 * T - S
    # window in units of b: [b(2T' - 1)/(2M), b(2T'+1)/(2M)) mod b, as integer x-ranges
    lo = (b * (2 * Tp - 1)) // (2 * M) + 1   # x/b >= (2T'-1)/(2M) <=> x > b(2T'-1)/(2M) (never equal, see proof)
    hi = (b * (2 * Tp + 1)) // (2 * M) + 1
    q0 = lo // b
    if hi <= (q0 + 1) * b:
        windows = [(lo - q0 * b, hi - q0 * b)]
    else:
        windows = [(lo - q0 * b, b), (0, hi - (q0 + 1) * b)]
    return dict(a=a, b=b, primes=primes, theta=theta, windows=windows, M=M)


def qc_count_crt(inst):
    """Diagnostic: enumerate the 2^k CRT roots and count those in the windows."""
    if inst is None:
        return 0
    b, theta = inst["b"], inst["theta"]
    k = len(theta)
    cnt = 0
    for mask in range(1 << k):
        x = sum((t if (mask >> i) & 1 else -t) for i, t in enumerate(theta)) % b
        assert (x * x - inst["a"]) % b == 0
        if any(lo <= x < hi for lo, hi in inst["windows"]):
            cnt += 1
    return cnt


def qc_count_floorsums(inst, oracle=floor_sum_brute):
    """The count through identity (W) with f(x) = x^2 - a, D = b, [u, v) = [0, 1): 2 calls/window."""
    if inst is None:
        return 0
    f = [-inst["a"], 0, 1]
    return sum(window_count(oracle, f, inst["b"], lo, hi, 0, 1) for lo, hi in inst["windows"])


# ----------------------------------------------------------------------------- (C) Taylor class

def taylor_numerator(N, x, j):
    """P_j(h) = sum_{i<=j} (-1)^i N h^i / x^{i+1} = f(h)/x^{j+1}, f integer polynomial."""
    return [(-1) ** i * N * x ** (j - i) for i in range(j + 1)], x ** (j + 1)


def ladder_block_length(N, x, j, c=(1, 2)):
    B = int(gmpy2.iroot(N, j + 2)[0])
    if B ** (j + 2) < N:
        B += 1
    return (c[0] * x) // (c[1] * B)


def embedding_t(A, b, cp, j, extra=1):
    """Smallest t (doubling) with t > 8 A c'^3 / b, t >= 2c'/b and ladder H(bt) >= c' for N = A b^2 t^3."""
    t = max(1, (8 * A * cp ** 3) // b + 1, (2 * cp) // b + 1) * extra
    while ladder_block_length(A * b * b * t ** 3, b * t, j) < cp:
        t *= 2
    return t


class TaylorEmbedding:
    """Prop C: N = A b^2 t^3 + e, x = bt; P_j = Taylor of N/(x+h).

    For 0 <= h < c': P_j(h) = Abt^2 - Ath + A h^2/b + R(h) + e*(...), with |R + e-part - e/x| < 1/(4b).
    Window for residue r: {P_j(h) - e/x + 1/(4b)} in [r/b, r/b + 1/(2b)).
    Implemented exactly over the common denominator 4 b x^{j+1}.
    """

    def __init__(self, A, b, cp, j, t=None, e=0):
        self.A, self.b, self.cp, self.j = A, b, cp, j
        self.t = t if t is not None else embedding_t(A, b, cp, j)
        self.x = b * self.t
        self.e = e
        self.N = A * b * b * self.t ** 3 + e
        f, D = taylor_numerator(self.N, self.x, j)
        self.H = ladder_block_length(self.N, self.x, j)
        # g/E = P_j - e/x + 1/(4b), E = 4 b x^{j+1}
        self.E = 4 * b * D
        g = [4 * b * c for c in f]
        g[0] += -4 * b * e * x_pow(self.x, j) + D
        self.g = g

    def window_count(self, oracle, r, h1, h2):
        """#{h in [h1,h2): {g/E} in [r/b, r/b + 1/(2b))} by identity (W): 2 oracle calls."""
        u = (r * self.E) // self.b
        v = ((2 * r + 1) * self.E) // (2 * self.b)
        # E divisible by 4b so these are exact
        assert u * self.b == r * self.E and v * 2 * self.b == (2 * r + 1) * self.E
        return window_count(oracle, self.g, self.E, h1, h2, u, v)


def x_pow(x, j):
    return x ** j


class PrefixOracle:
    """Brute-force oracle for a fixed polynomial family on [0, L): caches floor values per shift.
    Counts model calls. Diagnostic only."""

    def __init__(self, L):
        self.L = L
        self.cache = {}
        self.calls = 0

    def __call__(self, f, D, h1, h2):
        self.calls += 1
        key = (tuple(f), D)
        if key not in self.cache:
            pre = [0]
            for h in range(self.L):
                pre.append(pre[-1] + poly_eval(f, h) // D)
            self.cache[key] = pre
        pre = self.cache[key]
        return pre[h2] - pre[h1]


def factor_by_sqrt1_bisection(m, count):
    """count(lo, hi) = #{h in [lo,hi): h^2 = 1 mod m}. Finds a nontrivial square root of 1 by bisection
    on [2, m-1) and returns gcd(h - 1, m). Requires m odd with >= 2 distinct prime factors."""
    lo, hi = 2, m - 1
    if count(lo, hi) == 0:
        return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if count(lo, mid) > 0:
            hi = mid
        else:
            lo = mid
    g = gcd(lo - 1, m)
    return g if 1 < g < m else None


def factor_via_taylor_class(m, j=2):
    """Factor m with O(log m) window counts of the Taylor class at N = m^2 t^3, x = mt (A = 1, b = m).
    Returns (factor, oracle_calls, N)."""
    emb = TaylorEmbedding(1, m, m, j)
    orc = PrefixOracle(m)
    g = factor_by_sqrt1_bisection(m, lambda lo, hi: emb.window_count(orc, 1, lo, hi))
    return g, orc.calls, emb.N


def find_prime_embedding(A, b, cp, j, rng, tries=10000):
    """Prime variant: N = A b^2 t^3 + e prime with 0 <= e < Emax = b t^2/(16 c'), e random.
    t is increased so that Emax is large; extra conditions from the proof are re-checked exactly."""
    t = embedding_t(A, b, cp, j, extra=4)
    Emax = b * t * t // (16 * cp)
    for _ in range(tries):
        e = rng.randrange(0, Emax)
        N = A * b * b * t ** 3 + e
        if gmpy2.is_prime(N):
            return TaylorEmbedding(A, b, cp, j, t=t, e=e)
    return None
