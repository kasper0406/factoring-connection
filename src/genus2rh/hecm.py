"""Hyperelliptic curve method (genus 2), OPERATIONAL code.

Inputs: N, public parameters (v, w), a seed. Never uses p, q or group orders.
No primality test or factorisation routine is called on N. ``primes_upto`` sieves the
small primes r <= w used in the exponent (this is public and independent of N).

One curve (Algorithm C in notes/genus2rh_results.md):
  1. draw a, b, f1..f5 uniformly in Z/N; f0 := b^2 - sum_{i>=1} f_i a^i, so (a, b) lies
     on y^2 = f(x);
  2. g = gcd(f5, N): proper -> output; g = N -> give up on this curve;
  3. pass to the monic model F(X) = f5^4 f(X/f5), point (f5 a, f5^2 b);
  4. g = gcd(disc F, N): proper -> output; g = N -> give up;
  5. D = (X - A, B); for each prime r <= w in increasing order, e(r) times: D <- [r]D,
     where e(r) = max{m : r^m <= (isqrt(v) + 2)^4}  (>= the exponent of r in any
     #J(F_p) with p <= v, by Weil's bound #J <= (sqrt p + 1)^4);
  6. a NonUnit exception anywhere yields a proper divisor; otherwise return None.
"""
import random
from math import gcd, isqrt

from .cantor import Jac, NonUnit

__all__ = ["primes_upto", "exponent_schedule", "sample_curve", "one_curve", "hecm"]


def primes_upto(w):
    if w < 2:
        return []
    s = bytearray([1]) * (w + 1)
    s[0] = s[1] = 0
    for i in range(2, isqrt(w) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(w + 1) if s[i]]


def exponent_schedule(v, w):
    """List of primes r <= w, increasing, r repeated e(r) times."""
    bound = (isqrt(v) + 2) ** 4
    sched = []
    for r in primes_upto(w):
        m, x = 0, r
        while x <= bound:
            m, x = m + 1, x * r
        sched.extend([r] * m)
    return sched


def _det_bareiss(M):
    """Exact integer determinant (fraction-free Bareiss)."""
    M = [row[:] for row in M]
    n, sign, prev = len(M), 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for i in range(k + 1, n):
                if M[i][k] != 0:
                    M[k], M[i] = M[i], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def resultant_F_dF(F):
    """Res(F, F') as an integer for integer coefficients F (low-first, degree 5)."""
    dF = [i * F[i] for i in range(1, len(F))]
    m, n = len(F) - 1, len(dF) - 1
    size = m + n
    S = []
    for i in range(n):
        S.append([0] * i + F[::-1] + [0] * (size - m - 1 - i))
    for i in range(m):
        S.append([0] * i + dF[::-1] + [0] * (size - n - 1 - i))
    return _det_bareiss(S)


def sample_curve(N, rng):
    a = rng.randrange(N)
    b = rng.randrange(N)
    f = [0] + [rng.randrange(N) for _ in range(5)]
    f[0] = (b * b - sum(f[i] * pow(a, i, N) for i in range(1, 6))) % N
    return f, a, b


def monic_model(N, f, a, b):
    """F(X) = f5^4 f(X/f5) = X^5 + f4 X^4 + f3 f5 X^3 + f2 f5^2 X^2 + f1 f5^3 X + f0 f5^4."""
    c = f[5]
    F = [(f[i] * pow(c, 4 - i, N)) % N for i in range(5)] + [1]
    return F, (c * a) % N, (c * c * b) % N


def one_curve(N, v, w, rng, sched=None, record=None):
    """Return (divisor or None, stage) with stage in {'f5','disc','main',None,'discard'}."""
    f, a, b = sample_curve(N, rng)
    if record is not None:
        record.update(f=f, a=a, b=b)
    g = gcd(f[5], N)
    if g == N:
        return None, "discard"
    if g != 1:
        return g, "f5"
    F, A, B = monic_model(N, f, a, b)
    g = gcd(resultant_F_dF(F) % N, N)
    if g == N:
        return None, "discard"
    if g != 1:
        return g, "disc"
    J = Jac(N, F)
    D = ((-A % N, 1), (B % N,) if B % N else ())
    if sched is None:
        sched = exponent_schedule(v, w)
    try:
        for r in sched:
            D = J.smul(r, D)
    except NonUnit as e:
        if record is not None:
            record.update(ops=dict(J.ops))
        return e.d, "main"
    if record is not None:
        record.update(final=D, F=F, A=A, B=B, ops=dict(J.ops))
    return None, None


def hecm(N, seed, max_stage=None, curves_per_stage=None):
    """Toy driver with the doubling guess v = 2^j (correctness use only).

    curves_per_stage(v) and w(v) default to tiny values; the asymptotic parameter
    choice w = L_v[1/2, 1], h = L_v[1/2, 1 + eps] is in the notes.
    """
    from math import exp, log, sqrt
    rng = random.Random(seed)
    j = 4
    while max_stage is None or j <= max_stage:
        v = 2 ** j
        lv = log(v)
        w = max(3, int(exp(sqrt(lv * log(lv)))))
        h = curves_per_stage(v) if curves_per_stage else max(4, w)
        sched = exponent_schedule(v, w)
        for _ in range(h):
            d, _stage = one_curve(N, v, w, rng, sched)
            if d is not None:
                return d
        j += 1
    return None
