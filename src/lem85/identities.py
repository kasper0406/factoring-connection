"""Exact checks of the elementary identities (E2)-(E5) of notes/lem85_results.md (toy sizes)."""
from fractions import Fraction
from math import gcd

from sympy import factorint, isprime, jacobi_symbol, primerange, totient


def largest_pf(n):
    return max(factorint(n)) if n > 1 else 1


def smallest_pf(n):
    return min(factorint(n)) if n > 1 else float("inf")


def is_smooth(n, y):
    return n >= 1 and largest_pf(n) <= y


def smooth_list(x, y):
    return [n for n in range(1, int(x) + 1) if is_smooth(n, y)]


def psi_ap(x, y, q, a):
    return sum(1 for n in smooth_list(x, y) if (n - a) % q == 0)


def psi_coprime(x, y, q):
    return sum(1 for n in smooth_list(x, y) if gcd(n, q) == 1)


def E(x, y, q, a):
    """E(x,y;q,a) = Psi(x,y;q,a) - Psi_q(x,y)/phi(q), exact rational."""
    return Fraction(psi_ap(x, y, q, a)) - Fraction(psi_coprime(x, y, q), int(totient(q)))


def lift_average(x, y, ell, s, r):
    """Return (average of E(x,y;ell*s, r') over the ell-1 lifts r' of r mod s, E(x,y;s,r)/(ell-1))."""
    assert isprime(ell) and ell > y and s % ell != 0 and gcd(r, s) == 1
    t = ell * s
    lifts = [rp for rp in range(t) if rp % s == r % s and rp % ell != 0]
    assert len(lifts) == ell - 1
    avg = sum((E(x, y, t, rp) for rp in lifts), Fraction(0)) / (ell - 1)
    return avg, E(x, y, s, r) / (ell - 1)


def squarefree_divisors(n):
    ps = list(factorint(n)) if n > 1 else []
    out = [(1, 1)]
    for p in ps:
        out += [(d * p, -mu) for d, mu in out]
    return out  # (d, mu(d))


def mobius_identity(x, y, r, q):
    """chi* = Jacobi symbol (./r) (r odd squarefree, conductor r up to sign); chi = chi* times principal mod q.
    Return (sum_{n<=x y-smooth} chi(n), sum_{d | rad(q_y)} mu(d) chi*(d) Psi(x/d,y;chi*))."""
    chi_star = lambda n: jacobi_symbol(n % r, r) if gcd(n, r) == 1 else 0
    lhs = sum(chi_star(n) for n in smooth_list(x, y) if gcd(n, q) == 1)
    q_y = 1
    for p in (factorint(q) if q > 1 else {}):
        if p <= y:
            q_y *= p
    rhs = 0
    for d, mu in squarefree_divisors(q_y):
        rhs += mu * chi_star(d) * sum(chi_star(n) for n in smooth_list(x // d, y))
    return lhs, rhs


def mk_decompositions(n, y, M1):
    """All (m, k) with mk = n, M1 < m <= P+(m) M1 (i.e. m/P+(m) <= M1 < m), P+(m) <= P-(k), k y-smooth."""
    out = []
    for m in range(2, n + 1):
        if n % m:
            continue
        k = n // m
        p = largest_pf(m)
        if m > M1 and m <= p * M1 and p <= smallest_pf(k) and is_smooth(k, y) and p <= y:
            out.append((m, k))
    return out


def sum_inv_phi(X):
    return sum(Fraction(1, int(totient(j))) for j in range(1, X + 1))
