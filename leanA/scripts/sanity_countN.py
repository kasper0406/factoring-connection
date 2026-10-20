#!/usr/bin/env python3
"""Non-vacuity sanity check for the trusted definitions of LeanA/DefsNFS.lean.

A faithful Python transcription of Params.{Lv,Bsm,Xsp,jJ,Yj,Jset,Ssz,Vl,Mb,Hb,Wb,gSet,fields,cands,countN}
and of gbForm / FA (LeanA/Defs.lean), in exact integer arithmetic (floating point only to evaluate
L = exp(sqrt(log N · max(1, log log N))) and the real powers N^μ etc. to the nearest integer, then
floored/ceiled as in Lean; the check below does not depend on the last digit of these values).

It samples random fields A and primes ℓ ∈ J for a split input N and counts candidates, to show that
the filtered count N_ℓ(A) — the object of the internal hypotheses FirstMoment / SecondMoment — is not
identically zero (the hypotheses are not vacuous because of a definitional slip). This is a smoke
test only; it says nothing about the asymptotic claims.

Usage: .venv/bin/python -I leanA/scripts/sanity_countN.py [seed]
"""
import math, random, sys
from fractions import Fraction as Q
import gmpy2

def theta(d, k, mu, eta):
    nu = Q(1, 4 * k)
    hstar = max(2 * nu, mu - (d - 2) * nu, 1 - d * mu + nu)
    return 2 * mu + hstar + eta + (d - 1) * nu, nu, hstar + eta

def gb_coeffs(g0, g1, T, d):
    """Lemma GB: a_d..a_1 digits (valMinAbs residues of T_i * g0^{-i} mod g1), a_0 = T_0."""
    inv = pow(g0, -1, g1)
    a = [0] * (d + 1)
    Ti = T
    for i in range(d, 0, -1):
        r = (Ti % g1) * pow(inv, i, g1) % g1
        if r > g1 // 2:          # valMinAbs: representative in (-g1/2, g1/2]
            r -= g1
        a[i] = r
        num = Ti - r * g0 ** i
        assert num % g1 == 0
        Ti = num // g1
    a[0] = Ti
    return a

def form(a, d, x, y):
    return sum(a[i] * x ** i * y ** (d - i) for i in range(d + 1))

def smooth(n, B, primorial):
    """|n| nonzero and B-smooth (Bernstein-style: n | primorial^(2^e))."""
    n = abs(n)
    if n == 0:
        return False
    g = gmpy2.gcd(n, primorial)
    while g > 1:
        while n % g == 0:
            n //= g
        g = gmpy2.gcd(n, g)
    return n == 1

def main(seed=1):
    rnd = random.Random(seed)
    d, k, mu, eta = 3, 3, Q(7, 24), Q(1, 100)          # the admissible tuple quoted in [pp] §sec:dlog
    th, nu, h = theta(d, k, mu, eta)
    assert th < 1 and nu + eta <= mu < 1 and Q(d + 1, 4) < k
    a_exp = float(th) / math.sqrt(2)
    p = gmpy2.next_prime(rnd.getrandbits(50) | (1 << 49)); q = gmpy2.next_prime(rnd.getrandbits(50) | (1 << 49))
    N = int(p * q)
    lam = math.sqrt(math.log(N) * max(1.0, math.log(math.log(N))))
    L = math.exp(lam)
    B = int(math.floor(L ** a_exp))
    assert N % 2 == 1 and min(p, q) > B                  # SplitInput: odd, no prime factor <= B, not a prime power
    X = math.isqrt(N // 2)
    j = 0
    while 2 ** (k * (j + 1)) <= X:
        j += 1
    Y = 2 ** (j - 1)
    S = math.ceil(L ** (float(th) / (4 * a_exp)))
    M = int(math.floor(N ** float(mu))); H = int(math.floor(N ** float(h)))
    W = math.ceil(N ** float(mu + nu + eta / 2))
    print(f"N={N} (~2^{N.bit_length()}), theta={th}, a={a_exp:.4f}, B={B}, S={S}, Y={Y}, M={M}, H={H}, W={W}")
    primorial = gmpy2.primorial(B)
    hits = trials = 0
    for _ in range(40):
        # field A = ((g0, g1), w, c): uniform on gSet x primes(W, 2W] x C_H
        while True:
            g0 = rnd.randint((M + 1) // 2, M); g1 = rnd.randint((M + 1) // 2, M)
            if math.gcd(g0, g1) == 1:
                break
        while True:
            w = rnd.randint(W + 1, 2 * W)
            if gmpy2.is_prime(w):
                break
        c = [rnd.randint(-H, H) for _ in range(d)]
        while True:
            ell = rnd.randint(Y + 1, 2 * Y)
            if gmpy2.is_prime(ell) and g1 % ell != 0:
                break
        V = math.isqrt(ell * S * S)                       # floor(sqrt(ell) * S)
        acoef = gb_coeffs(g0, g1, N * w, d)
        assert form(acoef, d, g0, g1) == N * w            # Lemma GB identity
        trials += 1
        cnt = 0
        # Lambda_ell(g) ∩ box: b prime in (V/2, V], a ≡ g0 b g1^{-1} (mod ell), |a| <= V
        ginv = pow(g1, -1, ell)
        for b in range(V // 2 + 1, V + 1):
            if not gmpy2.is_prime(b):
                continue
            a0 = g0 * b * ginv % ell
            for aa in range(a0 - ((a0 + V) // ell) * ell, V + 1, ell):
                if aa < -V or math.gcd(aa, b) != 1:
                    continue
                G = g1 * aa - g0 * b
                assert G % ell == 0
                if 32 * abs(G) < M * V:
                    continue
                if not smooth(G // ell, B, primorial):
                    continue
                Fv = form(acoef, d, aa, b) + G * sum(c[i] * aa ** (d - 1 - i) * b ** i for i in range(d))
                if smooth(Fv, B, primorial):
                    cnt += 1
        if N % ell == 0:
            cnt = 0
        hits += cnt > 0
        print(f"  field {trials}: ell={ell}, V={V}, N_ell(A)={cnt}")
    print(f"fraction of sampled (A, ell) with N_ell(A) >= 1: {hits}/{trials}")
    return hits

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
