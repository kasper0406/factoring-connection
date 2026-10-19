"""HARNESS (knows the factors). Instance generators and hit classification.

Never imported by operational modules (polyprod, normbox, copper, smoothcopp).
Primality testing (gmpy2) is used here only to build instances.
"""
from __future__ import annotations

import random
from fractions import Fraction
from math import gcd, isqrt
from typing import Dict, List, Optional, Tuple

import gmpy2


def _rng(tag: str, seed: int) -> random.Random:
    return random.Random(f"specialfam|{tag}|{seed}")


def is_prime(x: int) -> bool:
    return x > 1 and bool(gmpy2.is_prime(x, 50))


def rand_prime_in(rng: random.Random, lo: int, hi: int) -> int:
    while True:
        x = rng.randrange(lo, hi + 1) | 1
        if lo <= x <= hi and is_prime(x):
            return x


def hl_prime(rng: random.Random, kp: int, eta: Fraction) -> Tuple[int, int, int]:
    """Prime p in (2^{kp-1}, 2^kp] with p = a^2 + b^2, 1 <= b <= p^eta (b uniform)."""
    lo, hi = (1 << (kp - 1)) + 1, 1 << kp
    while True:
        bmax = int(gmpy2.iroot(gmpy2.mpz(lo) ** eta.numerator, eta.denominator)[0])
        b = rng.randint(1, max(1, bmax))
        alo = isqrt(max(lo - b * b, 0))
        ahi = isqrt(hi - b * b)
        a = rng.randint(alo, ahi)
        p = a * a + b * b
        if lo <= p <= hi and is_prime(p) and b ** eta.denominator <= p ** eta.numerator:
            return p, a, b


def normbox_instance(n: int, beta: Fraction, eta: Fraction, seed: int) -> Dict:
    """N = p q, p as in hl_prime with kp = round(beta n), q a random prime with q > p."""
    rng = _rng(f"normbox|{n}|{beta}|{eta}", seed)
    kp = max(8, round(beta * n))
    while True:
        p, a, b = hl_prime(rng, kp, eta)
        kq = n - kp
        if kq < kp:
            continue
        if beta == Fraction(1, 2):
            q = rand_prime_in(rng, p + 1, 2 * p - 1)
        else:
            q = rand_prime_in(rng, 1 << (kq - 1), (1 << kq) - 1)
        N = p * q
        if q != p:
            return {"N": N, "p": p, "q": q, "a": a, "b": b, "kp": kp, "n": N.bit_length()}


def control_instance(n: int, seed: int) -> Dict:
    """Balanced semiprime with no imposed structure (control)."""
    rng = _rng(f"control|{n}", seed)
    kp = n // 2
    p = rand_prime_in(rng, (1 << (kp - 1)) + 1, 1 << kp)
    q = rand_prime_in(rng, p + 1, 2 * p - 1)
    return {"N": p * q, "p": p, "q": q, "n": (p * q).bit_length()}


def sumsq_rep(p: int) -> Optional[Tuple[int, int]]:
    """Harness: (a, b) with a >= b >= 1 and a^2 + b^2 = p (Cornacchia via sympy-free search)."""
    if p % 4 != 1:
        return None
    # find sqrt(-1) mod p by trying bases
    for c in range(2, 200):
        t = pow(c, (p - 1) // 4, p)
        if t * t % p == p - 1:
            break
    else:
        return None
    r0, r1 = p, t
    lim = isqrt(p)
    while r1 > lim:
        r0, r1 = r1, r0 % r1
    a = r1
    b = isqrt(p - a * a)
    if a * a + b * b != p:
        return None
    return (max(a, b), min(a, b))


def smooth_gcd_instance(n: int, g_fac: Dict[int, int], signs_p: str, signs_q: str,
                        seed: int) -> Dict:
    """Balanced N = p q with g | p^2 - 1 and g | q^2 - 1.

    g_fac: {ell: e}. signs_*: string over {'+','-'} per prime (sorted), choosing
    p == +-1 mod ell^e independently (mixed CRT sign patterns).
    """
    rng = _rng(f"smooth|{n}|{sorted(g_fac.items())}|{signs_p}|{signs_q}", seed)
    items = sorted(g_fac.items())
    g = 1
    for ell, e in items:
        g *= ell ** e

    def residue(signs: str) -> int:
        r = 0
        for (ell, e), s in zip(items, signs):
            q = ell ** e
            a = 1 if s == "+" else q - 1
            Q = g // q
            r = (r + a * Q * pow(Q, -1, q)) % g
        return r

    rp, rq = residue(signs_p), residue(signs_q)
    kp = n // 2
    lo, hi = (1 << (kp - 1)) + 1, 1 << kp
    if g > hi // 64:
        raise ValueError("g too large for the size")
    while True:
        p = rng.randrange(lo // g, hi // g + 1) * g + rp
        if lo <= p <= hi and is_prime(p):
            break
    while True:
        q = rng.randrange((p + 1) // g, (2 * p - 1) // g + 1) * g + rq
        if p < q < 2 * p and is_prime(q):
            break
    return {"N": p * q, "p": p, "q": q, "g": g, "n": (p * q).bit_length()}


def common_gcd_instance(n: int, gbits: int, seed: int, pm: str = "++") -> Dict:
    """Balanced N = p q with a prime-ish (not smooth) common factor g of p -+ 1, q -+ 1.

    g = 2 * (random prime of gbits-1 bits). pm: signs for p, q ('+' means == 1 mod g).
    """
    rng = _rng(f"common|{n}|{gbits}|{pm}", seed)
    kp = n // 2
    lo, hi = (1 << (kp - 1)) + 1, 1 << kp
    while True:
        g = 2 * rand_prime_in(rng, 1 << (gbits - 2), (1 << (gbits - 1)) - 1)
        if g < hi // 64:
            break
    rp = 1 if pm[0] == "+" else g - 1
    rq = 1 if pm[1] == "+" else g - 1
    while True:
        p = rng.randrange(lo // g, hi // g + 1) * g + rp
        if lo <= p <= hi and is_prime(p):
            break
    for _ in range(100000):
        q = rng.randrange((p + 1) // g, (2 * p - 1) // g + 1) * g + rq
        if p < q < 2 * p and is_prime(q):
            return {"N": p * q, "p": p, "q": q, "g": g, "n": (p * q).bit_length(),
                    "gcd_pm1": gcd(p - 1, q - 1)}
    raise RuntimeError("no q")


def classify_normbox(inst: Dict, res: Dict, eta: Fraction) -> str:
    """STRUCTURAL-PROMISED if the hit is in the band/block containing the promised (a, b);
    STRUCTURAL-OTHER if another leaf of the box is divisible by the found prime;
    TRIAL for a small-prime trial hit."""
    if res.get("how") == "trial":
        return "TRIAL"
    f = res.get("factor")
    if not f:
        return "FAIL"
    k = res.get("band")
    blk = res.get("block")
    m = res.get("m")
    a, b, p = inst["a"], inst["b"], inst["p"]
    if f in (inst["p"], inst["q"]) and blk is not None and m:
        v, bb = blk
        if f == p and bb == b and v * m <= a < (v + 1) * m:
            return "STRUCTURAL-PROMISED"
    return "STRUCTURAL-OTHER"
