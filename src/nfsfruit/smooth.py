"""OPERATIONAL: factor bases and exact smoothness tests (no floating point)."""
from __future__ import annotations
from math import gcd, isqrt
from typing import Dict, List, Optional
import gmpy2

def primes_upto(B: int) -> List[int]:
    s = bytearray([1]) * (B + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, isqrt(B) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(B + 1) if s[i]]

class FactorBase:
    def __init__(self, B: int):
        self.B = B
        self.primes = primes_upto(B)
        P = 1
        for p in self.primes:
            P *= p
        self.P = gmpy2.mpz(P)

    def factor(self, z: int) -> Optional[Dict[int, int]]:
        """Exact: return {prime: exponent} if |z| is B-smooth (z != 0), else None."""
        z = abs(int(z))
        if z == 0:
            return None
        r = gmpy2.mpz(z)
        g = gmpy2.gcd(r, self.P)
        while g > 1:
            r //= g
            g = gmpy2.gcd(r, g)
        if r != 1:
            return None
        out = {}
        for p in self.primes:
            if z % p == 0:
                e = 0
                while z % p == 0:
                    z //= p
                    e += 1
                out[p] = e
                if z == 1:
                    break
        return out
