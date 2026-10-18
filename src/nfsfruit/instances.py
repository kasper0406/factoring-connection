"""HARNESS: seeded composite generators with residue-class constraints mod 4.

Seeding follows src/common/instances.py: random.Random(f"{family}|{params}|{seed}").
Classes:  blum  : p = q = 3 (mod 4)
          mixed : p = 1, q = 3 (mod 4) (roles assigned at random)
          both1 : p = q = 1 (mod 4)
          tri   : three distinct primes of ~bits/3 bits, no congruence constraint
"""
from __future__ import annotations
import random
import gmpy2

def _rng(cls, bits, seed):
    return random.Random(f"nfsfruit|{cls}|bits={bits}|{seed}")

def _prime(rng, b, res=None):
    while True:
        x = rng.getrandbits(b) | (1 << (b - 1)) | (1 << (b - 2)) | 1
        if res is not None and x % 4 != res:
            continue
        if gmpy2.is_prime(x):
            return int(x)

def composite(cls: str, bits: int, seed: int):
    """Return (N, sorted factor list). HARNESS ONLY."""
    rng = _rng(cls, bits, seed)
    while True:
        if cls == "tri":
            b = [bits // 3, bits // 3, bits - 2 * (bits // 3)]
            fs = [_prime(rng, x) for x in b]
        else:
            hb = [(bits + 1) // 2, bits // 2]
            if cls == "blum":
                res = [3, 3]
            elif cls == "both1":
                res = [1, 1]
            elif cls == "mixed":
                res = [1, 3] if rng.random() < 0.5 else [3, 1]
            else:
                raise ValueError(cls)
            fs = [_prime(rng, hb[0], res[0]), _prime(rng, hb[1], res[1])]
        if len(set(fs)) != len(fs):
            continue
        N = 1
        for f in fs:
            N *= f
        if N.bit_length() == bits:
            return N, sorted(fs)
