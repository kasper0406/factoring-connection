"""Seeded instance generators (shared infrastructure, owner: infraA).  HARNESS code.

Every generator returns an `Instance`.  ``inst.N`` (and ``inst.public()``) is
what a candidate algorithm may see.  ``inst.labels`` holds the factorisation and
other truth; it is for the harness only and is hidden from ``repr``.

Sampling procedure (applies to every family)
--------------------------------------------
* Randomness: ``random.Random(f"{family}|{params}|{seed}")`` -- Python's
  string seeding hashes the string with SHA-512, so instances are a
  deterministic function of (family, params, seed) and are independent of the
  factors (no factor is encoded in a seed or filename).
* ``random_prime(rng, b)``: draw a uniform b-bit odd integer with the top bit
  set, test with ``gmpy2.is_prime`` (BPSW + Miller-Rabin), redraw until prime.
  This is the standard "uniform b-bit prime" sampler (density-biased toward
  primes after large gaps, as usual).
* Balanced semiprime of n bits: p, q independent ``random_prime`` of
  ceil(n/2), floor(n/2) bits with the top TWO bits set (so N has exactly n
  bits), resampled until p != q.  Labels: p < q sorted.
* Unbalanced semiprime: p of `small_bits` bits, q of n - small_bits bits (top
  two bits set on q, resample until N has n bits).
* Prime power p^k of exactly n bits, p uniform prime in the admissible window.
* Products of `r` distinct primes (r = 3, 4) of ~n/r bits each, resampled until
  distinct and N has exactly n bits.
* Carmichael: the fixed list ``CARMICHAEL_SMALL`` (all Carmichael numbers below
  600000 we list, verified by Korselt's criterion in tests) or Chernick numbers
  (6k+1)(12k+1)(18k+1) with all three factors prime, k drawn uniformly from a
  range and incremented until valid.
* Primes: ``random_prime`` of n bits.
* Diagnostic families (labelled ``diag_*``, NOT the intended hard distribution):
  - ``diag_near_equal``: p random prime of n/2 bits, q = next prime after p.
  - ``diag_smooth_pm1``: p = 2*prod(random small primes <= B)+1 prime,
    q a random prime; so p-1 is B-smooth.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field
from math import prod
from typing import Any, Dict, List, Optional

import gmpy2


@dataclass
class Instance:
    N: int
    family: str
    params: Dict[str, Any]
    seed: int
    labels: Dict[str, Any] = field(repr=False, default_factory=dict)

    def public(self) -> Dict[str, Any]:
        """What a candidate may see: N plus the public sampling configuration."""
        return {"N": self.N, "bits": self.N.bit_length(), "family": self.family,
                "params": dict(self.params), "seed": self.seed}

    @property
    def factors(self) -> List[int]:
        """Harness only: sorted prime factors with multiplicity."""
        return list(self.labels["factors"])


def _rng(family: str, params: Dict[str, Any], seed: int) -> random.Random:
    key = ",".join(f"{k}={params[k]}" for k in sorted(params))
    return random.Random(f"{family}|{key}|{seed}")


def is_prime(x: int) -> bool:
    """Primality (harness / instance generation only)."""
    return bool(gmpy2.is_prime(x, 50))


def random_prime(rng: random.Random, bits: int, top2: bool = False) -> int:
    if bits < 2:
        raise ValueError("bits >= 2")
    if bits == 2:
        return rng.choice([2, 3]) if not top2 else 3
    while True:
        x = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if top2:
            x |= 1 << (bits - 2)
        if is_prime(x):
            return x


def _mk(N: int, family: str, params: Dict[str, Any], seed: int, factors: List[int], **extra) -> Instance:
    factors = sorted(factors)
    assert prod(factors) == N
    labels = {"factors": factors, "p": factors[0], "distinct": sorted(set(factors))}
    labels.update(extra)
    return Instance(N=N, family=family, params=params, seed=seed, labels=labels)


def balanced_semiprime(bits: int, seed: int) -> Instance:
    params = {"bits": bits}
    rng = _rng("balanced", params, seed)
    while True:
        p = random_prime(rng, (bits + 1) // 2, top2=True)
        q = random_prime(rng, bits // 2, top2=True)
        N = p * q
        if p != q and N.bit_length() == bits:
            p, q = sorted((p, q))
            return _mk(N, "balanced", params, seed, [p, q], q=q)


BALANCED_BITS = (16, 24, 32, 48, 64, 96, 128)


def unbalanced_semiprime(bits: int, small_bits: int, seed: int) -> Instance:
    params = {"bits": bits, "small_bits": small_bits}
    rng = _rng("unbalanced", params, seed)
    while True:
        p = random_prime(rng, small_bits)
        q = random_prime(rng, bits - small_bits, top2=True)
        N = p * q
        if p != q and N.bit_length() == bits:
            return _mk(N, "unbalanced", params, seed, [p, q])


def prime_power(bits: int, k: int, seed: int) -> Instance:
    """p^k with exactly `bits` bits: p uniform prime in [ceil(2^((bits-1)/k)), floor((2^bits-1)^(1/k))]."""
    params = {"bits": bits, "k": k}
    rng = _rng("prime_power", params, seed)
    lo = int(gmpy2.iroot(gmpy2.mpz(1) << (bits - 1), k)[0])
    while lo ** k < (1 << (bits - 1)):
        lo += 1
    hi = int(gmpy2.iroot(gmpy2.mpz((1 << bits) - 1), k)[0])
    if lo > hi:
        raise ValueError("no prime power of that shape")
    for _ in range(100000):
        p = rng.randint(lo, hi)
        if is_prime(p):
            return _mk(p ** k, "prime_power", params, seed, [p] * k)
    raise ValueError("no prime found in window")


def multi_prime(bits: int, r: int, seed: int) -> Instance:
    """Product of r distinct primes, N of exactly `bits` bits."""
    params = {"bits": bits, "r": r}
    rng = _rng("multi_prime", params, seed)
    sizes = [bits // r + (1 if i < bits % r else 0) for i in range(r)]
    while True:
        ps = [random_prime(rng, max(2, s), top2=(s >= 3)) for s in sizes]
        N = prod(ps)
        if len(set(ps)) == r and N.bit_length() == bits:
            return _mk(N, "multi_prime", params, seed, ps)


CARMICHAEL_SMALL = [561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841, 29341,
                    41041, 46657, 52633, 62745, 63973, 75361, 101101, 115921,
                    126217, 162401, 172081, 188461, 252601, 278545, 294409,
                    314821, 334153, 340561, 399001, 410041, 449065, 488881,
                    512461]


def _small_factor(n: int) -> List[int]:
    out, d = [], 2
    while d * d <= n:
        while n % d == 0:
            out.append(d)
            n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def carmichael_small(index: int) -> Instance:
    N = CARMICHAEL_SMALL[index]
    return _mk(N, "carmichael", {"index": index}, index, _small_factor(N))


def carmichael_chernick(kbits: int, seed: int) -> Instance:
    """(6k+1)(12k+1)(18k+1) with all three prime; k starts at a random kbits-bit value."""
    params = {"kbits": kbits}
    rng = _rng("chernick", params, seed)
    k = rng.getrandbits(kbits) | (1 << (kbits - 1))
    while True:
        ps = [6 * k + 1, 12 * k + 1, 18 * k + 1]
        if all(is_prime(x) for x in ps):
            return _mk(prod(ps), "carmichael", params, seed, ps, k=k)
        k += 1


def prime_instance(bits: int, seed: int) -> Instance:
    params = {"bits": bits}
    rng = _rng("prime", params, seed)
    p = random_prime(rng, bits)
    return Instance(N=p, family="prime", params=params, seed=seed,
                    labels={"factors": [p], "p": p, "distinct": [p]})


def diag_near_equal(bits: int, seed: int) -> Instance:
    params = {"bits": bits}
    rng = _rng("diag_near_equal", params, seed)
    while True:
        p = random_prime(rng, bits // 2, top2=True)
        q = int(gmpy2.next_prime(p))
        if (p * q).bit_length() == bits:
            return _mk(p * q, "diag_near_equal", params, seed, [p, q])


_SMALL_PRIMES = [x for x in range(3, 2000) if all(x % d for d in range(2, int(x ** 0.5) + 1))]


def diag_smooth_pm1(bits: int, seed: int, B: int = 50) -> Instance:
    params = {"bits": bits, "B": B}
    rng = _rng("diag_smooth_pm1", params, seed)
    small = [x for x in _SMALL_PRIMES if x <= B]
    hb = bits // 2
    while True:
        m = 2
        while m.bit_length() < hb - 1:
            m *= rng.choice(small)
        p = m + 1
        if p.bit_length() > hb or not is_prime(p):
            continue
        q = random_prime(rng, bits - p.bit_length() + 1)
        N = p * q
        if p != q:
            return _mk(N, "diag_smooth_pm1", params, seed, [p, q], smooth_side=p)


FAMILIES = {
    "balanced": balanced_semiprime,
    "unbalanced": unbalanced_semiprime,
    "prime_power": prime_power,
    "multi_prime": multi_prime,
    "carmichael_small": carmichael_small,
    "carmichael_chernick": carmichael_chernick,
    "prime": prime_instance,
    "diag_near_equal": diag_near_equal,
    "diag_smooth_pm1": diag_smooth_pm1,
}


def make(family: str, *args, **kwargs) -> Instance:
    return FAMILIES[family](*args, **kwargs)


def batch(family: str, count: int, *args, seed0: int = 0, **kwargs) -> List[Instance]:
    return [FAMILIES[family](*args, seed=seed0 + i, **kwargs) for i in range(count)]
