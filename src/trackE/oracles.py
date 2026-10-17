"""DIAGNOSTIC ORACLES (HARNESS ONLY: these know the factors).

Each oracle takes an `Instance` (from src.common.instances) and returns ONLY the
partial-information object of an exact format; the decoders in coppersmith.py /
formats.py receive N and that object, never the instance.  Operational code
must never import this module.

Convention: "p" is the smaller prime factor (inst.labels['p']), "q" the other.
"""
from __future__ import annotations

import random
from typing import Dict, List, Optional, Sequence

from src.common.instances import Instance


def _pq(inst: Instance):
    p = inst.labels["p"]
    return p, inst.N // p


# --- (a) residue / low bits -------------------------------------------------

def oracle_low_bits(inst: Instance, k: int, flip_positions: Sequence[int] = ()) -> Dict:
    """{k, r}: r = p mod 2^k, with the listed bit positions flipped (noise model)."""
    p, _ = _pq(inst)
    r = p % (1 << k)
    for pos in flip_positions:
        r ^= 1 << pos
    return {"k": k, "r": r}


def oracle_residue(inst: Instance, M: int) -> Dict:
    p, _ = _pq(inst)
    return {"M": M, "r": p % M}


# --- (b) high bits / interval ----------------------------------------------

def oracle_high_bits(inst: Instance, h: int, noise: int = 0) -> Dict:
    """Top h bits of p as an interval [A, A+W]; `noise` = additive error on A
    (the oracle's A is off by `noise`, and it does NOT tell the decoder)."""
    p, _ = _pq(inst)
    b = p.bit_length()
    sh = b - h
    A = (p >> sh) << sh
    return {"A": A + noise, "W": (1 << sh) - 1, "h": h, "bits_p": b}


# --- derived formats ----------------------------------------------------------

def oracle_sum_mod(inst: Instance, M: int) -> Dict:
    """s = p + q mod M (equivalently phi(N) = N + 1 - s mod M)."""
    p, q = _pq(inst)
    return {"M": M, "s": (p + q) % M}


def oracle_residues(inst: Instance, primes: Sequence[int], wrong: Sequence[int] = (),
                    rng: Optional[random.Random] = None) -> Dict:
    """p mod l for each l (labelled consistently: always the smaller factor).
    Indices in `wrong` get a uniformly random *incorrect* residue."""
    p, _ = _pq(inst)
    rng = rng or random.Random(0)
    res = []
    for i, l in enumerate(primes):
        r = p % l
        if i in wrong:
            r = rng.choice([x for x in range(l) if x != r])
        res.append(r)
    return {"primes": list(primes), "residues": res}


def oracle_legendre(inst: Instance, primes: Sequence[int]) -> Dict:
    """(p / l) in {+1,-1} for odd primes l not dividing p."""
    p, _ = _pq(inst)
    return {"primes": list(primes), "legendre": [1 if pow(p, (l - 1) // 2, l) == 1 else -1 for l in primes]}


def oracle_power_residue_class(inst: Instance, l: int, e: int) -> Dict:
    """Index of p mod l in (F_l^*)/(F_l^*)^e, e | l-1, via a fixed generator: the
    value  (p^{(l-1)/e} mod l)  (an e-th root of unity) -- log2(e) bits."""
    p, _ = _pq(inst)
    assert (l - 1) % e == 0
    return {"l": l, "e": e, "chi": pow(p, (l - 1) // e, l)}


def oracle_N_mod_pm1(inst: Instance) -> Dict:
    """e_p = N mod (p-1) = q mod (p-1): the linear complexity found by Track B/D's D3."""
    p, _ = _pq(inst)
    return {"e": inst.N % (p - 1)}


def oracle_random_bits(inst: Instance, delta: float, nbits: int, rng: random.Random) -> Dict:
    """Each of the low `nbits` bits of p and of q revealed independently with prob delta."""
    p, q = _pq(inst)
    kp = {i: (p >> i) & 1 for i in range(nbits) if rng.random() < delta}
    kq = {i: (q >> i) & 1 for i in range(nbits) if rng.random() < delta}
    return {"nbits": nbits, "known_p": kp, "known_q": kq}
