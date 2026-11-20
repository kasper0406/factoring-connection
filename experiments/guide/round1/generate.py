"""Enumerate and label exact NFS candidates. Run with the repository's PARI environment."""
import concurrent.futures
import hashlib
import json
import math
import os
import random
import time
from pathlib import Path

import gmpy2
import numpy as np
import sympy as sp
from cypari2 import Pari

from settings import BITS, BOUND as B, POWER, REFS, SEED

Y = 31
LOCAL = list(sp.primerange(Y + 1))
PRIMORIAL = gmpy2.primorial(B)
SPLITS = ["fit", "cal", "select", "test"]
MASK = (1 << 64) - 1


def polynomial(deg, bits, rng, pari):
    x = sp.Symbol("x")
    while True:
        primes = [int(gmpy2.next_prime(rng.randrange(3 << (bits // 2 - 2), (1 << (bits // 2)) - 1))) for _ in range(2)]
        p, q = primes
        n = p*q
        if p == q or n.bit_length() != bits or max(p, q).bit_length() != bits // 2:
            continue
        m = int(gmpy2.iroot(n, deg)[0])
        co, t = [], n
        while t:
            t, digit = divmod(t, m)
            co.append(digit)
        if len(co) != deg + 1 or not sp.Poly.from_list(co[::-1], gens=x, domain=sp.ZZ).is_irreducible:
            continue
        assert pari(p).isprime() and pari(q).isprime()
        assert sum(c*m**i for i, c in enumerate(co)) == n
        return p, q, m, co


def coordinates(m, split, rng):
    residue = (m + 15) % 31 - 15
    rows, seen = [], set()
    block_size = 1024 if POWER else REFS
    for _ in range(REFS // block_size):
        block = []
        for k in range(block_size):
            g = k % 2
            while True:
                b, i = rng.randrange(513, 769), rng.randrange(-224, 224)
                a = 31*i + residue*b if g else i
                if math.gcd(a, b) == 1 and (a, b) not in seen:
                    seen.add((a, b))
                    block.append((a, b, i, g))
                    break
        rng.shuffle(block)
        rows.extend(block)
    for g in range(2):
        for b in range(1 + 256*g, 257 + 256*g):
            if ((b-1)//8) % 4 != split:
                continue
            for i in range(-224, 224):
                a = 31*i + residue*b if g else i
                if math.gcd(a, b) == 1:
                    rows.append((a, b, i, g))
    return np.array(rows, dtype=np.int32)


def independent_label(c):
    t = gmpy2.mpz(c)
    while t != 1:
        g = gmpy2.gcd(t, PRIMORIAL)
        if g == 1:
            return False
        t //= g
    return True


def run(task):
    idx, deg, split, bits = task
    stamp = time.monotonic()
    cpu = time.process_time()
    rng = random.Random(SEED + 1009*idx)
    pari = Pari()
    p, q, m, co = polynomial(deg, bits, rng, pari)
    rows = coordinates(m, split, rng)
    # Nonzero for these irreducible polynomials and the enumerated rational ranges.
    val = np.zeros((len(rows), 2, len(LOCAL)), dtype=np.uint8)
    words = np.zeros((len(rows), 2, 2), dtype=np.uint64)
    labels = np.zeros((len(rows), 2), dtype=np.uint8)
    for k, (a, b, i, g) in enumerate(rows.tolist()):
        assert not g or (a-m*b) % 31 == 0
        norms = [abs(a-m*b), abs(sum(v*a**j*b**(deg-j) for j, v in enumerate(co)))]
        for side, n in enumerate(norms):
            assert n > 0
            c = n
            for j, prime in enumerate(LOCAL):
                if c % prime == 0:
                    e = 0
                    while c % prime == 0:
                        c //= prime
                        e += 1
                    val[k, side, j] = e
            assert c*math.prod(prime**int(e) for prime,e in zip(LOCAL, val[k, side]) if e) == n
            assert c.bit_length() <= 128
            words[k, side] = c & MASK, c >> 64
            factors = pari(c).factor() if POWER else pari(c).factor(B)
            smooth = all(int(t) <= B for t in factors[0])
            assert smooth == independent_label(c), (idx,k,side,c,smooth)
            labels[k,side] = smooth
    path = Path(f"data/n{idx:02}.npz")
    np.savez_compressed(path, coords=rows, valuation=val, cofactor_words=words, sides=labels)
    info = dict(id=idx, degree=deg, split=SPLITS[split], bits=bits, p=str(p), q=str(q), N=str(p*q), m=m,
                coefficients_ascending=co, reference_count=REFS, targets=len(rows)-REFS,
                label_audits=int(2*len(rows)), seconds=time.monotonic()-stamp, cpu_seconds=time.process_time()-cpu,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    Path(f"data/n{idx:02}.json").write_text(json.dumps(info, indent=2)+"\n")
    return info


if __name__ == "__main__":
    tasks = [(d,s,bits) for d in (3,4) for s, grid in enumerate(BITS) for bits in grid]
    tasks = [(j,*t) for j,t in enumerate(tasks)]
    start = time.monotonic()
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        infos = []
        for info in pool.map(run, tasks, chunksize=1):
            infos.append(info)
            print(info["id"], info["split"], info["targets"], round(info["seconds"],1), flush=True)
    assert len({v["N"] for v in infos}) == len(infos)
    assert sum(v["targets"] for v in infos) >= 10**6
    Path("data/manifest.json").write_text(json.dumps(infos,indent=2)+"\n")
    Path("logs/generation_time.json").write_text(json.dumps(dict(wall=time.monotonic()-start,
        worker_cpu=sum(x["cpu_seconds"] for x in infos))))
