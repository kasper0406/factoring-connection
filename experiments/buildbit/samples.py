from pathlib import Path
from random import Random
from gmpy2 import is_prime

SEED = 20261126

def prime(rng, bits, residue, modulus):
    low, high = 1 << (bits-1), 1 << bits
    start = low + (residue-low) % modulus
    count = (high-1-start)//modulus + 1
    while True:
        v = start + modulus*rng.randrange(count)
        if is_prime(v):
            return v

def collect(topic, sizes, count):
    tag = {'cf': 1, 'cubic': 2, 'carry': 3}[topic]
    for n in sizes:
        rng = Random(SEED*10000 + n*10 + tag)
        rows = []
        for label in (0, 1):
            for i in range(count):
                if topic == 'cf':
                    a = b = 1 if label else 3
                    m = h = 4
                elif topic == 'cubic':
                    a = rng.choice([2,3,4,5] if label else [1,6])
                    b = pow(a, -1, 7)
                    m = h = 7
                else:
                    a, b = (9 if label else 1), 3
                    m, h = 16, 8
                while True:
                    p = prime(rng, n//2, a, m)
                    q = prime(rng, n//2, b, h)
                    if p != q and (p*q).bit_length() == n:
                        break
                rows.append((p*q, p, q, label))
        rng.shuffle(rows)
        Path(f'{topic}{n}.tsv').write_text(''.join('\t'.join(map(str,row))+'\n' for row in rows))

if __name__ == '__main__':
    collect('cf', (40,52,64,80), 2048)
    collect('cubic', (40,60,80), 512)
    collect('carry', (40,60,80), 256)
