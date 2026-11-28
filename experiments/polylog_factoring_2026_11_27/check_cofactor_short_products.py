"""Cross-check signed low-product lists, then split planted products without their factors."""

from math import isqrt
from random import Random

from sympy import isprime, nextprime

if __package__:
    from .cofactor_short_products import (centered, factor_from_near_products, factor_from_short_products,
                                         near_product_classes, near_word_modulus, short_product_classes, word_modulus)
else:
    from cofactor_short_products import (centered, factor_from_near_products, factor_from_short_products,
                                        near_product_classes, near_word_modulus, short_product_classes, word_modulus)


def check_classes():
    checked = 0
    for k in range(2, 10):
        m = 1 << k
        for n in range(1, m, 2):
            for budget in (0, 1, 3):
                head = centered(n, m)
                expected = {a for a in range(1, m, 2)
                            if abs(centered(a, m)*centered(n*pow(a, -1, m), m)-head) <= budget*m}
                assert set(short_product_classes(n, m, budget)) == expected
                near = set()
                for a in expected:
                    x = abs(centered(a, m))
                    y = abs(centered(n*pow(a, -1, m), m))
                    if (x+y)//2-(isqrt(x*y-1)+1) < 2:
                        near.add(a)
                assert set(near_product_classes(n, m, budget, 2)) == near
                checked += 1
    return checked


def planted(bits, rng):
    low = 1 << (bits-1)
    m = word_modulus(low*low)
    while True:
        half = isqrt(m)//2
        x, y = (rng.randrange(half, 2*half) | 1 for _ in range(2))
        if min(x, y) < 3 or (x-y) % 16 == 0 or (x+y) % 16 == 0:
            continue
        start_p = low+rng.randrange(low//8)
        start_q = 3*low//2+rng.randrange(low//8)
        p = start_p+(x-start_p) % m
        q = start_q+(-y-start_q) % m
        while not isprime(p):
            p += m
        while not isprime(q):
            q += m
        n = p*q
        if word_modulus(n) != m or p >= q or q > 2*p:
            continue
        carry = (centered(p, m)*centered(q, m)-centered(n, m))//m
        return p, q, carry


def check_cases():
    rng = Random(581915)
    data = []
    for bits in (20, 28, 36):
        p, q, carry = planted(bits, rng)
        n = p*q
        stats = {}
        found = factor_from_short_products(n, 1, stats=stats)
        assert abs(carry) <= 1 and found in (p, q)
        data.append((bits, n, found, stats["modulus"], carry, stats["classes"], stats["method"]))
    controls = []
    for bits in (16, 32, 48):
        hits = 0
        for _ in range(100):
            p, q = (int(nextprime((1 << (bits-1))+rng.getrandbits(bits-2))) for _ in range(2))
            n = p*q
            m = word_modulus(n)
            promise = abs(centered(p, m)*centered(q, m)-centered(n, m)) <= m
            assert (p % m in short_product_classes(n, m, 1)) == promise
            hits += promise
        controls.append((bits, hits))
    return data, controls


def check_near():
    rows = []
    rng = Random(710325)
    for bits in (48, 64, 96, 128):
        low = 1 << (bits-1)
        m = near_word_modulus(low*low)
        x = (5*isqrt(m)//8) | 1
        delta = 8*(isqrt(x)//16)+2
        if (2*x+delta) % 16 == 0:
            delta += 8
        start_p = low+rng.randrange(low//8)
        start_q = 3*low//2+rng.randrange(low//8)
        p = start_p+(x-start_p) % m
        q = start_q+(x+delta-start_q) % m
        while not isprime(p):
            p += m
        while not isprime(q):
            q += m
        n = p*q
        assert m == near_word_modulus(n) and p < q < 2*p and p.bit_length() == q.bit_length() == bits
        stats = {}
        found = factor_from_near_products(n, 1, 1, stats=stats)
        assert found in (p, q) and (p-q) % 16 and (p+q) % 16
        a, c = (p+q)//2, (q-p)//2
        direct = min(1+(a-isqrt(n-1)-1)//(c & -c)**2, 1+c//(a & -a)**2)
        assert direct > 1
        rows.append((bits, n, found, delta, stats["classes"], stats["affine_windows"], direct))
    return rows


if __name__ == "__main__":
    print("Exact low-product lists:", check_classes())
    data, controls = check_cases()
    print("Planted input-only splits (prime bits, N, divisor, m, actual carry, classes, method):")
    for row in data:
        print(row)
    print("True factor residue in list among 100 unplanted inputs per size:", controls)
    print("Nested Fermat splits (prime bits, N, divisor, suffix gap, classes, windows, direct signed budget):")
    for row in check_near():
        print(row)
