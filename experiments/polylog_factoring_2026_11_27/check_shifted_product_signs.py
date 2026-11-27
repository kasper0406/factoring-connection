"""Decode a sum residue from supplied, possibly corrupted, shifted-product signs."""

from math import isqrt

from sympy import nextprime


def sign_table(prime):
    signs = [-1] * prime
    for x in range(prime):
        signs[x*x % prime] = 1
    return signs  # Convention: zero gets +1.


def decode_residue(n, prime, observations):
    if len(observations) != prime or any(v not in (-1, 1) for v in observations):
        raise ValueError("Expected one binary sign per shift")
    signs = sign_table(prime)
    tolerance = prime // 8
    for s in range(prime):
        if (s*s-4*n) % prime == 0:
            continue
        errors = sum(signs[(c*c-s*c+n) % prime] != observed
                     for c, observed in enumerate(observations))
        if errors <= tolerance:
            return s
    return None


def recover(n, observe):
    prime, modulus, total = 97, 1, 0
    tried = 1
    while modulus <= n and tried <= n*n:
        prime = int(nextprime(prime))
        tried *= prime
        if n % prime == 0:
            return prime, n // prime
        s = decode_residue(n, prime, observe(prime))
        if s is not None:
            total += modulus * ((s-total) * pow(modulus, -1, prime) % prime)
            modulus *= prime
    if modulus <= n:
        raise ValueError("Insufficient consistent residues")
    delta = total*total-4*n
    if delta < 0 or isqrt(delta)**2 != delta:
        raise ValueError("Decoded residues do not give a factor sum")
    p = (total-isqrt(delta)) // 2
    if p <= 1 or n % p:
        raise ValueError("Invalid factor candidate")
    return p, n // p


if __name__ == "__main__":
    for p, q in [(1009, 1000003), (1009, 2221), (42697, 70019)]:
        def planted(prime):
            signs = sign_table(prime)
            result = [signs[(p-c)*(q-c) % prime] for c in range(prime)]
            for i in range(prime // 8):
                result[(17*i+3) % prime] *= -1
            return result

        assert sorted(recover(p*q, planted)) == sorted((p, q))
        print(p*q, p, q)
