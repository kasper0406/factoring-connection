"""Exact tests for the elementary seed lemmas; small instances only."""
from math import comb, gcd, isqrt


def least_divisor_reference(value: int) -> int:
    if value < 2:
        raise ValueError("value must be at least 2")
    for divisor in range(2, isqrt(value) + 1):
        if value % divisor == 0:
            return divisor
    return value


def check_seed_lemmas(limit: int = 200) -> int:
    if not 4 <= limit <= 2000:
        raise ValueError("use a small test limit between 4 and 2000")
    composites = 0
    for modulus in range(2, limit + 1):
        first = least_divisor_reference(modulus)
        # Only interior coefficients: the endpoints of F_N are removed.
        coeffs = [comb(modulus, k) % modulus
                  for k in range(1, modulus)]
        if first == modulus:
            assert not any(coeffs), (modulus, "prime identity")
            continue

        composites += 1
        assert not any(coeffs[:first - 1])
        assert coeffs[first - 1] == modulus // first
        for k, coefficient in enumerate(coeffs, start=1):
            assert k * coefficient % modulus == 0
            if coefficient:
                assert 1 < gcd(k, modulus) < modulus

        # product = (K-1)! mod modulus at each iteration.
        product = 1
        for cutoff in range(1, isqrt(modulus) + 2):
            prefix_zero = not any(coeffs[:cutoff - 1])
            assert prefix_zero == (cutoff <= first)
            assert prefix_zero == (gcd(product, modulus) == 1)
            product = product * cutoff % modulus

    for point in range(561):
        assert (pow(1 + point, 561, 561)
                - 1 - pow(point, 561, 561)) % 561 == 0
    return composites



# ---------------------------------------------------------------------------
# Extensions (infraA).  Everything above this line is BRIEF section 9 verbatim.
# Reference / harness code: uses trial division and explicit coefficients.
# ---------------------------------------------------------------------------
def factorize_reference(value: int) -> dict:
    """Prime factorisation by trial division (reference only)."""
    out, d = {}, 2
    while d * d <= value:
        while value % d == 0:
            out[d] = out.get(d, 0) + 1
            value //= d
        d += 1
    if value > 1:
        out[value] = out.get(value, 0) + 1
    return out


def digits(value: int, base: int) -> list:
    out = []
    while value:
        out.append(value % base)
        value //= base
    return out


def dominated(j: int, m: int, base: int) -> bool:
    """True iff every base-`base` digit of j is <= the matching digit of m."""
    while j or m:
        if j % base > m % base:
            return False
        j //= base
        m //= base
    return True


def carries(a: int, b: int, base: int) -> int:
    """Number of carries when adding a and b in base `base` (Kummer)."""
    c = count = 0
    while a or b or c:
        s = a % base + b % base + c
        c = 1 if s >= base else 0
        count += c
        a //= base
        b //= base
    return count


def kummer_support(modulus: int) -> set:
    """Interior k with C(N,k) != 0 mod N, via Kummer: exists p^e || N with
    carries_p(k, N-k) < e."""
    fac = factorize_reference(modulus)
    return {k for k in range(1, modulus)
            if any(carries(k, modulus - k, p) < e for p, e in fac.items())}


def lucas_support_squarefree_part(modulus: int, p: int) -> set:
    """For p || N: {k : C(N,k) != 0 mod p} = {p*j : j digit-dominated by N/p in base p}."""
    m = modulus // p
    return {p * j for j in range(1, m) if dominated(j, m, p)}


def check_extended(limit: int = 1000) -> dict:
    """Seed lemmas + Lucas/Kummer support description for all composites <= limit."""
    stats = {"composites": 0, "semiprimes_pq": 0, "prime_powers": 0,
             "three_factor": 0, "scalar_zero_composites": []}
    for modulus in range(4, limit + 1):
        fac = factorize_reference(modulus)
        if len(fac) == 1 and sum(fac.values()) == 1:
            continue
        stats["composites"] += 1
        first = min(fac)
        coeffs = [comb(modulus, k) % modulus for k in range(1, modulus)]
        support = {k for k, c in enumerate(coeffs, start=1) if c}
        # 2.1 first surviving coefficient
        assert min(support) == first and coeffs[first - 1] == modulus // first
        # 2.2 k a_k = 0 and gcd test
        for k in support:
            assert k * coeffs[k - 1] % modulus == 0 and 1 < gcd(k, modulus) < modulus
        # 2.4 prefix <=> factorial unit, all K up to N
        product = 1
        for cutoff in range(1, modulus + 1):
            prefix_zero = all(k >= cutoff for k in support)
            assert prefix_zero == (cutoff <= first) == (gcd(product, modulus) == 1)
            product = product * cutoff % modulus
        # 2.5 F_N' = 0: coefficient k*a_k of X^(k-1), plus N*(1)... all zero mod N
        assert all(k * c % modulus == 0 for k, c in enumerate(coeffs, start=1))
        # Kummer description of the support (general N)
        assert support == kummer_support(modulus), modulus
        # Lucas description for each prime appearing to the first power
        for p, e in fac.items():
            if e == 1:
                sp = {k for k in support if comb(modulus, k) % p}
                assert sp == lucas_support_squarefree_part(modulus, p), (modulus, p)
        if len(fac) == 2 and all(e == 1 for e in fac.values()):
            stats["semiprimes_pq"] += 1
            p, q = sorted(fac)
            # mod-q part: {q*j : 1 <= j < p}; mod-p part: {p*j : j <=_p q}; disjoint
            assert lucas_support_squarefree_part(modulus, q) == {q * j for j in range(1, p)}
            size = (p - 1)
            prod_digits = 1
            for d in digits(q, p):
                prod_digits *= d + 1
            size += prod_digits - 2
            assert len(support) == size, (modulus, len(support), size)
            assert len(support) >= p - 1  # 2.6 support lower bound
        if len(fac) == 1:
            stats["prime_powers"] += 1
            p = first
            assert support == set(range(p, modulus, p)), modulus
        if len(fac) == 3 and all(e == 1 for e in fac.values()):
            stats["three_factor"] += 1
            union = set()
            for p in fac:
                union |= lucas_support_squarefree_part(modulus, p)
            assert support == union
        # 2.6 scalar-zero trap: record composites with F_N(a) = 0 for all a
        if all((pow(1 + a, modulus, modulus) - 1 - pow(a, modulus, modulus)) % modulus == 0
               for a in range(modulus)):
            stats["scalar_zero_composites"].append(modulus)
    return stats


if __name__ == "__main__":
    count = check_seed_lemmas()
    print(f"Seed checks passed for {count} small composites.")
    print("Extended checks:", check_extended())
