import math
import random

def compute_E_exact(max_a: int):
    """Exact integer jet sequence E_a = D^a s(0) for s'' = -2 s^3, s(0)=0, s'(0)=1."""
    # Nonzero E_a only at a = 1 mod 4; S2_k = [u^k/k!](s^2) nonzero only at k = 2 mod 4.
    E = {1: 1}
    S2 = {}
    # Precompute binomials
    C = [[0] * (max_a + 1) for _ in range(max_a + 1)]
    for n in range(max_a + 1):
        C[n][0] = 1
        for k in range(1, n + 1):
            C[n][k] = C[n - 1][k - 1] + C[n - 1][k]

    for a in range(5, max_a + 1, 4):
        k_s2 = a - 3  # k_s2 == 2 mod 4
        # Compute S2[k_s2] = sum_{i+j = k_s2} C[k_s2][i] E[i] E[j]
        val_s2 = 0
        for i in range(1, k_s2, 4):
            j = k_s2 - i
            val_s2 += C[k_s2][i] * E[i] * E[j]
        S2[k_s2] = val_s2

        n = a - 2  # E_{n+2} = -2 sum_{k+m = n} C[n][k] S2[k] E[m]
        val_e = 0
        for k in range(2, n, 4):
            m = n - k
            val_e += C[n][k] * S2[k] * E[m]
        E[a] = -2 * val_e
    return E

def compute_E_mod_N(target_a: int, N: int):
    """Pure mod-N computation of E_a mod N (zero divisions, zero rational arithmetic)."""
    E = {1: 1 % N}
    S2 = {}
    # Rolling / row binomials up to target_a - 2 mod N
    max_n = target_a - 2
    C = [[0] * (max_n + 1) for _ in range(max_n + 1)]
    for n in range(max_n + 1):
        C[n][0] = 1
        for k in range(1, n + 1):
            C[n][k] = (C[n - 1][k - 1] + C[n - 1][k]) % N

    for a in range(5, target_a + 1, 4):
        k_s2 = a - 3
        val_s2 = 0
        for i in range(1, k_s2, 4):
            j = k_s2 - i
            val_s2 = (val_s2 + C[k_s2][i] * E[i] * E[j]) % N
        S2[k_s2] = val_s2

        n = a - 2
        val_e = 0
        for k in range(2, n, 4):
            m = n - k
            val_e = (val_e + C[n][k] * S2[k] * E[m]) % N
        E[a] = (-2 * val_e) % N
    return E

def vp(n: int, p: int) -> int:
    if n == 0:
        return 10**9
    n = abs(n)
    c = 0
    while n % p == 0:
        n //= p
        c += 1
    return c

def vp_fact(a: int, p: int) -> int:
    c = 0
    while a >= p:
        a //= p
        c += a
    return c

def odd_power_sum(a: int, p: int) -> int:
    s = 0
    pk = p
    while pk <= a:
        s += a // pk
        pk *= p * p
    return s

def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True

def main():
    MAX_A = 401
    E = compute_E_exact(MAX_A)
    a_list = list(range(1, MAX_A + 1, 4))

    # --- P1: Inert primes ---
    inert_primes = [7, 11, 19, 23, 31]
    p1a_total, p1a_pass = 0, 0
    p1b_total, p1b_pass = 0, 0
    for p in inert_primes:
        for a in a_list:
            v = vp(E[a], p)
            if a >= p:
                p1a_total += 1
                if v >= 1:
                    p1a_pass += 1
            lb = odd_power_sum(a, p) - (int(math.log(a, p)) if a >= 1 else 0) - 1
            p1b_total += 1
            if v >= lb:
                p1b_pass += 1

    v7_401 = vp(E[401], 7)
    v7_fact_401 = vp_fact(401, 7)
    odd7_401 = odd_power_sum(401, 7)

    print(f"P1a (p | E_a for all a>=p, p inert): {p1a_pass}/{p1a_total}")
    print(f"P1b (v_p(E_a) >= sum_{{j>=0}} floor(a/p^(2j+1)) - floor(log_p a) - 1): {p1b_pass}/{p1b_total}")
    print(f"At a=401, p=7: v_7(E_401)={v7_401}, odd_sum={odd7_401}, v_7(401!)={v7_fact_401} "
          f"(theory: 401*7/48 = {401*7/48:.2f})")

    # --- P2: Split primes ---
    split_primes = [5, 13, 17, 29, 37]
    p2_total, p2_pass = 0, 0
    p2_unit_total, p2_unit_pass = 0, 0
    p2_aq_match = 0
    for q in split_primes:
        Aq = E[q] % q
        m = (q - 1) // 4
        Aq_gauss = ((-1) ** m * math.comb(2 * m, m)) % q
        if Aq == Aq_gauss and Aq != 0:
            p2_aq_match += 1
        for r in range(1, q, 4):
            k = 0
            while r + k * (q - 1) <= MAX_A:
                a = r + k * (q - 1)
                p2_total += 1
                if E[a] % q == (pow(Aq, k, q) * (E[r] % q)) % q:
                    p2_pass += 1
                if r == 1:
                    p2_unit_total += 1
                    if E[a] % q != 0:
                        p2_unit_pass += 1
                k += 1

    print(f"P2 Gauss Aq match & nonzero: {p2_aq_match}/{len(split_primes)}")
    print(f"P2 Kummer periodicity E_{{r+k(q-1)}} == Aq^k E_r (mod q): {p2_pass}/{p2_total}")
    print(f"P2 q nmid E_{{1+k(q-1)}}: {p2_unit_pass}/{p2_unit_total}")

    # --- P3: Archimedean Bost-slope ledger vs single-prime ratio ---
    varpi = math.gamma(0.25) ** 2 / (4.0 * math.sqrt(2.0 * math.pi))
    log_E401 = math.log(abs(E[401]))
    phF_401 = 401 * math.log(401) - 401
    emp_neg_log_om = (log_E401 - phF_401) / 401.0
    theo_neg_log_om = -math.log(varpi)

    max_copp_ratio = 0.0
    for p in inert_primes:
        for a in a_list:
            if a > 1:
                ratio = (vp(E[a], p) * math.log(p)) / math.log(abs(E[a]))
                if ratio > max_copp_ratio:
                    max_copp_ratio = ratio

    print(f"P3a archimedean slope at a=401: empirical={emp_neg_log_om:.6f}, "
          f"-log(varpi)={theo_neg_log_om:.6f}, diff={abs(emp_neg_log_om - theo_neg_log_om):.6f}")
    print(f"P3b max single-prime Coppersmith ratio v_p(E_a)log(p)/log|E_a|: {max_copp_ratio:.6f} (< 0.15)")

    # --- P4: Mod-N hidden-factor test on 30 balanced semiprimes (14-20 bits, seed 42) ---
    rng = random.Random(42)
    semiprimes = []
    seen = set()
    while len(semiprimes) < 30:
        bits = rng.choice([14, 16, 18, 20])
        lo = int(2 ** ((bits - 1) / 2))
        hi = int(2 ** (bits / 2))
        cands_p = [x for x in range(max(7, lo), hi + 1) if is_prime(x) and x % 4 == 3]
        cands_q = [x for x in range(max(5, lo), hi + 1) if is_prime(x) and x % 4 == 1]
        p = rng.choice(cands_p)
        q = rng.choice(cands_q)
        N = p * q
        if N not in seen:
            seen.add(N)
            semiprimes.append((N, p, q))

    low_hits = 0
    high_p_div = 0
    high_exact_splits = 0
    for N, p, q in semiprimes:
        a_low = 4 * int((N ** 0.25) / 4.0) + 1
        a_high = 4 * math.ceil(math.sqrt(2.0 * N) / 4.0) + 1
        assert a_low < p and a_high > p
        E_mod = compute_E_mod_N(a_high, N)
        g_low = math.gcd(E_mod[a_low], N)
        g_high = math.gcd(E_mod[a_high], N)
        if g_low == p:
            low_hits += 1
        if g_high % p == 0:
            high_p_div += 1
        if g_high == p:
            high_exact_splits += 1

    print(f"P4a low scale a_low ~ N^(1/4) exact factor p hits: {low_hits}/30 (criterion <= 2/30)")
    print(f"P4b high scale a_high > sqrt(N): p | E_{{a_high}} in {high_p_div}/30, "
          f"exact gcd=p in {high_exact_splits}/30 (criterion >= 25/30)")

if __name__ == "__main__":
    main()
