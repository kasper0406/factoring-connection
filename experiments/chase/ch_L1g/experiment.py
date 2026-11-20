"""Pre-registered exact-arithmetic experiments for Chase Lemma 1 (H / branch selection).
Run with: python3 /Volumes/git/factoring-connection/codex/scratch/ch_L1g/experiment.py
"""
import json
import math
import random
from collections import Counter
from pathlib import Path
from time import perf_counter

OUT_DIR = Path(__file__).resolve().parent
W_BUDGET = 4000


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def next_prime(n: int) -> int:
    c = n + 1 if n % 2 == 0 else n + 2
    if n <= 1:
        return 2
    if n == 2:
        return 3
    c = n if n % 2 == 1 else n + 1
    while not is_prime(c):
        c += 2
    return c


def primes_upto(B: int) -> list[int]:
    if B < 2:
        return []
    sieve = bytearray(b"\x01") * (B + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, math.isqrt(B) + 1):
        if sieve[i]:
            sieve[i * i : B + 1 : i] = b"\x00" * (((B - i * i) // i) + 1)
    return [i for i, v in enumerate(sieve) if v]


def tonelli_shanks(a: int, p: int) -> list[int]:
    a %= p
    if a == 0:
        return [0]
    if p == 2:
        return [a]
    if pow(a, (p - 1) // 2, p) != 1:
        return []
    if p % 4 == 3:
        r = pow(a, (p + 1) // 4, p)
        return sorted({r, p - r})
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m = s
    c = pow(z, q, p)
    t = pow(a, q, p)
    r = pow(a, (q + 1) // 2, p)
    while t != 1:
        i = 1
        t2 = (t * t) % p
        while t2 != 1:
            t2 = (t2 * t2) % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m = i
        c = (b * b) % p
        t = (t * c) % p
        r = (r * b) % p
    assert (r * r - a) % p == 0
    return sorted({r, p - r})


def reduce_form(a: int, b: int, c: int, N: int) -> tuple[int, int, int]:
    assert b * b - 4 * a * c == -4 * N and a > 0
    while True:
        if a > c or (a == c and b < 0):
            a, b, c = c, -b, a
        if abs(b) > a or (b == -a):
            q = (b + a) // (2 * a) if b >= 0 else -((-b + a - 1) // (2 * a))
            b -= 2 * q * a
            c = (b * b + 4 * N) // (4 * a)
        else:
            break
    assert b * b - 4 * a * c == -4 * N and abs(b) <= a <= c
    return a, b, c


def compose_coprime(f1: tuple[int, int, int], f2: tuple[int, int, int], N: int) -> tuple[int, int, int]:
    a1, b1, _ = f1
    a2, b2, _ = f2
    assert math.gcd(a1, a2) == 1
    k = (((b2 - b1) // 2) * pow(a1, -1, a2)) % a2
    A = a1 * a2
    B = b1 + 2 * a1 * k
    if B > A:
        B -= 2 * A * ((B + A - 1) // (2 * A))
    C = (B * B + 4 * N) // (4 * A)
    return reduce_form(A, B, C, N)


def make_fixture(bits: int, primes: list[int], rng: random.Random) -> tuple[int, tuple[int, int, int]]:
    lo = 3 * (1 << (bits // 2)) // 4
    hi = (1 << (bits // 2)) - 1000
    while True:
        p = next_prime(rng.randrange(lo, hi))
        q = next_prime(rng.randrange(lo, hi))
        if p != q and (p * q).bit_length() == bits:
            N = p * q
            break
    gens = []
    for pr in primes:
        if pr > 2 and N % pr != 0:
            rts = tonelli_shanks((-N) % pr, pr)
            if len(rts) == 2:
                w = rts[0]
                gens.append((pr, 2 * w, (w * w + N) // pr))
        if len(gens) >= 32:
            break
    f = (1, 0, N)
    for step in range(2000):
        g = rng.choice(gens)
        if math.gcd(f[0], g[0]) != 1:
            continue
        f = compose_coprime(f, g, N)
        if step >= 24 and 16 * f[0] * f[0] >= N and f[0] * f[0] <= 4 * N:
            return N, f
    raise RuntimeError("Failed to generate balanced reduced form")


def eval_q(form: tuple[int, int, int], u: int, v: int) -> int:
    a, b, c = form
    return a * u * u + b * u * v + c * v * v


def good_split_primes(form: tuple[int, int, int], N: int, primes: list[int]) -> list[tuple[int, tuple[int, int]]]:
    a, b, _ = form
    out = []
    for p in primes:
        if p == 2 or a % p == 0 or N % p == 0:
            continue
        rts = tonelli_shanks((-N) % p, p)
        if len(rts) == 2:
            inv2a = pow(2 * a, -1, p)
            r0 = ((-b + 2 * rts[0]) * inv2a) % p
            r1 = ((-b + 2 * rts[1]) * inv2a) % p
            if r0 > r1:
                r0, r1 = r1, r0
            assert r0 != r1 and eval_q(form, r0, 1) % p == 0 and eval_q(form, r1, 1) % p == 0
            out.append((p, (r0, r1)))
    return out


def smooth_part_and_flag(val: int, primorial: int) -> tuple[int, bool]:
    orig = val
    g = math.gcd(val, primorial % val)
    while g > 1:
        while val % g == 0:
            val //= g
        g = math.gcd(val, g)
        if val == 1:
            return orig, True
        g = math.gcd(val, primorial % val)
    return orig // val, val == 1


def is_smooth(val: int, primorial: int) -> bool:
    if val <= 1:
        return val == 1
    g = math.gcd(val, primorial % val)
    while g > 1:
        val //= g
        while val % g == 0:
            val //= g
        if val == 1:
            return True
        g = math.gcd(val, primorial % val)
    return val == 1


def reconstruct_uv(r: int, d: int, K: int) -> tuple[int, int] | None:
    old, new = d, r % d
    old_v, new_v = 0, 1
    while new > K:
        q, rem = divmod(old, new)
        old, new = new, rem
        old_v, new_v = new_v, old_v - q * new_v
    if new > 0 and 1 <= new_v <= K and math.gcd(new, new_v) == 1:
        return new, new_v
    return None


def gauss_lattice_2d(r: int, d: int) -> tuple[tuple[int, int], tuple[int, int]]:
    a = (d, 0)
    b = (r % d, 1)
    while True:
        aa = a[0] * a[0] + a[1] * a[1]
        bb = b[0] * b[0] + b[1] * b[1]
        if aa < bb:
            a, b = b, a
            aa, bb = bb, aa
        inner = a[0] * b[0] + a[1] * b[1]
        if 2 * abs(inner) <= bb:
            if b[1] < 0 or (b[1] == 0 and b[0] < 0):
                b = (-b[0], -b[1])
            return b, a
        t = (2 * inner + bb) // (2 * bb) if inner >= 0 else -((-2 * inner + bb - 1) // (2 * bb))
        a = (a[0] - t * b[0], a[1] - t * b[1])


def min_qnorm_on_lattice(form: tuple[int, int, int], N: int, r: int, d: int) -> tuple[int, tuple[int, int]]:
    """Exact minimum of Q(u,v)//d over nonzero (u,v) in Lambda_{d,r} via Gauss form reduction."""
    a, b, c = form
    r %= d
    # Basis e1 = (d, 0), e2 = (r, 1): Q(d X + r Y, Y) / d = (a*d)*X^2 + (2*a*r + b)*X*Y + (Q(r,1)//d)*Y^2
    A = a * d
    B = 2 * a * r + b
    C = (a * r * r + b * r + c) // d
    assert B * B - 4 * A * C == -4 * N
    # Track unimodular matrix along reduction to get exact minimizing (u,v)
    u1, v1 = d, 0
    u2, v2 = r, 1
    while True:
        if A > C or (A == C and B < 0):
            A, B, C = C, -B, A
            u1, v1, u2, v2 = u2, v2, -u1, -v1
        if abs(B) > A or B == -A:
            q = (B + A) // (2 * A) if B >= 0 else -((-B + A - 1) // (2 * A))
            B -= 2 * q * A
            C = (B * B + 4 * N) // (4 * A)
            u2 -= q * u1
            v2 -= q * v1
        else:
            break
    if v1 < 0 or (v1 == 0 and u1 < 0):
        u1, v1 = -u1, -v1
    assert eval_q(form, u1, v1) == A * d
    return A, (u1, v1)


def full_box_ground_truth(form: tuple[int, int, int], K: int, primes: list[int], good: list[tuple[int, tuple[int, int]]]) -> dict:
    a, b, c = form
    prim_pts = 0
    witnesses = set()
    main_shell_witnesses = set()
    good_map = {p: rts for p, rts in good}
    # Row-by-row exact sieve over all primes <= B
    for v in range(1, K + 1):
        row = [0] * (K + 1)
        for u in range(1, K + 1):
            if math.gcd(u, v) == 1:
                row[u] = a * u * u + b * u * v + c * v * v
                prim_pts += 1
        for p in primes:
            if p in good_map:
                if v % p == 0:
                    continue
                r0, r1 = good_map[p]
                s0 = (r0 * v) % p
                if s0 == 0:
                    s0 = p
                for u in range(s0, K + 1, p):
                    val = row[u]
                    if val:
                        while val % p == 0:
                            val //= p
                        row[u] = val
                s1 = (r1 * v) % p
                if s1 == 0:
                    s1 = p
                for u in range(s1, K + 1, p):
                    val = row[u]
                    if val:
                        while val % p == 0:
                            val //= p
                        row[u] = val
            else:
                for u in range(1, K + 1):
                    val = row[u]
                    if val and val % p == 0:
                        while val % p == 0:
                            val //= p
                        row[u] = val
        for u in range(1, K + 1):
            if row[u] == 1:
                witnesses.add((u, v))
                if max(u, v) > K // 4:
                    main_shell_witnesses.add((u, v))
    return {
        "primitive_points": prim_pts,
        "witnesses": witnesses,
        "main_shell_witnesses": main_shell_witnesses,
    }


def check_p1_bipartite_identity() -> list[dict]:
    """Exhaustive 40-bit check of P1: exact equality between winning branches and divisor fibers."""
    B, K = 128, 32
    primes = primes_upto(B)
    primorial = math.prod(primes)
    good_limit_lo = 2 * K * K
    good_limit_hi = 2 * B * K * K
    records = []
    for rep in range(4):
        rng = random.Random(910000 + rep)
        N, form = make_fixture(40, primes, rng)
        good = good_split_primes(form, N, primes)
        good_primes_only = [p for p, _ in good]
        good_primorial = math.prod(good_primes_only)
        truth = full_box_ground_truth(form, K, primes, good)
        # Build prime-power Hensel pieces up to good_limit_hi
        a, b, _ = form
        pieces = []
        for p, (r0, r1) in good:
            powers = []
            m = p
            rts = [r0, r1]
            while m <= good_limit_hi:
                powers.append((m, rts))
                if m * p > good_limit_hi:
                    break
                rts = [
                    (r + m * ((-(eval_q(form, r, 1) // m) * pow(2 * a * r + b, -1, p)) % p)) % (m * p)
                    for r in rts
                ]
                m *= p
            pieces.append((p, powers))

        winning_branches = 0
        total_branches = 0
        in_box_branches = 0

        def dfs(idx: int, mod: int, residues: list[int]) -> None:
            nonlocal winning_branches, total_branches, in_box_branches
            if mod > good_limit_lo:
                total_branches += len(residues)
                for r in residues:
                    pt = reconstruct_uv(r, mod, K)
                    if pt is not None:
                        in_box_branches += 1
                        if pt in truth["witnesses"]:
                            winning_branches += 1
            for j in range(idx, len(pieces)):
                p, powers = pieces[j]
                if mod * p > good_limit_hi:
                    break
                for m, rts in powers:
                    if mod * m > good_limit_hi:
                        break
                    inv = pow(mod, -1, m)
                    comb = [r + mod * (((z - r) * inv) % m) for z in rts for r in residues]
                    dfs(j + 1, mod * m, comb)

        dfs(0, 1, [0])

        # Sum exact divisor counts tau_{(2K^2, 2BK^2]}(Q_good(u,v)) over witnesses
        fiber_sum = 0
        for u, v in truth["witnesses"]:
            qval = eval_q(form, u, v)
            qgood, _ = smooth_part_and_flag(qval, good_primorial)
            # Enumerate all divisors of qgood in (good_limit_lo, good_limit_hi]
            divs = [1]
            tmp = qgood
            for p in good_primes_only:
                if tmp % p == 0:
                    mult = []
                    pw = 1
                    while tmp % p == 0:
                        tmp //= p
                        pw *= p
                        mult.append(pw)
                    divs = [d * mp for d in divs for mp in ([1] + mult)]
            assert tmp == 1
            fiber_sum += sum(1 for d in divs if good_limit_lo < d <= good_limit_hi)

        assert winning_branches == fiber_sum
        records.append({
            "rep": rep,
            "N": N,
            "form": list(form),
            "witnesses": len(truth["witnesses"]),
            "total_branches": total_branches,
            "in_box_branches": in_box_branches,
            "winning_branches": winning_branches,
            "fiber_sum": fiber_sum,
            "exact_match": winning_branches == fiber_sum,
        })
    return records


def run_rnd_branch(form, K, good, primorial, rng, mode="uniform"):
    target = 2 * K * K
    work = 0
    seen = set()
    hits = set()
    pool = [item for item in good if item[0] > (max(p for p, _ in good) // 4)] if mode == "upper" else good
    weights = [p ** (-0.85) for p, _ in pool] if mode == "tilt" else None
    t0 = perf_counter()
    while work < W_BUDGET:
        mod = 1
        r = 0
        used = set()
        while mod <= target and work < W_BUDGET:
            if weights is None:
                p, rts = rng.choice(pool)
            else:
                p, rts = rng.choices(pool, weights=weights, k=1)[0]
            if p in used:
                continue
            used.add(p)
            z = rts[rng.randrange(2)]
            if mod == 1:
                mod, r = p, z
            else:
                inv = pow(mod, -1, p)
                r = r + mod * (((z - r) * inv) % p)
                mod *= p
            work += 1
        if mod <= target:
            break
        pt = reconstruct_uv(r, mod, K)
        work += 1
        if pt is not None and pt not in seen and work < W_BUDGET:
            seen.add(pt)
            work += 1
            if is_smooth(eval_q(form, pt[0], pt[1]), primorial):
                hits.add(pt)
    elapsed = perf_counter() - t0
    main_hits = {pt for pt in hits if max(pt) > K // 4}
    return {
        "work": work,
        "tested": len(seen),
        "hits": len(hits),
        "main_hits": len(main_hits),
        "ms": elapsed * 1000.0,
    }


def run_beam(form, N, K, good, primorial, rng, use_qnorm=False):
    target = 2 * K * K
    work = 0
    seen = set()
    hits = set()
    coords = []
    t0 = perf_counter()
    while work < W_BUDGET:
        batch = rng.sample(good, min(24, len(good)))
        mod = 1
        branches = [0]
        for p, rts in batch:
            if work >= W_BUDGET:
                break
            nmod = mod * p
            inv = 1 if mod == 1 else pow(mod, -1, p)
            ext = []
            for r in branches:
                for z in rts:
                    if work >= W_BUDGET:
                        break
                    x = z if mod == 1 else r + mod * (((z - r) * inv) % p)
                    if use_qnorm:
                        score, (u, v) = min_qnorm_on_lattice(form, N, x, nmod)
                        work += 2
                    else:
                        (u, v), _ = gauss_lattice_2d(x, nmod)
                        score = max(abs(u), abs(v))
                        work += 1
                    ext.append((score, x, u, v))
            if not ext:
                break
            ext.sort(key=lambda item: item[0])
            best = ext[:4]
            branches = [item[1] for item in best]
            mod = nmod
            if mod > target:
                for _, x, _, _ in best:
                    if work >= W_BUDGET:
                        break
                    pt = reconstruct_uv(x, mod, K)
                    work += 1
                    if pt is not None and pt not in seen and work < W_BUDGET:
                        seen.add(pt)
                        coords.append(max(pt))
                        work += 1
                        if is_smooth(eval_q(form, pt[0], pt[1]), primorial):
                            hits.add(pt)
                break
    elapsed = perf_counter() - t0
    main_hits = {pt for pt in hits if max(pt) > K // 4}
    k_eff = sorted(coords)[int(0.9 * (len(coords) - 1))] if coords else max(8, K // 8)
    return {
        "work": work,
        "tested": len(seen),
        "hits": len(hits),
        "main_hits": len(main_hits),
        "k_eff": k_eff,
        "ms": elapsed * 1000.0,
    }


def run_tiny_box_ctrl(form, K, k_eff, primorial, rng):
    k_eff = max(4, min(K, k_eff))
    all_pts = [(u, v) for v in range(1, k_eff + 1) for u in range(1, k_eff + 1) if math.gcd(u, v) == 1]
    rng.shuffle(all_pts)
    work = 0
    seen = set()
    hits = set()
    t0 = perf_counter()
    for pt in all_pts:
        if work >= W_BUDGET:
            break
        seen.add(pt)
        work += 1
        if is_smooth(eval_q(form, pt[0], pt[1]), primorial):
            hits.add(pt)
    elapsed = perf_counter() - t0
    main_hits = {pt for pt in hits if max(pt) > K // 4}
    return {
        "work": work,
        "tested": len(seen),
        "hits": len(hits),
        "main_hits": len(main_hits),
        "k_eff": k_eff,
        "ms": elapsed * 1000.0,
    }


def run_mitm_collide(form, K, good, primorial, rng):
    target = 2 * K * K
    left_pool = good[0::2]
    right_pool = good[1::2]
    work = 0
    t0 = perf_counter()
    table_L = set()
    table_R = set()
    half_budget = (W_BUDGET * 4) // 10
    for pool, table in ((left_pool, table_L), (right_pool, table_R)):
        limit_w = half_budget if table is table_L else 2 * half_budget
        while work < limit_w and len(pool) >= 4:
            mod = 1
            r = 0
            used = set()
            while mod <= target and work < limit_w:
                p, rts = rng.choice(pool)
                if p in used:
                    continue
                used.add(p)
                z = rts[rng.randrange(2)]
                mod, r = (p, z) if mod == 1 else (mod * p, r + mod * (((z - r) * pow(mod, -1, p)) % p))
                work += 1
            if mod <= target:
                break
            pt = reconstruct_uv(r, mod, K)
            work += 1
            if pt is not None:
                table.add(pt)
                work += 1  # hash insert
    collisions = list(table_L & table_R)
    non_collisions = list((table_L ^ table_R))
    rng.shuffle(non_collisions)
    seen = set()
    hits = set()
    collision_hits = 0
    for pt in collisions:
        if work >= W_BUDGET:
            break
        seen.add(pt)
        work += 1
        if is_smooth(eval_q(form, pt[0], pt[1]), primorial):
            hits.add(pt)
            collision_hits += 1
    for pt in non_collisions:
        if work >= W_BUDGET:
            break
        if pt not in seen:
            seen.add(pt)
            work += 1
            if is_smooth(eval_q(form, pt[0], pt[1]), primorial):
                hits.add(pt)
    elapsed = perf_counter() - t0
    main_hits = {pt for pt in hits if max(pt) > K // 4}
    return {
        "work": work,
        "collisions_found": len(collisions),
        "collision_hits": collision_hits,
        "tested": len(seen),
        "hits": len(hits),
        "main_hits": len(main_hits),
        "ms": elapsed * 1000.0,
    }


def eval_branch_crt(factors: list[tuple[int, int]]) -> tuple[int, int]:
    mod = 1
    r = 0
    for p, z in factors:
        if mod == 1:
            mod, r = p, z
        else:
            r = r + mod * (((z - r) * pow(mod, -1, p)) % p)
            mod *= p
    return mod, r


def run_mcmc_proxy(form, K, good, primorial, small_primorial, rng):
    target = 2 * K * K
    work = 0
    seen = set()
    hits = set()
    t0 = perf_counter()
    # Initialize onto first valid in-box branch
    cur_factors = None
    cur_score = 0.0
    while work < W_BUDGET and cur_factors is None:
        facs = []
        mod = 1
        used = set()
        while mod <= target and work < W_BUDGET:
            p, rts = rng.choice(good)
            if p in used:
                continue
            used.add(p)
            facs.append((p, rts[rng.randrange(2)]))
            mod *= p
            work += 1
        if mod <= target:
            break
        _, r = eval_branch_crt(facs)
        pt = reconstruct_uv(r, mod, K)
        work += 1
        if pt is not None and work < W_BUDGET:
            seen.add(pt)
            qval = eval_q(form, pt[0], pt[1])
            work += 1
            if is_smooth(qval, primorial):
                hits.add(pt)
            cof = qval // mod
            cur_score = math.log2(math.gcd(cof, small_primorial))
            cur_factors = facs

    # Metropolis-Hastings single-prime swaps
    while work < W_BUDGET and cur_factors is not None:
        idx = rng.randrange(len(cur_factors))
        used_p = {p for i, (p, _) in enumerate(cur_factors) if i != idx}
        base_mod = math.prod(p for i, (p, _) in enumerate(cur_factors) if i != idx)
        p_new, rts = rng.choice(good)
        if p_new in used_p:
            continue
        cand = list(cur_factors)
        cand[idx] = (p_new, rts[rng.randrange(2)])
        # Top up if needed to stay > 2K^2
        mod = base_mod * p_new
        while mod <= target and work < W_BUDGET:
            p_ex, rts_ex = rng.choice(good)
            if p_ex in {p for p, _ in cand}:
                continue
            cand.append((p_ex, rts_ex[rng.randrange(2)]))
            mod *= p_ex
            work += 1
        # Trim redundant tail primes if possible
        while len(cand) > 3 and mod // cand[-1][0] > target:
            mod //= cand[-1][0]
            cand.pop()
        mod, r = eval_branch_crt(cand)
        work += len(cand)
        pt = reconstruct_uv(r, mod, K)
        work += 1
        if pt is None:
            continue
        qval = eval_q(form, pt[0], pt[1])
        cof = qval // mod
        new_score = math.log2(math.gcd(cof, small_primorial))
        work += 1
        if pt not in seen and work < W_BUDGET:
            seen.add(pt)
            work += 1
            if is_smooth(qval, primorial):
                hits.add(pt)
        if new_score >= cur_score or rng.random() < math.exp((new_score - cur_score) / 2.0):
            cur_factors = cand
            cur_score = new_score
    elapsed = perf_counter() - t0
    main_hits = {pt for pt in hits if max(pt) > K // 4}
    return {
        "work": work,
        "tested": len(seen),
        "hits": len(hits),
        "main_hits": len(main_hits),
        "ms": elapsed * 1000.0,
    }


def enumerate_sublattice_box(r: int, d: int, K: int) -> list[tuple[int, int]]:
    b1, b2 = gauss_lattice_2d(r, d)
    # Bound integer combinations x*b1 + y*b2 in [1, K]^2 using inverse of 2x2 basis (det = +-d)
    det = b1[0] * b2[1] - b1[1] * b2[0]
    assert abs(det) == d
    corners = ((1, 1), (1, K), (K, 1), (K, K))
    xs = [(c[0] * b2[1] - c[1] * b2[0]) / det for c in corners]
    ys = [(-c[0] * b1[1] + c[1] * b1[0]) / det for c in corners]
    x_lo, x_hi = math.floor(min(xs)) - 1, math.ceil(max(xs)) + 1
    y_lo, y_hi = math.floor(min(ys)) - 1, math.ceil(max(ys)) + 1
    pts = []
    for x in range(x_lo, x_hi + 1):
        u0 = x * b1[0]
        v0 = x * b1[1]
        for y in range(y_lo, y_hi + 1):
            u = u0 + y * b2[0]
            v = v0 + y * b2[1]
            if 1 <= u <= K and 1 <= v <= K and math.gcd(u, v) == 1:
                pts.append((u, v))
    return pts


def run_tree_sieve_and_scramble(form, B, K, good, primes, primorial, rng):
    lo_d = (2 * K * K) // B
    hi_d = 2 * K * K
    work = 0
    seen = set()
    hits = set()
    scramble_pairs_checked = 0
    scramble_violations = 0
    high_good = [(p, rts) for p, rts in good if p > B // 2]
    t0 = perf_counter()
    while work < W_BUDGET:
        mod = 1
        r = 0
        used = set()
        while mod <= lo_d and work < W_BUDGET:
            p, rts = rng.choice(good)
            if p in used or mod * p > hi_d:
                if mod * p > hi_d and mod > lo_d:
                    break
                if len(used) > 8:
                    break
                continue
            used.add(p)
            z = rts[rng.randrange(2)]
            mod, r = (p, z) if mod == 1 else (mod * p, r + mod * (((z - r) * pow(mod, -1, p)) % p))
            work += 1
        if not (lo_d < mod <= hi_d):
            continue
        pts = enumerate_sublattice_box(r, mod, K)
        work += 1 + len(pts)
        # Audit P5 (CRT-Scramble) on up to 20 pairs per sub-branch (uncharged diagnostic)
        if len(pts) >= 2 and scramble_pairs_checked < 500:
            for i in range(min(len(pts), 6)):
                for j in range(i + 1, min(len(pts), 6)):
                    u1, v1 = pts[i]
                    u2, v2 = pts[j]
                    scramble_pairs_checked += 1
                    for p_hi, (rt0, rt1) in high_good:
                        if p_hi in used:
                            continue
                        if (u1 - rt0 * v1) % p_hi == 0 and (u2 - rt0 * v2) % p_hi == 0:
                            scramble_violations += 1
                        if (u1 - rt1 * v1) % p_hi == 0 and (u2 - rt1 * v2) % p_hi == 0:
                            scramble_violations += 1
        # Sieve sublattice points in coordinates (x, y) or directly by evaluating cofactor and stripping good primes
        if not pts or work >= W_BUDGET:
            continue
        # 1-D/2-D sublattice sieve: map each point mod p via index, or charge exact prime-line incidences
        rems = [eval_q(form, u, v) // mod for u, v in pts]
        for p in primes:
            if work >= W_BUDGET:
                break
            for idx, (u, v) in enumerate(pts):
                if rems[idx] % p == 0:
                    work += 1
                    while rems[idx] % p == 0:
                        rems[idx] //= p
        for (u, v), rem in zip(pts, rems):
            if (u, v) not in seen:
                seen.add((u, v))
                if rem == 1:
                    hits.add((u, v))
    elapsed = perf_counter() - t0
    main_hits = {pt for pt in hits if max(pt) > K // 4}
    return {
        "work": work,
        "tested": len(seen),
        "hits": len(hits),
        "main_hits": len(main_hits),
        "scramble_pairs": scramble_pairs_checked,
        "scramble_violations": scramble_violations,
        "ms": elapsed * 1000.0,
    }


def run_box_sieve_ctrl(form, K, primes, good, rng):
    a, b, c = form
    good_map = {p: rts for p, rts in good}
    rows = list(range(1, K + 1))
    rng.shuffle(rows)
    work = 0
    tested = 0
    hits = set()
    t0 = perf_counter()
    for v in rows:
        if work >= W_BUDGET:
            break
        row = [0] * (K + 1)
        for u in range(1, K + 1):
            if math.gcd(u, v) == 1:
                row[u] = a * u * u + b * u * v + c * v * v
                tested += 1
        work += K
        for p in primes:
            if p in good_map:
                if v % p == 0:
                    continue
                r0, r1 = good_map[p]
                for r in (r0, r1):
                    s = (r * v) % p or p
                    for u in range(s, K + 1, p):
                        if row[u]:
                            work += 1
                            while row[u] % p == 0:
                                row[u] //= p
            else:
                for u in range(1, K + 1):
                    if row[u] and row[u] % p == 0:
                        work += 1
                        while row[u] % p == 0:
                            row[u] //= p
            if work >= W_BUDGET:
                break
        if work < W_BUDGET:
            for u in range(1, K + 1):
                if row[u] == 1:
                    hits.add((u, v))
    elapsed = perf_counter() - t0
    main_hits = {pt for pt in hits if max(pt) > K // 4}
    return {
        "work": work,
        "tested": tested,
        "hits": len(hits),
        "main_hits": len(main_hits),
        "ms": elapsed * 1000.0,
    }


def measure_p4_correlation(form, N, K, good, primorial, rng, max_samples=160) -> list[tuple[float, float]]:
    """Sample valid main-shell branches (max(u,v) > K//4) and pair step-1 partial norm score with terminal smooth fraction."""
    target = 2 * K * K
    pairs = []
    attempts = 0
    while len(pairs) < max_samples and attempts < 4000:
        attempts += 1
        mod = 1
        r = 0
        used = set()
        step1_score = None
        while mod <= target:
            p, rts = rng.choice(good)
            if p in used:
                continue
            used.add(p)
            z = rts[rng.randrange(2)]
            mod, r = (p, z) if mod == 1 else (mod * p, r + mod * (((z - r) * pow(mod, -1, p)) % p))
            if step1_score is None:
                qmin1, _ = min_qnorm_on_lattice(form, N, r, mod)
                step1_score = -math.log2(qmin1)
        pt = reconstruct_uv(r, mod, K)
        if pt is None or max(pt) <= K // 4:
            continue
        cof = eval_q(form, pt[0], pt[1]) // mod
        if cof <= 1:
            frac = 1.0
        else:
            sp, _ = smooth_part_and_flag(cof, primorial)
            frac = math.log2(sp) / math.log2(cof)
        pairs.append((step1_score, frac))
    return pairs


def pearson_r(pairs: list[tuple[float, float]]) -> float:
    n = len(pairs)
    if n < 3:
        return 0.0
    mx = sum(x for x, _ in pairs) / n
    my = sum(y for _, y in pairs) / n
    num = sum((x - mx) * (y - my) for x, y in pairs)
    dx = math.sqrt(sum((x - mx) ** 2 for x, _ in pairs))
    dy = math.sqrt(sum((y - my) ** 2 for _, y in pairs))
    return num / (dx * dy) if dx > 0 and dy > 0 else 0.0


def main():
    p1_records = check_p1_bipartite_identity()
    schedule = [
        (48, 256, 256),
        (60, 384, 384),
        (72, 512, 512),
        (84, 1024, 768),
        (96, 2048, 1536),
    ]
    small_primorial = math.prod(p ** 8 for p in primes_upto(31))
    all_rows = []
    pooled_corr_72_96 = []

    out_jsonl = OUT_DIR / "results.jsonl"
    with out_jsonl.open("w") as f_out:
        for rec in p1_records:
            f_out.write(json.dumps({"experiment": "P1_bipartite", **rec}) + "\n")

        for bits, B, K in schedule:
            primes = primes_upto(B)
            primorial = math.prod(primes)
            for rep in range(6):
                seed = 920000 + 1000 * bits + rep
                rng = random.Random(seed)
                N, form = make_fixture(bits, primes, rng)
                good = good_split_primes(form, N, primes)
                truth = full_box_ground_truth(form, K, primes, good)

                res = {
                    "experiment": "main_comparison",
                    "bits": bits,
                    "rep": rep,
                    "N": N,
                    "form": list(form),
                    "B": B,
                    "K": K,
                    "good_primes": len(good),
                    "primitive_points": truth["primitive_points"],
                    "truth_hits": len(truth["witnesses"]),
                    "truth_main_hits": len(truth["main_shell_witnesses"]),
                }
                res["RND_BRANCH"] = run_rnd_branch(form, K, good, primorial, random.Random(seed + 1), "uniform")
                res["RND_UPPER"] = run_rnd_branch(form, K, good, primorial, random.Random(seed + 2), "upper")
                res["PROXY_TILT"] = run_rnd_branch(form, K, good, primorial, random.Random(seed + 3), "tilt")
                res["BEAM_SV"] = run_beam(form, N, K, good, primorial, random.Random(seed + 4), False)
                res["TINY_BOX_CTRL"] = run_tiny_box_ctrl(
                    form, K, res["BEAM_SV"]["k_eff"], primorial, random.Random(seed + 5)
                )
                res["BEAM_QNORM"] = run_beam(form, N, K, good, primorial, random.Random(seed + 6), True)
                res["MITM_COLLIDE"] = run_mitm_collide(form, K, good, primorial, random.Random(seed + 7))
                res["MCMC_PROXY"] = run_mcmc_proxy(
                    form, K, good, primorial, small_primorial, random.Random(seed + 8)
                )
                res["TREE_SIEVE"] = run_tree_sieve_and_scramble(
                    form, B, K, good, primes, primorial, random.Random(seed + 9)
                )
                res["BOX_SIEVE"] = run_box_sieve_ctrl(form, K, primes, good, random.Random(seed + 10))

                corr_pairs = measure_p4_correlation(form, N, K, good, primorial, random.Random(seed + 11))
                res["step1_corr_r"] = pearson_r(corr_pairs)
                res["step1_corr_n"] = len(corr_pairs)
                if bits >= 72:
                    pooled_corr_72_96.extend(corr_pairs)

                all_rows.append(res)
                f_out.write(json.dumps(res) + "\n")
                f_out.flush()

    # Summarize all pre-registered metrics
    methods = [
        "RND_BRANCH",
        "RND_UPPER",
        "PROXY_TILT",
        "BEAM_SV",
        "TINY_BOX_CTRL",
        "BEAM_QNORM",
        "MITM_COLLIDE",
        "MCMC_PROXY",
        "TREE_SIEVE",
        "BOX_SIEVE",
    ]
    by_bits = {}
    for bits, B, K in schedule:
        sub = [r for r in all_rows if r["bits"] == bits]
        by_bits[str(bits)] = {
            "truth_hits": sum(r["truth_hits"] for r in sub),
            "truth_main_hits": sum(r["truth_main_hits"] for r in sub),
            "methods": {
                m: {
                    "hits": sum(r[m]["hits"] for r in sub),
                    "main_hits": sum(r[m]["main_hits"] for r in sub),
                    "tested": sum(r[m]["tested"] for r in sub),
                    "work": sum(r[m]["work"] for r in sub),
                    "ms": round(sum(r[m]["ms"] for r in sub), 2),
                }
                for m in methods
            },
        }

    large_rows = [r for r in all_rows if r["bits"] >= 72]
    pooled_large = {
        m: {
            "hits": sum(r[m]["hits"] for r in large_rows),
            "main_hits": sum(r[m]["main_hits"] for r in large_rows),
            "tested": sum(r[m]["tested"] for r in large_rows),
            "work": sum(r[m]["work"] for r in large_rows),
        }
        for m in methods
    }
    p5_pairs = sum(r["TREE_SIEVE"]["scramble_pairs"] for r in all_rows)
    p5_viol = sum(r["TREE_SIEVE"]["scramble_violations"] for r in all_rows)
    pooled_r = pearson_r(pooled_corr_72_96)

    summary = {
        "P1_all_exact": all(r["exact_match"] for r in p1_records),
        "P1_records": p1_records,
        "by_bits": by_bits,
        "pooled_72_84_96": pooled_large,
        "P4_pooled_r_72_96": round(pooled_r, 5),
        "P4_pooled_n_72_96": len(pooled_corr_72_96),
        "P5_scramble_pairs": p5_pairs,
        "P5_scramble_violations": p5_viol,
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
