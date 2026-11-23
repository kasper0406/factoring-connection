"""Independent arithmetic, population, contingency and interval audits."""

import math, itertools
from collections import Counter
import gmpy2 as g
from scipy.stats import binomtest
from controls import *

REASONS = [
    "integer_gate",
    "lucas_ladder",
    "CF_state",
    "Fermat",
    "Lehman",
    "p_minus_1",
    "p_plus_1",
]


def direct(p, q, n):
    N = p * q
    flags = 0
    gate = []

    def check(x, bit):
        nonlocal flags
        if 1 < math.gcd(int(x), N) < N:
            flags |= bit

    def add(v, mod=False):
        v = int(v)
        gate.append(v)
        check(v, 1)
        if mod:
            check(v - 1, 1)
            check(v + 1, 1)

    add(N)
    for e in (N - 1, N + 1, (N - 1) // 2):
        check(e, 1)
    for x in (2, 3, 5, 7):
        for e in (N, (N - 1) // 2, N - 1, N + 1):
            add(pow(x, e, N), True)
    for P, Q in ((3, 1), (5, 1), (3, -1)):
        add(g.lucasv_mod(P, Q, N, N), True)
        ix = 0
        for digit in bin(N)[2:]:
            ix = 2 * ix + int(digit)
            for j in (ix, ix + 1):
                check(g.lucasv_mod(P, Q, j, N), 2)
    for x in (2, 3):
        for c in (-1, 1):
            shifted = pow(x, N, N) + c
            check(shifted, 1)
            add(pow(shifted, N, N), True)
    for k in (3, 5, 7, 11, 13, 16, 31):
        a, b = divmod(N, k)
        add(a)
        add(b)
    root = math.isqrt(N)
    add(root)
    for u in (N, 2 * N * N, 3 * N * N):
        scaled = math.isqrt(u * 2**32)
        check(scaled, 1)
        check(scaled // 65536, 1)
        add(scaled % 65536)
    m = 0
    d = 1
    a = root
    for _ in range(4):
        m = d * a - m
        d = (N - m * m) // d
        check(m, 4)
        check(d, 4)
        check(root + m, 4)
        a = (root + m) // d
        add(a)
    u, v = N, 2 ** (n // 2)
    for i in range(5):
        a, b = divmod(u, v) if v else (0, 0)
        check(b, 4)
        if i:
            add(a)
        else:
            check(a, 4)
        u, v = v, b

    def square_scan(k, start, steps, bit):
        for a in range(start, start + steps):
            b = math.isqrt(a * a - k)
            if b * b == a * a - k:
                check(a - b, bit)
                check(a + b, bit)

    square_scan(N, 1 + math.isqrt(N - 1), 256, 8)
    for k in range(2, 33):
        square_scan(4 * k * N, 1 + math.isqrt(4 * k * N - 1), 8, 16)
    E = 1
    for prime in (e for e in range(2, 129) if g.is_prime(e)):
        pp = prime
        while pp * prime <= 128:
            pp *= prime
        E *= pp
        check(pow(2, E, N) - 1, 32)
        check(g.lucasv_mod(3, 1, E, N) - 2, 64)
    enc = []
    for v in gate:
        j = g.jacobi(v, N)
        enc.append(
            v % 256
            | ((v % 256).bit_count() % 2) << 8
            | ((v % 65536).bit_count() % 2) << 9
            | (v.bit_count() % 2) << 10
            | (int(j < 0) << 11)
            | (int(j == 0) << 12)
        )
    return flags, enc


def tvalues(catalog_ids, unary):
    vals = []
    for i in catalog_ids:
        op, a, b = SPECS[i]
        if op == 0:
            v = unary[i]
        elif op == 1:
            v = unary[a] ^ unary[b]
        elif op == 2:
            v = unary[a] & unary[b]
        else:
            v = unary[a] | unary[b]
        vals.append(v)
    return np.array(vals)


def run():
    auditcases = 0
    scorecases = 0
    data_stats = []
    overlap = []
    split_stats = []
    signatures = []
    overlap_stats = []
    d_hits = []
    allprime = 0
    totalrows = 0
    for n in SIZES:
        sets = {}
        primes = []
        for s in SPLITS:
            arr = rows(n, s)
            M = len(arr)
            totalrows += M
            assert arr.min() >= 2 ** (n // 2 - 1) and arr.max() < 2 ** (n // 2)
            assert np.all(arr[:, 0] < arr[:, 1]) and np.all(arr[:, 2] < arr[:, 3])
            assert len(arr) == SPLITS[s]
            p, q, r, t = arr.T
            y, keep, pref = context(n, s)
            assert pref.min() >= 128
            assert np.all(y[17])  # 4<=floor(8p/q)<=7 for balanced factors.
            primes.append(np.unique(arr))
            unique = np.ascontiguousarray(arr[:, :2]).view("V16").ravel()
            sets[s] = np.unique(unique)
            dups = M - len(sets[s])
            same = int(np.all(arr[:, :2] == arr[:, 2:], axis=1).sum())
            data_stats.append(
                dict(
                    n=n,
                    split=s,
                    primary=M,
                    refs=M,
                    distinct_primary=len(sets[s]),
                    repeated_draws=dups,
                    identical_reference=same,
                    occupied_prefixes=len(np.unique(pref)),
                )
            )
            flags = np.fromfile(str(stem(n, s)) + ".f", "u1").reshape(-1, 2)
            for col, name in enumerate(["primary", "reference"]):
                f = flags[:, col]
                z = Counter(map(int, f))
                overlap_stats.extend(
                    dict(n=n, split=s, member=name, mask=k, count=v)
                    for k, v in sorted(z.items())
                )
                split_stats.append(
                    dict(
                        n=n,
                        split=s,
                        member=name,
                        m=M,
                        flagged=int((f > 0).sum()),
                        **{
                            key: int(((f & 1 << i) != 0).sum())
                            for i, key in enumerate(REASONS)
                        },
                        both_unflagged=int((flags == 0).all(axis=1).sum()),
                    )
                )
            columns = np.memmap(str(stem(n, s)) + ".g", "<u2", "r", shape=(M, 50))
            rng = np.random.default_rng(seed("audit", n, s))
            ix = set(map(int, rng.choice(M, size=24, replace=False)))
            for i in range(7):
                has = np.flatnonzero(flags[:, 0] & (1 << i))
                if len(has):
                    ix.add(int(has[0]))
            for row in sorted(ix):
                for j in (0, 1):
                    fp, fq = map(int, arr[row, 2 * j : 2 * j + 2])
                    flag, gates = direct(fp, fq, n)
                    assert flag == flags[row, j], (
                        n,
                        s,
                        row,
                        j,
                        fp,
                        fq,
                        flag,
                        flags[row, j],
                    )
                    if j == 0:
                        assert np.array_equal(gates, columns[row]), (
                            n,
                            s,
                            row,
                            gates,
                            columns[row],
                        )
                    auditcases += 1
                for j, ell in enumerate((3, 5, 7, 11)):
                    assert y[j + 9, row] == (g.legendre(ell, int(p[row])) < 0)
            un = unaries(n, s)
            a = np.load(R / f"{n}_{s}_counts.npz")
            cs = a["counts"]
            tasks = a["pairs"]
            if s == "D":
                sig = tvalues(range(len(CAT)), un[:, :32])
                signatures.append(sig.copy())
                d_hits.append(
                    dict(
                        n=n,
                        matched=int((agreement_p(cs[:, 2]) <= THRESH).sum()),
                        matched_clean=int((agreement_p(cs[:, 3]) <= THRESH).sum()),
                    )
                )
            stats = np.load(R / f"{n}_{s}_stats.npz")["stats"]
            for i in rng.choice(len(tasks), size=32, replace=False):
                test, j = tasks[i]
                testbits = np.unpackbits(
                    tvalues([test], un)[0].view("u1"), bitorder="little"
                )[:M].astype(bool)
                for v in range(4):
                    ix = keep[j, v]
                    table = np.bincount(
                        testbits[ix].astype(int) * 2 + y[j, ix], minlength=4
                    )
                    assert np.array_equal(cs[i, v], table)
                    m = int(ix.sum())
                    hit = int(table[0] + table[3])
                    ss = stats[i, v]
                    if m:
                        interval = binomtest(hit, m).proportion_ci(method="exact")
                        assert np.allclose(
                            ss[1:3],
                            [interval.low - 0.5, interval.high - 0.5],
                            atol=2e-10,
                            rtol=0,
                        )
                        joint = table.reshape(2, 2) / m
                        frac = joint.sum(1)[:, None] * joint.sum(0)
                        jj = joint > 0
                        mi = (joint[jj] * np.log2(joint[jj] / frac[jj])).sum()
                        assert (
                            abs(ss[3] - mi) < 1e-12
                            and ss[4] - 1e-12 <= mi <= ss[5] + 1e-12
                        )
                    else:
                        assert np.isnan(ss[0]) and np.isnan(ss[3])
                    scorecases += 1
        primes = np.unique(np.concatenate(primes))
        allprime += len(primes)
        assert all(g.is_probab_prime(int(p)) == 2 for p in primes)
        for a, b in itertools.combinations(SPLITS, 2):
            inter = len(np.intersect1d(sets[a], sets[b], assume_unique=True))
            raw = rows(n, b)
            ub = np.ascontiguousarray(raw[:, :2]).view("V16").ravel()
            overlap.append(
                dict(
                    n=n,
                    earlier=a,
                    later=b,
                    shared_distinct=inter,
                    later_draws_seen_earlier=int(np.isin(ub, sets[a]).sum()),
                    later_rows=len(raw),
                )
            )
    full = np.ascontiguousarray(np.column_stack(signatures))
    unique = len(np.unique(full.view(np.dtype(("V", full.shape[1] * 8))).ravel()))
    p24 = [p for p in range(2048, 4096) if g.is_prime(p)]
    population = sum(p * q >= 2**23 for p, q in itertools.combinations(p24, 2))
    output = dict(
        primary_rows=totalrows,
        reference_rows=totalrows,
        prime_values_verified=allprime,
        arithmetic_reference_cases=auditcases,
        contingency_and_interval_checks=scorecases,
        admissible_semiprime_witness_rows=2048 * len(SIZES),
        distinct_on_semiprime_witnesses=unique,
        exact_24bit_population=population,
    )
    for name, data in [
        ("sample_counts", data_stats),
        ("overlap", overlap),
        ("split_flags", split_stats),
        ("flag_mask_counts", overlap_stats),
        ("discovery_by_size", d_hits),
    ]:
        write_rows(name + ".tsv", data)
    json.dump(output, open(R / "audit.json", "w"), indent=2)
    print(output)


if __name__ == "__main__":
    run()
