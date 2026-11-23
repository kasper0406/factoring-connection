from config import *
import math, hmac, random, json, tempfile, sys
import gmpy2 as g

sys.path.append("/Volumes/git/factoring-connection/.venv/lib/python3.12/site-packages")
from cypari2 import Pari
from flint import fmpz_mat
from screen import Dataset, OPS, BASE, TARGETS, pack, pair_bins, CAT, RES
from confirm_eval import locks, jobs

pari = Pari()


def direct_cf(p, q, k, length):
    N = p * q
    D = k * N
    s = math.isqrt(D)
    a = s
    P = 0
    Q = 1
    mx = even = odd = 0
    kr = 1
    A0, A = 0, 1
    hit = False

    def audit(*args):
        nonlocal hit
        hit |= any(1 < math.gcd(int(y), N) < N for y in args)

    for i in range(1, length + 1):
        v = (a * A) % N
        A0, A = A, (v + A0) % N
        audit(v, A)
        P = Q * a - P
        Q = (D - P * P) // Q
        a = (s + P) // Q
        audit(P, Q, s + P, a)
        mx = max(mx, a)
        z = Q
        for t in [2, 3, 5, 7]:
            while z % t == 0:
                z //= t
        if math.isqrt(z) ** 2 == z:
            if i % 2:
                odd += 1
            else:
                even += 1
        r = math.isqrt(Q)
        if r * r == Q:
            audit(A - r, A + r)
        kr *= int(g.legendre(Q, 11))
    vals = [Q, mx, even + 128 * odd, kr]
    audit(*vals)
    return vals, int(hit)


def ext(N, d, rel, base, e):
    def mul(a, b):
        t = [0] * (2 * d - 1)
        for i in range(d):
            for j in range(d):
                t[i + j] += a[i] * b[j]
        for i in range(len(t) - 1, d - 1, -1):
            for j, f in enumerate(rel):
                t[i - d + j] += t[i] * f
        return [x % N for x in t[:d]]

    out = [1] + [0] * (d - 1)
    while e:
        if e % 2:
            out = mul(base, out)
        base = mul(base, base)
        e //= 2
    return out


def pointmul(P, e, prime, a2, a4):
    def add(P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        x, y = P
        u, v = Q
        if x == u and (y + v) % prime == 0:
            return None
        slope = (
            (3 * x * x + 2 * a2 * x + a4) * pow(2 * y, -1, prime)
            if P == Q
            else (v - y) * pow(u - x, -1, prime)
        ) % prime
        w = (slope * slope - a2 - x - u) % prime
        return w, (slope * (x - w) - y) % prime

    out = None
    while e:
        if e % 2:
            out = add(out, P)
        P = add(P, P)
        e //= 2
    return out


def ringcheck(n, a):
    src = DATA / f"audit{n}.input"
    a.tofile(src)
    ids = [
        0,
        1,
        2,
        3,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        15,
        16,
        18,
        19,
        21,
        23,
        24,
        25,
        27,
        29,
        31,
        32,
        35,
        43,
        44,
        46,
    ]
    res = run(
        ["./engine", "dump", ",".join(map(str, ids)), src, "/dev/null"],
        capture_output=True,
        text=True,
    )
    vals = np.array(
        [[int(z) for z in l.split()] for l in res.stdout.splitlines()], dtype=object
    ).reshape(len(ids), len(a), 6)
    ct = 0
    infra = []
    for ix, (pp, qq) in enumerate(a):
        p, q = map(int, (pp, qq))
        N = p * q
        s = math.isqrt(N)
        for v, vv in zip(ids, vals[:, ix]):
            o = list(map(int, vv[:4]))
            sing = vv[5]
            x = pow(2, N, N)
            y = pow(3, N, N)
            if v == 0:
                assert o[0] == pow(x + pow(3, (N - 1) // 2, N) + 1, N, N)
            elif v == 1:
                assert o[0] == (pow(x + 1, N, N) + pow(y - 1, N, N)) % N
            elif v == 2:
                assert o[0] == (pow(pow(x + 1, N, N) + 2, N, N) - y) % N
            elif v == 3 and not sing:
                assert o[0] == pow(pow(x + 1, -1, N) + 2, N, N)
            elif v == 6:
                assert o[0] == pow(pow(x + 1, math.lcm(*range(1, 65)), N) + 3, N, N)
            elif v == 7:
                assert o[0] == (pow(2, 2**10, N) + 3) % N
            elif v in (8, 9, 10, 11, 23, 24, 25):
                d = {8: 3, 9: 4, 10: 5, 11: 2, 23: 2, 24: 2, 25: 4}[v]
                rel = {
                    8: [1, 1, 0],
                    9: [1, 1, 0, 0],
                    10: [1, 0, 0, 0, 0],
                    11: [-3, 5],
                    23: [-3, 1],
                    24: [-5, 3],
                    25: [pow(3, -1, N), -1, 2, 0],
                }[v]
                b = [2 if v == 10 else 0, 1] + [0] * (d - 2)
                expect = ext(N, d, rel, b, N + (1 if v == 23 else -1 if v == 24 else 0))
                assert o[: min(4, d)] == expect[: min(4, d)], (n, v)
            elif v in (12, 13, 15):
                k = 3 if v == 13 else 1
                e = N if v == 15 else 1 << (n // 2)
                l = 3
                while True:
                    if g.is_prime(l):
                        roots = [r for r in range(1, l) if (r * r + k * N) % l == 0]
                        if roots:
                            break
                    l += 2
                b = 2 * min(roots)
                f = pari.Qfb(l, b, (b * b + 4 * k * N) // (4 * l))
                f = pari.qfbnupow(f, e)
                assert o[:3] == [int(f[i]) for i in range(3)], (n, v, o, f)
            elif v in (16, 18, 19, 46):
                k = {16: 1, 18: 5, 19: 7, 46: 11}[v]
                f, h = direct_cf(p, q, k, 16 if v == 16 else 64)
                assert o == f and vv[4] == h, (n, v, o, f, h, vv[4])
            elif v == 21:
                aa, bb, cc = o[:3]
                assert bb * bb - 4 * aa * cc == 4 * N
                if n == 40:
                    with g.context(g.get_context(), precision=192):
                        sd = g.sqrt(4 * N)
                        a0, b0, c0 = 1, 2 * s, s * s - N
                        dist = 0
                        for _ in range(20000):
                            step = g.log(abs((b0 + sd) / (b0 - sd))) / 2
                            if dist + step > math.isqrt(s):
                                break
                            dist += step
                            t = max(2 * s, abs(c0))
                            rem = (t + b0) % abs(2 * c0)
                            quo = (t + b0 - rem) // (2 * c0)
                            a0, b0, c0 = c0, t - rem, a0 - quo * (b0 - quo * c0)
                        assert (abs(aa), bb, abs(cc)) == (abs(a0), b0, abs(c0)), (
                            o[:3],
                            (a0, b0, c0),
                            n,
                        )
            elif v in (27, 29):
                X = math.isqrt(s)
                A = (
                    [[N, 0, 0], [s, X, 0], [0, s * X, X**2]]
                    if v == 27
                    else [[N, 0, 0], [x, X, 0], [y, 0, X]]
                )
                F = fmpz_mat(A).lll()
                b = [int(F[0, i]) for i in range(3)]
                assert (
                    abs(o[0]) == abs(b[0])
                    and o[1] == sum(i * i for i in b)
                    and o[2] == sum(map(abs, b))
                )
            elif v in (31, 32, 35, 43, 44) and not sing:
                for prime in (p, q):
                    if v >= 43:
                        z = pointmul(
                            (0, pow(2, -1, prime)), N + (v == 44), prime, 0, -1
                        )
                        expect = None if z is None else z[0]
                        value = o[3]
                    else:
                        x0 = 3 if v == 32 else 2
                        A = 4 * ((N >> 7) % 17) + 6 if v == 35 else 6
                        B = (x0**3 + A * x0**2 + x0) % prime
                        z = pointmul(
                            (B * x0 % prime, B * B % prime),
                            N + (v == 32),
                            prime,
                            A * B % prime,
                            B * B % prime,
                        )
                        expect = None if z is None else z[0] * pow(B, -1, prime) % prime
                        value = o[2]
                    assert expect == value % prime, (n, v, prime, expect, value % prime)
            ct += 1
    return ct


def external(p, q, n):
    N = p * q
    fl = int(q - p <= math.isqrt(math.isqrt(N << 32)))
    for k in range(2, 33):
        a = 1 + math.isqrt(4 * k * N - 1)
        for u in range(a, a + 8):
            b = math.isqrt(u * u - 4 * k * N)
            if b * b == u * u - 4 * k * N and any(
                1 < math.gcd(v, N) < N for v in (u - b, u + b)
            ):
                fl |= 2

    def smooth(z):
        for t in range(2, 65):
            while z % t == 0:
                z //= t
        return z == 1

    if smooth(p - 1) or smooth(q - 1):
        fl |= 4
    if smooth(p + 1) or smooth(q + 1):
        fl |= 8
    return fl


def main():
    arithmetic = 0
    sample = 0
    sc = 0
    target_checks = 0
    for n in (40, 48, 56, 64, 80, 96):
        a = (
            load(n, "D")[:6]
            if n < 64
            else np.memmap(jobs(n)[0], "<u8", "r").reshape(-1, 2)[:6]
        )
        arithmetic += ringcheck(n, a)
        for fn in jobs(n)[:: max(1, len(jobs(n)) // 8)]:
            a = np.memmap(fn, "<u8", "r").reshape(-1, 2)
            b = np.fromfile(str(fn) + ".budget", "u1")
            locked = np.memmap(str(fn) + ".locked", "<u4", "r").reshape(-1, 4)
            rng = np.random.default_rng(seed("independent-audit", n, fn.name))
            for i in rng.choice(len(a), 48, replace=False):
                p, q = map(int, a[i])
                N = p * q
                assert g.is_probab_prime(p) == 2 and g.is_probab_prime(q) == 2
                for prime in (p, q):
                    digest = hmac.digest(
                        b"ent-c-R2-v1/prime-pool", prime.to_bytes(8, "little"), "sha256"
                    )
                    assert int.from_bytes(digest[:2], "big") >= 39322
                assert (b[i] & 15) == external(p, q, n)
                sample += 1
                cache = {}
                for j, row in enumerate(locks(n)):
                    t = int(row["test"])
                    op, v, slot, bit = OPS[t]
                    assert op == 0 and slot == 2
                    if v not in cache:
                        cache[v] = direct_cf(p, q, {18: 5, 19: 7, 46: 11}[v], 64)
                    vals, hit = cache[v]
                    z = vals[slot]
                    code = (
                        (z & 63)
                        | (((z % 256).bit_count() % 2) << 6)
                        | (z.bit_count() % 2 << 7)
                        | (int(g.jacobi(z, N) < 0) << 8)
                        | (int(g.jacobi(z, N) == 0) << 9)
                    )
                    code |= ((8 * (z % N)) // N) << 10
                    code |= ((z % N).bit_count() > n // 2) << 13
                    code |= (
                        (
                            z & ((0xE1990AAFA5BABC65 << 64) | 0xD5A65BB99D13A6AD)
                        ).bit_count()
                        % 2
                    ) << 15
                    assert (locked[i, 1] >> j & 1) == (code >> bit & 1) and (
                        locked[i, 2] >> j & 1
                    ) == hit
                    sc += 1
        if n <= 56:
            ds = Dataset(n, "D")
            l, r = pair_bins(ds.a, n, "D", ds.clean)
            counts = np.load(RES / f"D{n}_screen.npz")["counts"]
            for t in np.random.default_rng(seed("count-audit", n)).choice(
                len(CAT), 64, replace=False
            ):
                xx = ds.x(t).astype("i1")
                for j in range(len(TARGETS)):
                    yy = ds.y[j].astype("i1")
                    ok = ds.ok[j, l] & ds.ok[j, r]
                    di = (xx[l] - xx[r]) * (yy[l] - yy[r])
                    assert list(counts[t, j]) == [
                        int((ok & (di == 1)).sum()),
                        int((ok & (di == -1)).sum()),
                    ]
                    target_checks += 1
    res = dict(
        independent_variant_cases=arithmetic,
        budget_and_prime_audits=sample,
        locked_prediction_and_path_checks=sc,
        independent_discovery_contingencies=target_checks,
    )
    (RES / "verification.json").write_text(json.dumps(res, indent=2))
    print(res)


if __name__ == "__main__":
    main()
