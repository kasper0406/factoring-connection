from inference import *


def witness(N, p, q):
    D = 7 * N
    s = math.isqrt(D)
    a = s
    P = 0
    Q = 1
    A0, A = 0, 1
    B0, B = 1, 0
    found = []
    for i in range(1, 65):
        A0, A = A, a * A + A0
        B0, B = B, a * B + B0
        P = a * Q - P
        Q = (D - P * P) // Q
        a = (P + s) // Q
        sign = 1 if i % 2 == 0 else -1
        assert A * A - D * B * B == sign * Q
        h = Q
        core = 1
        for prime in (2, 3, 5, 7):
            e = 0
            while h % prime == 0:
                h //= prime
                e += 1
            if e % 2:
                core *= prime
        if math.isqrt(h) ** 2 == h:
            assert (
                pow((sign * core) % p, (p - 1) // 2, p) == 1
                and pow((sign * core) % q, (q - 1) // 2, q) == 1
            )
            found.append((i, core))
    assert sum(i % 2 == 0 for i, core in found) % 2 == 1
    return found


def main():
    out = []
    for n in (64, 80, 96):
        count = 0
        freq = Counter()
        total = 0
        for fn in jobs(n):
            h, x, fl, c = (
                np.memmap(str(fn) + ".locked", dtype="<u4", mode="r").reshape(-1, 4).T
            )
            b = np.memmap(str(fn) + ".budget", dtype="u1", mode="r")
            a = np.memmap(fn, dtype="<u8", mode="r").reshape(-1, 2)
            good = (
                (x & 1 != 0) & (fl & 1 == 0) & (RESIDUES[h // 1024] % 8 == 1) & (b == 0)
            )
            for p, q in a[np.flatnonzero(good)[: 256 - count]]:
                p, q = int(p), int(q)
                w = witness(p * q, p, q)
                for i, core in w:
                    if i % 2 == 0:
                        freq[core] += 1
                        total += 1
                count += 1
            if count >= 256:
                break
        for core, ct in sorted(freq.items()):
            out.append(
                dict(
                    n=n,
                    positive_rows_checked=count,
                    even_certificates=total,
                    squarefree_core=core,
                    count=ct,
                )
            )
    write_records("certificate_checks.tsv", out)
    print("CF identity and residue-square checks passed", flush=True)


if __name__ == "__main__":
    main()
