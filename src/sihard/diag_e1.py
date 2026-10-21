"""POST-HOC (unregistered) diagnostics for sihard E1. Run: .venv/bin/python -m src.sihard.diag_e1"""
import random, statistics
from src.common.instances import balanced_semiprime
from src.sihard.core import smooth_mask, primes_upto

def part1():
    y = 1 << 10; ps = primes_upto(y)
    for seed in (0, 1):
        N = balanced_semiprime(48, seed).N
        h = 18000
        K = 1500
        near = [int(smooth_mask(k*N, h+1, ps).sum()) for k in range(1, K+1)]
        half = [int(smooth_mask(k*N + N//2, h+1, ps).sum()) for k in range(1, K+1)]
        rng = random.Random(seed)
        ctrl = [int(smooth_mask(rng.randrange(N, (K+1)*N), h+1, ps).sum()) for _ in range(K)]
        # size-matched control: x = kN + random offset in [0, N)
        ctrl2 = [int(smooth_mask(k*N + rng.randrange(N - h), h+1, ps).sum()) for k in range(1, K+1)]
        m = lambda a: statistics.mean(a)
        print(seed, "near", m(near), "kN+N/2", m(half), "ctrl", m(ctrl), "ctrl2(kN+rand)", m(ctrl2))
        # dispersion within narrow k-bins (remove size heterogeneity)
        for name, arr in (("near", near), ("ctrl2", ctrl2)):
            ids = []
            for b in range(0, K, 100):
                seg = arr[b:b+100]; mu = m(seg)
                if mu > 0: ids.append(statistics.pvariance(seg)/mu)
            print("   binned ID", name, round(m(ids), 3))

def part2():
    y = 1 << 12; ps = primes_upto(y); h = 1509; K = 3000
    rng = random.Random(7)
    mods = [("semiprime0", balanced_semiprime(48, 0).N), ("semiprime5", balanced_semiprime(48, 5).N)]
    mods += [("random_odd%d" % i, rng.randrange(2**47, 2**48) | 1) for i in range(3)]
    mods += [("random_even", rng.randrange(2**47, 2**48) & ~1)]
    for name, N in mods:
        near = [int(smooth_mask(k*N, h+1, ps).sum()) for k in range(1, K+1)]
        c2 = [int(smooth_mask(k*N + rng.randrange(N - h), h+1, ps).sum()) for k in range(1, K+1)]
        a, b = statistics.mean(near), statistics.mean(c2)
        se = (a/K + b/K) ** 0.5
        print(f"{name:12s} N%2={N%2} near={a:.3f} ctrl2={b:.3f} ratio={a/b:.3f} z={(a-b)/se:.1f}", flush=True)

if __name__ == "__main__":
    part1()
    part2()
