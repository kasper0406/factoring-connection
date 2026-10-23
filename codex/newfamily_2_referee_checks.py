"""Independent referee checks. All output is confined to this scratch directory."""
import cmath
import gzip
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from sympy import primerange

from src.newfamily.recip import cf, dedekind6, jacobi
from src.newfamily.battery import arguments
from src.newfamily.diag import dedekind_local_part, factor_from_local_part
from src.eps.analysis import strat_z

ROOT = Path('/Volumes/git/factoring-connection')
out = {}
count = 0
for c in range(3, 401, 2):
    for a in range(1, c):
        if math.gcd(a, c) != 1:
            continue
        d = 2 * dedekind6(a, c)
        qs = cf(a, c)
        alt = sum(x if i % 2 == 0 else -x for i, x in enumerate(qs))
        assert Fraction(d, c) == alt + Fraction(a + pow(a, -1, c), c) - 2 + (-1) ** len(qs)
        assert (d - a - pow(a, -1, c)) % c == 0
        assert (d - c - 1 + 2 * jacobi(a, c)) % 8 == 0
        if c % 3:
            assert d % 3 == 0
        count += 1
out['BHK_and_congruences'] = count

count = 0
ps = list(primerange(3, 45))
for i, p in enumerate(ps):
    for q in ps[i + 1:]:
        N = p * q
        for a in range(1, N):
            if math.gcd(2 * a * (a * a + 1), N) != 1:
                continue
            L = dedekind_local_part(a, p, q)
            assert factor_from_local_part(N, a, L) in (p, q)
            U = Fraction(dedekind6(a, N), 6 * N) - L
            assert (6 * N * U).denominator == 1
            assert (int(6 * N * U) - (a * a + 1) * (p - 1) * (q - 1) * pow(2 * a, -1, N)) % N == 0
            count += 1
out['DSF_including_p3'] = count
out['DSF_excluded_a2_minus1'] = dict(N=65, a=8, L=str(dedekind_local_part(8, 5, 13)))

errs = []
for N in [15, 21, 33, 35, 39, 55, 77]:
    M = N * N
    units = [(x, pow(x, -1, M)) for x in range(M) if math.gcd(x, M) == 1]
    for a, b in [(1, 1), (4, 1), (4, 9), (2, 1)]:
        if math.gcd(a * b, N) != 1:
            continue
        val = sum(cmath.exp(2j * math.pi * ((a * x + b * xi) % M) / M) for x, xi in units)
        closed = N * sum(cmath.exp(4j * math.pi * v / M) for v in range(M) if (v * v - a * b) % M == 0)
        errs.append(abs(val - closed))
out['K_square_max_error'] = max(errs)

N = 193 * 233
out['gold_counterexample'] = dict(p=193, q=233, N=N, a=arguments(N)['gold'], gcd=math.gcd(arguments(N)['gold'], N), Nmod4=N % 4)

data = {}
for n in [32, 48, 64, 128]:
    summary = json.loads((ROOT / f'experiments/analysis/newfamily_main_n{n}.json').read_text())
    data[n] = {t: {k: summary[t][k] for k in ['n_inst', 'planted', 'shuffle_mean_z2', 'menu_mean_z2', 'menu_max_abs_z', 'mde_rho', 'hits', 'degenerate']} for t in ['eps', 'qr2']}
out['saved_validity'] = data

# Exact finite-prime population: one row per unordered, distinct pair, conditional on N = 1 mod 4.
ps = np.array(list(primerange(3 * 2**14, 2**16)), dtype=np.int64)
i, j = np.triu_indices(len(ps), k=1)
p, q = ps[i], ps[j]
mask = p % 4 == q % 4
p, q = p[mask], q[mask]
N = p * q
y = (p % 4 == 1).astype(float)
ds5 = {r: dedekind6(r, 5) for r in range(1, 5)}
# Reciprocity: 10 D(5,N) = 25 + N^2 + 1 - 15N - 2N D(N mod 5, 5).
D5 = np.array([(int(v)**2 + 26 - 15*int(v) - 2*int(v)*ds5[int(v) % 5]) // 10 for v in N], dtype=object)
stats = {'nc_lead16': N >> 16, 'nc_ded5_bits': np.array([int(x).bit_length() for x in abs(D5)])}
pop = {'n_primes': len(ps), 'n_pairs': len(N), 'eps_plus': float(y.mean()), 'statistics': {}}
for name, x in stats.items():
    z = strat_z(x, y, N % 256)
    pop['statistics'][name] = dict(z_population=z['z'], predicted_z_200000=z['z']*math.sqrt(200000/len(N)))
out['exact_D32_population'] = pop

raw = {}
for n in [32, 48, 64, 128]:
    nus = []
    expected = 0.0
    count = 0
    with gzip.open(ROOT / f'experiments/raw/newfamily_main_n{n}.jsonl.gz', 'rt') as fh:
        header = json.loads(next(fh))
        for line in fh:
            r = json.loads(line)
            count += 1
            expected += 1/r['p'] + 1/r['q'] - 2/r['N']
            nus.extend(r['nonunits'])
    raw[n] = dict(count=count, nonunit_events=nus, expected_uniform_unit_failures=expected)
out['raw_counts'] = raw
Path('referee_checks.json').write_text(json.dumps(out, indent=2))
print(json.dumps({k:v for k,v in out.items() if k != 'saved_validity'}, indent=2))
