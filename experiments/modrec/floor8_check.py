"""Factor-blind hull evaluator of the floor8 short-square query vs factor-backed lattice_sum."""
import sys, json, time, random
from math import isqrt, log
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "floor8")); sys.path.insert(0, str(HERE.parents[1]))
import os; os.chdir(HERE.parent / 'floor8')
from experiment import lattices, lattice_sum, prime, FORMS
import sympy as sp
from hull_count import short_square

def reference(p, q, D):
    n = p*q
    return lattice_sum(lattices(p, q, D), isqrt(n)+1, False)[0] % 8

rows = []
fixtures = [(33119086993,), (38942974393,)]
cases = []
for (n,) in fixtures:
    f = sp.factorint(n); p, q = sorted(f)
    cases.append(('fixture', p, q))
rng = random.Random(20261127)
for bits in (30, 34, 38, 42, 46, 50, 54):
    for _ in range(20 if bits <= 42 else (6 if bits <= 46 else 2)):
        while True:
            p, q = prime(rng, bits//2), prime(rng, bits//2)
            if p != q and (p*q).bit_length() == bits: break
        cases.append((bits, p, q))
agree = 0
for tag, p, q in cases:
    n = p*q
    for D, (b, c) in FORMS.items():
        st = {}
        t0 = time.time()
        s = short_square(n, b, c, 3, st)
        dt = time.time() - t0
        ref = reference(p, q, D)
        agree += (s == ref)
        rows.append(dict(tag=tag, n=n, D=D, S=s, ref=ref, vertices=st['vertices'], tests=st['tests'], sec=round(dt, 3)))
        print(tag, n.bit_length(), D, s, ref, st['vertices'], round(dt,2), flush=True)
print('agree', agree, 'of', len(rows))
Path(HERE / 'floor8_check.json').write_text(json.dumps(rows, indent=0))
# growth: vertices per N^{1/3}
import numpy as np
for D in (-23, -4):
    xs = np.array([log(r['n']) for r in rows if r['D'] == D]); ys = np.array([log(r['vertices']) for r in rows if r['D'] == D])
    print(D, 'slope log(vertices) vs log N:', round(float(np.polyfit(xs, ys, 1)[0]), 3))
for r in rows:
    if r['tag'] == 'fixture': print(r)
