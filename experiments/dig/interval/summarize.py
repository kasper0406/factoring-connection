"""Independently check recorded labels and witnesses and summarize the locked run."""
import csv
from collections import defaultdict
from fractions import Fraction
from math import comb, isqrt, sqrt
from pathlib import Path
import sys

from sympy import Rational, continued_fraction_convergents, continued_fraction_iterator, isprime


def tail(k, n):
    return sum(comb(n, i) for i in range(k, n+1)) / (1 << n)


def wilson(k, n):
    z = 1.959964
    p = k/n
    center = p + z*z/(2*n)
    radius = z * sqrt(p*(1-p)/n + z*z/(4*n*n))
    return tuple(round((center+s*radius)/(1+z*z/n), 4) for s in (-1,1))


def check(row):
    N, p, q, lo, X, truth, W = (row[k] for k in ('N','p','q','lo','X','truth','W'))
    assert p*q == N and p < q < 2*p and isprime(p) and isprime(q)
    assert row['bits'] == N.bit_length() and (p-lo)//X == truth in (0,1)
    assert lo+2*X <= isqrt(N) + 1
    alpha = Fraction(row['exponent'])
    assert X**alpha.denominator <= N**alpha.numerator < (X+1)**alpha.denominator
    start = lo + truth*X
    block = start + (X-W)//2
    assert block == row['block'+str(truth)] and W >= 1
    assert row['copp'] == (truth if block <= p < block+W else -1)
    c = start + X//2
    weights = {(1, 1)}
    for r in continued_fraction_convergents(continued_fraction_iterator(Rational(N, c*c))):
        u, v = int(r.p), int(r.q)
        if max(u, v) > 2*N.bit_length():
            break
        if u*v:
            weights.add((u,v))
    hit = any(abs((u*p+v*q) - (u*c+v*N//c)) <= 4 for u,v in weights)
    assert row['tangent'] == (truth if hit else -1)
    if row['soft0'] != row['soft1']:
        s = [Fraction(row[k]) if row[k] else float("inf") for k in ('soft0','soft1')]
        assert row['soft_pred'] == int(s[1] < s[0])
    assert row['pred'] == (truth if hit or row['copp'] >= 0 else row['soft_pred'])


def summarize(rows):
    correct = lambda data, key: sum(r[key] == r['truth'] for r in data)
    residual = [r for r in rows if r['copp'] == r['tangent'] == -1]
    soft = correct(residual, 'soft_pred')
    discordant = [r for r in residual if r['support_pred'] != r['soft_pred']]
    wins = correct(discordant, 'soft_pred')
    return dict(trials=len(rows), combined=correct(rows,'pred'),
                copp=sum(r['copp']>=0 for r in rows),
                tangent=sum(r['tangent']>=0 for r in rows),
                both=sum(min(r['tangent'],r['copp'])>=0 for r in rows),
                miss=len(residual), soft=soft, ci=wilson(soft,len(residual)),
                p=round(tail(soft,len(residual)),5), paired=round(tail(wins,len(discordant)),5),
                support=correct(rows,'support_pred'), support_miss=correct(residual,'support_pred'),
                coin=correct(rows,'coin'), left=sum(r['truth']==0 for r in rows),
                soft_ties=sum(r['soft0']==r['soft1'] for r in residual),
                expected_copp=round(sum(r['W']/r['X'] for r in rows),1))


groups = defaultdict(list)
holdout = sys.argv[1:] == ['--holdout']
with Path(__file__).with_name('holdout.csv' if holdout else 'trials.csv').open() as f:
    for raw in csv.DictReader(f):
        row = {k:v if k in ('exponent','soft0','soft1') else int(v) for k,v in raw.items()}
        check(row)
        for key in ((row['bits'],row['exponent']), ('all',row['exponent'])):
            groups[key].append(row)
assert len(groups) == (2 if holdout else 12)
for key, rows in groups.items():
    assert len(rows) == (3072 if holdout else 1152 if key[0] == 'all' else 384)
    print(key, summarize(rows))
