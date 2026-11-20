"""Verify reported counts and export per-box measurements from the exact records."""
import csv
import json
from collections import Counter
from pathlib import Path

root = Path(__file__).resolve().parent
sources = ['results.jsonl','lll_results.jsonl','sparse_results.jsonl'] + [
    f'diagonal_{n}.jsonl' for n in (60,72,84,96)]
counts = Counter()
rows = []
for name in sources:
    for raw in (root/name).read_text().splitlines():
        r = json.loads(raw)
        N,n,m = int(r['N']),r['n'],r['m']
        assert N.bit_length()==r['bits']
        assert int(r['span']['order_index'])%N==0
        assert r['span']['degrees'][-1]==2*n
        assert len(r['span']['degrees'])<=2*((2*n+m-1)//m)+1
        if r['mode']!='power':
            assert min(r['span']['square_dimensions'])>2*m-1
        caps = [r[k] for k in ('aggregate','gp','all')]
        for k in ('points','directions','smooth_directions','rank','new'):
            assert caps[0][k]<=caps[1][k]<=caps[2][k]
        assert all(c['new']<=c['new_within_cap']<=c['rank']<=c['smooth_directions']<=c['directions']<=c['points'] for c in caps)
        assert caps[-1]['points']+sum(r['nonunits'].values())==3**(2*m)-1
        for caps_label,c in zip(('aggregate','gp','all'),caps):
            rows.append(dict(N=N,bits=r['bits'],seed=r['seed'],fixture=r.get('kind','special'),
                             n=n,m=m,method=r['mode'],B=r['B'],cap=caps_label,
                             low_points=c['points'],directions=c['directions'],smooth=c['smooth_directions'],
                             parity_rank=c['rank'],incremental_rank=c['new'],mixed_span=c['mixed_k_span'],
                             finite_yield_exponent=c['eta']))
        counts[r['mode']] += 1
assert sum(counts.values())==128
with (root/'yield.csv').open('w',newline='') as f:
    writer = csv.DictWriter(f,rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
print(dict(counts))
