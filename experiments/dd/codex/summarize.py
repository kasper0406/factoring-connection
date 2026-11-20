"""Certified endpoint order statistics, fixed decimal input from measure.py."""
import csv
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
SCALE=10**9

def summarize(name,group,rows):
    out=dict(group=name,value=group,embeddings=len(rows))
    for quantity in ('loss','ratio'):
        lower=sorted(int(v[quantity+'_lo_1e9']) for v in rows)
        upper=sorted(int(v[quantity+'_hi_1e9']) for v in rows)
        for stat,idx in [('min',0),('median',(len(rows)-1)//2),('p90',(9*len(rows)+9)//10-1),('max',-1)]:
            out[f'{quantity}_{stat}_lo_1e9']=lower[idx]
            out[f'{quantity}_{stat}_hi_1e9']=upper[idx]
    for q in (25,50,100):
        cutoff=q*SCALE//100
        exact=sum(int(v['ratio_lo_1e9'])>=cutoff for v in rows)
        assert exact==sum(int(v['ratio_hi_1e9'])>=cutoff for v in rows)
        out[f'ge_{q}pct']=exact
    return out

def main():
    with (HERE/'embeddings.csv').open() as f:
        rows=list(csv.DictReader(f))
    fx=json.loads((HERE/'fixtures.json').read_text())
    cells=[json.loads(l) for l in (HERE/'cells.jsonl').read_text().splitlines()]
    grid={(int(bits),m,p) for bits in fx for m in range(3,9) for p in range(m,2*m+1)}
    assert grid=={(c['bits'],c['m'],c['p']) for c in cells} and len(grid)==len(cells)==195
    ids=[(int(r['bits']),int(r['m']),int(r['p']),r['side'],int(r['root'])) for r in rows]
    assert len(rows)==len(set(ids)) and set(ids)=={(c['bits'],c['m'],c['p'],side,i) for c in cells for side in 'PQ' for i in range(2*c['n'])}
    assert all(0<=int(r[q+'_hi_1e9'])-int(r[q+'_lo_1e9'])<1000 for r in rows for q in ('loss','ratio'))
    sums=[summarize('all','all',rows)]
    for col in ['m','bits','p','kind','side']:
        for val in sorted({r[col] for r in rows},key=lambda v:int(v) if v.isdigit() else v):
            sums.append(summarize(col,val,[r for r in rows if r[col]==val]))
    with (HERE/'summary.csv').open('w') as f:
        w=csv.DictWriter(f,list(sums[0]))
        w.writeheader()
        w.writerows(sums)
    with (HERE/'cell_summary.csv').open('w') as f:
        w=csv.DictWriter(f,['bits','m','p','n',*list(sums[0])[2:]])
        w.writeheader()
        for c in cells:
            group=[r for r in rows if all(int(r[k])==c[k] for k in ('bits','m','p'))]
            stats=summarize('cell',None,group)
            record={k:c[k] for k in ('bits','m','p','n')}
            record.update({k:stats[k] for k in list(stats)[2:]})
            w.writerow(record)
    for s in sums:
        print(s['group'],s['value'],s['embeddings'],*(f"{s[k+'_lo_1e9']/SCALE:.5f}" for k in ('loss_median','loss_max','ratio_median','ratio_p90','ratio_max')),s['ge_25pct'],s['ge_50pct'],s['ge_100pct'])
    print(max(rows,key=lambda r:int(r['ratio_lo_1e9'])))
    print(max(int(r[q+'_hi_1e9'])-int(r[q+'_lo_1e9']) for r in rows for q in ('loss','ratio')))
if __name__=='__main__':
    main()
