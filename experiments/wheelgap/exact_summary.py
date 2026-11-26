import csv
import json
import math
from collections import defaultdict
import numpy as np

groups=defaultdict(dict)
for line in open('trials.jsonl'):
    row=json.loads(line)
    groups[row['bits'],row['typ'],row['trial']][row['method']]=row
for line in open('exact.jsonl'):
    row=json.loads(line)
    groups[row['bits'],row['typ'],row['trial']]['exact']=row
assert len(groups)==384 and all(len(g)==4 for g in groups.values())
items=defaultdict(list)
for (n,g,j), rows in groups.items():
    exact=rows['exact']
    logW=math.log2(math.prod(rows['lll']['wheel']))
    minint=exact['h']
    assert 0<minint<=rows['lll']['h']
    items[n,g,'integer minimum'].append(math.log2(minint)/logW)
    items[n,g,'rational minimum'].append(math.log2(exact['a'])/logW)
    items[n,g,'integer ms'].append(exact['int_secs']*1000)
    items[n,g,'rational ms'].append(exact['rat_secs']*1000)
    items[n,g,'best registered integer'].append(math.log2(min(r['h'] for m,r in rows.items() if m!='exact' and r['h']))/logW)
    items[n,g,'best of XOR and LLL integer'].append(math.log2(min(r['h'] for m,r in rows.items() if m in ('xor','lll') and r['h']))/logW)

indices=np.random.default_rng(20261126).integers(0,32,(5000,32))
def ci(data):
    a=np.array(data)
    return (a.mean(),*np.percentile(a[indices].mean(axis=1),[2.5,97.5]))

with open('exact_summary.csv','w',newline='') as f:
    w=csv.writer(f);w.writerow(('bits','type','metric','mean','ci_lo','ci_hi'))
    for (n,g,metric),a in sorted(items.items()):
        w.writerow((n,g,metric,*ci(a)))
        if g=='P':
            contrast=np.array(a)-np.array(items[n,'R',metric])
            w.writerow((n,'P-R',metric,*ci(contrast)))
            print(n,metric,*(f'{x:.4f} [{l:.4f},{u:.4f}]' for x,l,u in (ci(a),ci(items[n,'R',metric]),ci(items[n,'S',metric]),ci(contrast))))
