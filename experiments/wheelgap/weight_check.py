import csv
import json
import math
from collections import defaultdict
import numpy as np
from benchmark import SEED, lattice, signature

rows=defaultdict(list)
with open('sqrt_weights.jsonl','w') as out:
    for line in open('trials.jsonl'):
        r=json.loads(line)
        if r['method']!='lll': continue
        k={42:11,82:18,120:24,168:31}[r['bits']]
        result=lattice(r['wheel'],r['t'],SEED+1000*k+r['trial'],sqrt_weights=True)
        assert result['h'] and signature(result['h'],r['wheel'])==r['t']
        assert signature(result['a']*result['b'],r['wheel'])==r['t']
        out.write(json.dumps(dict(bits=r['bits'],typ=r['typ'],trial=r['trial'],**result))+'\n')
        rows[r['bits'],r['typ']].append([math.log2(result['h'])/math.log2(math.prod(r['wheel'])),1000*result['secs']])
rng=np.random.default_rng(SEED)
ix=rng.integers(0,32,(5000,32))
with open('sqrt_weights.csv','w',newline='') as f:
    w=csv.writer(f)
    w.writerow(['bits','type','int_mean','lo','hi','ms_mean','lo','hi'])
    for key,vals in rows.items():
        a=np.array(vals)
        mean=a.mean(axis=0)
        lo,hi=np.percentile(a[ix].mean(axis=1),[2.5,97.5],axis=0)
        w.writerow([*key,*np.stack([mean,lo,hi],axis=1).flat])
        if key[0]==168: print(key,mean,lo,hi)
