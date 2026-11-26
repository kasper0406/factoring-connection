import csv
import json
import math
from collections import defaultdict
import numpy as np

data = defaultdict(list)
for line in open('trials.jsonl'):
    r=json.loads(line)
    data[r['bits'],r['typ'],r['method']].append(r)
rng=np.random.default_rng(20261126)
boot=rng.integers(0,32,(5000,32))

def ci(a):
    a=np.array(a)
    ix=boot if len(a)==32 else rng.integers(0,len(a),(5000,len(a)))
    return [float(a.mean()),*map(float,np.percentile(a[ix].mean(axis=1),[2.5,97.5]))]

def alpha(h,r):
    return math.log2(h)/math.log2(math.prod(r['wheel']))

def metrics(rows, cap, method):
    rh=[r for r in rows if r['h']]
    integer=ci([alpha(r['h'],r) for r in rh]) if rh else [None]*3
    rat=ci([alpha(max(r['a'],r['b']),r) for r in rh]) if rh else [None]*3
    if method != 'lll':
        restricted=ci([alpha(min(r['prefix'],cap),r) for r in rows])
    else: restricted=[None]*3
    return dict(hits=len(rh),integer_conditional=integer,rational_conditional=rat,
                restricted=restricted,cap=cap,time_ms=ci([1000*r['secs'] for r in rows]),
                overlap=np.mean([r['overlap'] for r in rows]))
result={}
with open('summary.csv','w',newline='') as f:
    w=csv.writer(f)
    w.writerow(['bits','type','method','hits','cap','overlap','int_mean','int_low','int_high','rat_mean','rat_low','rat_high','restricted_mean','restricted_low','restricted_high','ms_mean','ms_low','ms_high'])
    for n in (42,82,120,168):
        for m in ('scan','lll','xor'):
            cap={'scan':2**22,'lll':0,'xor':2**19}[m]
            for g in 'PRS':
                for r in data[n,g,m]:
                    if not r['h'] and m!='lll': cap=min(cap,r['prefix'])
            for g in 'PRS':
                a=metrics(data[n,g,m],cap,m)
                result[n,g,m]=a
                w.writerow([n,g,m,a['hits'],cap,a['overlap'],*a['integer_conditional'],*a['rational_conditional'],*a['restricted'],*a['time_ms']])
                def fmt(x): return ' --- ' if x[0] is None else '%.4f [%.4f,%.4f]' % tuple(x)
                print(n,g,m,a['hits'],'int',fmt(a['integer_conditional']),'rat',fmt(a['rational_conditional']),'restricted',fmt(a['restricted']),'time',fmt(a['time_ms']), 'overlap',a['overlap'])
            if m=='xor':
                ds=[alpha(min(l['prefix'],cap),l)-alpha(min(r['prefix'],cap),r) for l,r in zip(data[n,'P',m],data[n,'R',m])]
                print('PRIMARY P-R',n,cap,ci(ds))
print('maximum cpu', max(r['secs'] for rows in data.values() for r in rows))
