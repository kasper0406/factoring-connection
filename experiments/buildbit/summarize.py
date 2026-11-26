import json
import math
from collections import Counter
from pathlib import Path

TAIL = 0.01 / 40  # two tails, 20 registered rates

def log_cdf(k,n,p):
    if p == 0: return 0
    if p == 1: return 0 if k == n else -math.inf
    terms = [math.lgamma(n+1)-math.lgamma(j+1)-math.lgamma(n-j+1)
             +j*math.log(p)+(n-j)*math.log1p(-p) for j in range(k+1)]
    peak = max(terms)
    return peak + math.log(sum(math.exp(t-peak) for t in terms))

def upper(k,n):
    if k == n: return 1.
    if k == 0: return -math.expm1(math.log(TAIL)/n)
    low,high = k/n,1.
    for i in range(44):
        m = (low+high)/2
        if log_cdf(k,n,m)>math.log(TAIL): low = m
        else: high = m
    return high

def rate(k,n):
    return [k/n, 1-upper(n-k,n) if k else 0., upper(k,n)]

def score(pos,neg,scale):
    return [(pos[0]-neg[0])*scale,(pos[1]-neg[2])*scale,(pos[2]-neg[1])*scale]

def load(topic,n):
    return [json.loads(x) for x in Path(f'{topic}{n}_out.jsonl').read_text().splitlines()]

def cf(n):
    rows = load('cf',n)
    result = {}
    for label in (0,1):
        all_rows = [r for r in rows if r['label']==label]
        clean = [r for r in all_rows if not r['splits']]
        single = [len(r['individual']) for r in clean]
        result[label] = dict(n=len(clean),splits=len(all_rows)-len(clean),full=sum(r['full'] for r in clean),
             first=sum(r['first'] for r in clean),union=sum(h>0 for h in single),individual_sum=sum(single),
             multiple=sum(h>1 for h in single),all_full=sum(r['full'] for r in all_rows),
             splits_via=dict(Counter(v for r in all_rows for v in r['splits'])),
             individual_hits=dict(Counter(k for r in clean for k in r['individual'])))
    result['advantage'] = [v/4 for v in rate(result[1]['full'],result[1]['n'])]
    return result

def cubic(n):
    rows = load('cubic',n)
    result={}
    for label in (0,1):
        cell=[r for r in rows if r['label']==label]
        clean=[r for r in cell if not r['splits'] and not r.get('degenerate')]
        result[label] = dict(n=len(clean),splits=sum(bool(r['splits']) for r in cell),
           degenerate=sum(bool(r.get('degenerate')) for r in cell),
           split_via=dict(Counter(k for r in cell for k in r['splits'])),
           zero=sum(r['zero'] for r in clean),plus=sum(r['sign']==1 for r in clean),
           oracle=sum(r['oracle_square'] for r in clean))
    for k in ('zero','plus'):
        pos,neg = [rate(result[l][k],result[l]['n']) for l in (1,0)]
        result[k]=score(pos,neg,0.5)
    return result

if __name__=='__main__':
    stats={'cf':{n:cf(n) for n in (40,52,64,80)}, 'cubic':{n:cubic(n) for n in (40,60,80)}}
    Path('summary.json').write_text(json.dumps(stats,indent=2)+'\n')
    for name,section in stats.items():
        for n,v in section.items():
            print(name,n,v)
