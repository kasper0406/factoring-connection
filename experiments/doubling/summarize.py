import json
import math
from functools import lru_cache
from experiment import HERE

@lru_cache(maxsize=None)
def fisher(n,A,B,x):
    def weight(k):
        return math.comb(B,k)*math.comb(n-B,A-k)
    cutoff = weight(x)
    numer = sum(w for k in range(max(0,A+B-n),min(A,B)+1) if (w:=weight(k))<=cutoff)
    return numer/math.comb(n,A)

def association(xs,ys):
    n=len(xs)
    A,B=sum(xs),sum(ys)
    both=sum(x and y for x,y in zip(xs,ys))
    agree=n-A-B+2*both
    sign = (n*both>A*B)-(n*both<A*B)
    return dict(n=n,positive=A,target_positive=B,both=both,agree=agree,
                sign=sign,p=fisher(n,A,B,both))

def passed(banks):
    keys = list(banks[0])
    cutoff=0.01/(3*len(keys))
    return [k for k in keys if max(b[k]['p'] for b in banks)<cutoff
            and len({b[k]['sign'] for b in banks})==1]

def exacts(banks):
    return [k for k in banks[0] if all(b[k]['agree']==b[k]['n'] for b in banks)
            or all(b[k]['agree']==0 for b in banks)]

def scan():
    data = {'pre':[],'post':[],'hidden':[], 'splits':{}, 'jacobi':{}}
    for size in (32,48,64):
        rows=[json.loads(s) for s in (HERE/f'rows{size}.jsonl').read_text().splitlines()]
        assert len(rows)==256
        split = [r for r in rows if r['c']==1]
        labels = [r['c']==3 for r in rows]
        for name in ('pre','post'):
            data[name].append({k:association([r[name][k] for r in rows],labels) for k in rows[0][name]})
        hidden = dict(sign_p=[r['ap']<0 for r in split],sign_q=[r['aq']<0 for r in split],
                      sign_product=[r['ap']*r['aq']<0 for r in split],
                      square_p=[2*r['ap']**2>r['p'] for r in split],
                      square_q=[2*r['aq']**2>r['q'] for r in split])
        data['hidden'].append({f'{k}/{t}':association([r['post'][k] for r in split],ys)
                               for k in split[0]['post'] if k.startswith(('mod','jac'))
                               for t,ys in hidden.items()})
        data['splits'][size] = dict(input_count=sum(bool(r['splits']) for r in rows),
                                    features={k:sum(k in r['splits'] for r in rows)
                                              for k in sorted(set().union(*(r['splits'] for r in rows)))})
        data['jacobi'][size] = sum(r['post']['jac1'] for r in split)
    for name in ('pre','post','hidden'):
        data[name+'_pass'] = passed(data[name])
        if name != 'hidden':
            data[name+'_exact'] = exacts(data[name])
        print(name,len(data[name][0]),'pass',data[name+'_pass'],'exact',data.get(name+'_exact'))
    print('split counts',data['splits'],'positive Jacobi',data['jacobi'])
    h=data['hidden']
    print('minimum single hidden p:',[min(x['p'] for x in b.values()) for b in h])
    print('best replicated hidden p:',min(max(b[k]['p'] for b in h) for k in h[0] if len({b[k]['sign'] for b in h})==1))
    (HERE/'counts.json').write_text(json.dumps(data,indent=2)+'\n')

if __name__ == '__main__':
    scan()
