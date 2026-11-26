import json
import math
import numpy as np
from experiment import HERE

def standardize(X):
    X=X-X.mean(axis=0)
    den=np.sqrt(np.sum(X**2,axis=0))
    return X/np.where(den>0,den,1)

def rank(v):
    ordered=sorted(v)
    return [(ordered.index(x)+len(v)-1-ordered[::-1].index(x))/2 for x in v]

if __name__=='__main__':
    out=[]
    for size in (32,48,64):
        rows=[json.loads(l) for l in (HERE/f'rows{size}.jsonl').read_text().splitlines()]
        rows=[r for r in rows if r['c']==1]
        keys=[k for k in rows[0]['post'] if k.startswith(('mod','jac'))]
        X=standardize(np.array([[r['post'][k] for k in keys] for r in rows],dtype=float))
        cols=[[r[col] for r in rows] for col in ('ap','aq')]
        cols += [[abs(r[col])/math.sqrt(r[s]) for r in rows] for col,s in [('ap','p'),('aq','q')]]
        Y=standardize(np.array([rank(c) for c in cols]).T)
        corr=X.T@Y
        i,j=np.unravel_index(np.argmax(np.abs(corr)),corr.shape)
        observed=abs(corr[i,j])
        rng=np.random.default_rng(202611260+size)
        exceed=0
        for k in range(8191):
            exceed+=np.max(np.abs(X.T@Y[rng.permutation(len(Y))])) >= observed-1e-12
        out.append(dict(bits=size,feature=keys[i],target=['a_p','a_q','rel_p','rel_q'][j],
                        corr=corr[i,j],exceed=int(exceed),family_p=(int(exceed)+1)/8192))
    (HERE/'ranks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out)
