from experiment import HERE, fmpz
import random
import gmpy2
import math
import time
import json

def tree(values):
    layers = [list(map(gmpy2.mpz,values))]
    while len(layers[-1])>1:
        layer=layers[-1]
        layers.append([layer[i]*layer[i+1] if i+1<len(layer) else layer[i] for i in range(0,len(layer),2)])
    return layers

def remainders(P,moduli):
    levels=tree(moduli)
    rs=[P%levels[-1][0]]
    for layer in reversed(levels[:-1]):
        rs=[rs[i//2]%m for i,m in enumerate(layer)]
    return rs

if __name__=='__main__':
    L=1<<15
    rows=[json.loads(line) for line in (HERE/'rows32.jsonl').read_text().splitlines()]
    labels={r['N']:r['c'] for r in rows}
    pools={c:[s for s in range(L+c,2*L,4) if fmpz(s).is_prime()] for c in (1,3)}
    rng=random.Random(2026112691)
    while len(labels)<L:
        c=rng.choice([1,3])
        p,q=rng.sample(pools[c],2)
        if (p*q).bit_length()==32:
            labels[p*q]=c
    start=time.perf_counter()
    P=tree(range(L+3,2*L,4))[-1][0]
    residues=remainders(P,labels)
    gs=[math.gcd(N,int(z)) for N,z in zip(labels,residues)]
    elapsed=time.perf_counter()-start
    exact=sum(g==(N if c==3 else 1) for (N,c),g in zip(labels.items(),gs))
    splits=sum(1<g<N for N,g in zip(labels,gs))
    result=dict(K=L,product_bits=P.bit_length(),correct=exact,splits=splits,seconds=elapsed,
                by_class={c:list(labels.values()).count(c) for c in (1,3)})
    assert exact==L and splits==0
    (HERE/'batch.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result)
