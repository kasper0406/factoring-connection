import json
import math
from flint import fmpz
from hypergeom import universal
from experiment import HERE,exact

def main():
    small=[]
    primes=[s for s in range(3,500,2) if fmpz(s).is_prime()]
    for p in primes:
        for q in primes:
            N=p*q
            if p<q and N<1500 and p%4==q%4:
                small.append((N,p,q))
    big=[]
    for size in (32,48,64):
        for line in (HERE/f'rows{size}.jsonl').read_text().splitlines():
            r=json.loads(line)
            if r['index']==0:
                big.append((r['N'],r['p'],r['q']))
    data={}
    for name,rows in [('small',small),('big',big)]:
        splits={}
        b_splits={}
        for N,p,q in rows:
            T,g,inert=universal(N)
            assert inert==(p%4==3)
            if g>1:
                assert g==p
                splits[N]=g
            else:
                assert all((T%s==0)==inert for s in (p,q))
            b=math.gcd(exact(N)[0],N)
            if 1<b<N:
                b_splits[N]=b
        data[name]=dict(count=len(rows),denominator_splits=splits,binomial_splits=b_splits)
        print(name,len(rows),'denominator splits',len(splits),'B splits',len(b_splits))
    (HERE/'uniform.json').write_text(json.dumps(data,indent=2)+'\n')
if __name__=='__main__':
    main()
