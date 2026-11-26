from math import gcd
from experiment import HERE,fmpz
import json

def full_sum(N):
    v=pow(72,-1,N)
    u=1
    coeff=1
    total=1
    for j in range(1,N):
        coeff,rem=divmod(coeff*4*(4*j-3)*(4*j-1),j*j)
        assert rem==0
        u=u*v%N
        total=(total+u*coeff)%N
    return total

if __name__=='__main__':
    ps=[s for s in range(5,301,2) if fmpz(s).is_prime()]
    out=[]
    for p in ps:
        for q in ps:
            N=p*q
            if p<q<4*p and p%4==q%4 and N<1500:
                T=full_sum(N)
                g=gcd(T,N)
                assert (T==0)==(p%4==3)
                out.append(dict(N=N,T=T,g=g))
    (HERE/'full_truncation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(len(out),'full-length sums; proper gcds',[(r['N'],r['g']) for r in out if 1<r['g']<r['N']])
