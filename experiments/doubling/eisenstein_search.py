from flint import fmpq, fmpz
from functools import lru_cache
from math import factorial,gcd
from experiment import HERE
from eisenstein import local_form
import json
import gmpy2

def hurwitz(t):
    # For g2=4, g3=0 the coefficient of z^(4j-2) in the Weierstrass function.
    cs=[0,fmpq(1,5)]
    G={4:fmpq(1,240)}
    for j in range(2,t//4+1):
        k=4*j
        num=6*sum(cs[i]*cs[j-i] for i in range(1,j))
        cs.append(num/((k-2)*(k-3)-12))
        G[k]=cs[-1]*factorial(k-2)/(2*48**j)
    return G

@lru_cache(maxsize=None)
def bern_num(h):
    return int((-fmpq.bernoulli(h)/(2*h)).p)

if __name__=='__main__':
    G=hurwitz(1200)
    for k in range(4,401,4):
        a,b,_,_=local_form(k)
        assert G[k]==fmpq(a,b)
    ps=[i for i in range(5,200001,4) if fmpz(i).is_prime()]
    hits=[]
    eligible=[]
    all_sieved=0
    for p in ps:
        if p>1201:
            break
        hs={h for h in range(4,p-1,4) if gcd(h+1,p-1)==1 and int(G[h].p)%p==0 and bern_num(h)%p}
        if not hs:
            continue
        eligible.append([p,sorted(hs)])
        rem=gmpy2.mpz(int(G[p-1].p))
        for q in ps:
            if q>p and (q-1)%(p-1) in hs and rem%q==0 and bern_num(p-1)%q:
                hits.append([p,q])
        for l in range(2,10001):
            if fmpz(l).is_prime():
                while rem%l==0:
                    rem//=l
        if rem>p and rem%4==1 and int((rem-1)%(p-1)) in hs and gmpy2.is_probab_prime(rem) and fmpz(int(rem)).is_prime():
            if bern_num(p-1)%rem:
                hits.append([p,int(rem)])
        all_sieved+=1
    out=dict(eligible=eligible, hits=hits, fully_sieved=all_sieved)
    (HERE/'eisenstein_search.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out)
