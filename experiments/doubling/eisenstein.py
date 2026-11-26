from functools import lru_cache
from math import gcd
import json
from flint import fmpq, fmpz
from structural import pure_e4
from experiment import HERE

@lru_cache(maxsize=None)
def local_form(k):
    a,b = pure_e4(k)
    constant = -fmpq.bernoulli(k)/(2*k)
    value = constant*fmpq(a,b)
    return int(value.p),int(value.q),int(constant.p),int(constant.q)

def at_prime(s,t):
    h=(t-1)%(s-1)
    if not h:
        # Hasse power is a unit at this ordinary curve.
        return 'unit'
    a,b,c,d=local_form(h)
    assert b%s and d%s
    if c%s:
        return 'unit' if a%s else 'zero'
    return 'unresolved' if a%s==0 else 'pole'

if __name__=='__main__':
    s,q=101,269
    w=s*q-1
    a,b,c,d=local_form(68)
    fullc=-fmpq.bernoulli(w)/(2*w)
    assert fmpz(q).is_prime() and (w%100)==68
    assert (a*pow(b,-1,s))%s and c%s==0 and int(fullc.p)%s==0 and int(fullc.q)%s
    assert int(fullc.p)%s**2!=0
    primes=[n for n in range(5,5001,4) if fmpz(n).is_prime()]
    candidates=[]
    for p in primes:
        if p>401:
            break
        for q in primes:
            if q>p and at_prime(q,p)=='zero':
                candidates.append([p,q,at_prime(p,q)])
    result=dict(pole=dict(N=w+1,p=s,q=269,weight=w,weight_mod_100=68,
                          G68_CM=(a*pow(b,-1,s))%s,
                          constant_mod_101_squared=int(fullc.p)*pow(int(fullc.q),-1,s*s)%(s*s)),
                larger_prime_zeros=candidates)
    (HERE/'eisenstein.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result)
