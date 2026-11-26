import json
from collections import Counter
from math import gcd
from pathlib import Path
from cypari2 import Pari
from gmpy2 import is_prime
from sympy import factorint

def pencil(N):
    d = 4*N-1
    return -N*d, (N,N,(5*N-1)//4), 8*N-1

def check():
    pari = Pari()
    counts = Counter()
    for p in range(5,160,2):
        if not is_prime(p): continue
        for q in range(p+2,min(160,2*p),2):
            if not is_prime(q) or p%4 != q%4: continue
            N = p*q
            D,coeffs,m = pencil(N)
            f = pari.Qfb(*coeffs).qfbred()
            assert gcd(D,m)==1 and sum(t*u for t,u in zip(coeffs,(1,2,4)))==m
            identity = pari.Qfb(1,1,(1-D)//4)
            assert f.qfbpow(2) == identity
            data = pari.quadclassunit(D)
            orders, gens = data[1], data[2]
            # Enumerate the image of all reduced classes under squaring using
            # generators. This class-group computation is oracle validation.
            squares = [identity]
            for n,g in zip(orders,gens):
                gen=g.qfbpow(2)
                current=identity
                new=[]
                for i in range(int(n)//gcd(int(n),2)):
                    new.extend(y.qfbcomp(current) for y in squares)
                    current=current.qfbcomp(gen)
                squares=new
            assert (f in squares) == (p%4==1)
            kind='nonfundamental' if max(factorint(4*N-1).values())>1 else 'fundamental'
            counts[kind]+=1
    return dict(counts)

if __name__=='__main__':
    Path('genus_counts.json').write_text(json.dumps(check(),indent=2)+'\n')
