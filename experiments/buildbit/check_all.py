import json
from collections import Counter
from math import comb, gcd, prod
from pathlib import Path
from cypari2 import Pari
from cubic_carry import carries
from cf import columns
from sympy import factorint

def validate():
    pari=Pari()
    known=set()
    for path in sorted(Path('.').glob('*.tsv')):
        for line in path.read_text().splitlines():
            N,p,q,label=map(int,line.split())
            assert N==p*q and p!=q and max(p,q)<2*min(p,q)
            known.update([p,q])
    for n in (40,60,80):
        for line in Path(f'lift{n}.jsonl').read_text().splitlines():
            row=json.loads(line); known.update((row['p'],row['q']))
    assert all(pari.isprime(s) for s in known)
    # Independent prime factorizations to verify selected right kernels
    # exhaustively, including composite coprime atoms with square factors.
    groups=[[12,75,98,8,27,18,50],[101*103**2,101*107,103**4,107,101],
            [103**3*107,103*107**3,103**2*107**2,109**2]]
    checks=0
    for group in groups:
        row=columns(group); fac=[factorint(x) for x in group]
        for mask in range(1<<len(group)):
            direct=Counter(); reduced=0
            for j in range(len(group)):
                if mask>>j&1: direct.update(fac[j]); reduced^=row[j]
            assert (not reduced)==all(v%2==0 for v in direct.values())
            checks+=1
    result=dict(distinct_primes=len(known),kernel_checks=checks)
    for n in (40,60,80):
        count=Counter()
        for line in Path(f'cubic{n}.tsv').read_text().splitlines():
            N,p,q,label=map(int,line.split())
            g=prod(s for s in (p,q) if carries(N//7,s))
            state='unit' if g==1 else ('zero' if g==N else 'proper')
            assert state==('proper' if label else 'unit' if p%7==1 else 'zero')
            count[state]+=1
        result[n]=dict(count)
    direct=0
    for p in range(3,160,2):
        if not pari.isprime(p): continue
        for q in range(p+2,min(160,2*p),2):
            N=p*q
            if not pari.isprime(q) or N%7!=1: continue
            g=gcd(comb(2*(N//7),N//7),N)
            assert (1<g<N)==(p%7 not in (1,6))
            direct+=1
    result['direct_conductor7']=direct
    return result

if __name__=='__main__':
    Path('verification.json').write_text(json.dumps(validate(),indent=2)+'\n')
