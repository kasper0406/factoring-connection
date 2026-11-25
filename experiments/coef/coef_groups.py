"""Exact enumeration of candidate torsion images. Labels I, A, D throughout are counts."""
from itertools import product
from collections import Counter
from cypari2 import Pari

def groups(l):
    matrices=[g for g in product(range(l),repeat=4) if (g[0]*g[3]-g[1]*g[2])%l]
    ident=(1,0,0,1)
    def mul(g,h):
        a,b,c,d=g;e,f,j,k=h
        return ((a*e+b*j)%l,(a*f+b*k)%l,(c*e+d*j)%l,(c*f+d*k)%l)
    def close(seed):
        found={ident}; queue=[ident]
        while queue:
            g=queue.pop()
            for t in seed:
                h=mul(g,t)
                if h not in found:
                    queue.append(h);found.add(h)
        return frozenset(found)
    identity_group=frozenset([ident])
    bygroup={identity_group:()}; queue=[identity_group]
    while queue:
        H=queue.pop()
        for t in set(matrices)-H:
            seeds=(*bygroup[H],t);J=close(seeds)
            if J not in bygroup:
                bygroup[J]=seeds;queue.append(J)
    inv={g:next(h for h in matrices if mul(g,h)==ident) for g in matrices}
    def conj(g,x):
        return mul(mul(g,x),inv[g])
    def det(g):
        a,b,c,d=g;return (a*d-b*c)%l
    def h(g):
        return (1+det(g)-g[0]-g[3])%l
    seen=set();stats=Counter()
    for H in bygroup:
        if H in seen or {det(g) for g in H}!=set(range(1,l)):
            continue
        seen.update(frozenset(conj(g,x) for x in H) for g in matrices)
        comm=close({mul(conj(g,x),inv[x]) for g in H for x in H})
        tag = ('zero' if all(h(g)==0 for g in H) else 'local abelian' if
               all(h(g)==h(mul(g,c)) for g in H for c in comm) else 'nonabelian')
        stats[tag]+=1
    print(l,dict(stats))
def witnesses():
    pari=Pari();E=pari.ellinit([0,0,0,-1,1])
    f=pari('3*x^4-6*x^2+12*x-1')
    seen=set()
    for r in range(5,100):
        if r==23 or not pari.isprime(r):
            continue
        factors=pari.factormod(f,r)[0]
        shape=tuple(sorted(int(pari.poldegree(g)) for g in factors))
        t=int(pari.ellap(E,r))
        for tag,present in [('3-cycle',shape==(1,3)),('4-cycle',shape==(4,)),
                            ('minus-I via square', r%3==1 and t%3==0)]:
            if present and tag not in seen:
                print(tag,r,shape,t);seen.add(tag)
if __name__=='__main__':
    groups(2);groups(3);witnesses()
