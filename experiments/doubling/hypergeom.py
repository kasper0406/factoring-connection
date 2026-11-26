import ctypes
import json
import math
import time
from flint import nmod_poly
from experiment import HERE,W,native
from structural import Poly,ptr,init,clear,evaluate

fit=native('nmod_poly_fit_length',None,ptr,ctypes.c_long)

def block(N,first,last):
    if last-first==1:
        j=first
        a=nmod_poly([(4*j-3)*(4*j-1),32*j-16,16],N)
        c=nmod_poly([18*j*j,36*j,18],N)
        return a,a,c
    mid=(first+last)//2
    a,b,c=block(N,first,mid)
    e,f,g=block(N,mid,last)
    return e*a,f*a+g*b,g*c

def evaluations(F,N,points):
    tmp=Poly()
    init(ctypes.byref(tmp),N)
    try:
        values=(W*len(points))()
        coeffs=(W*len(F))(*map(int,F))
        fit(ctypes.byref(tmp),len(F))
        ctypes.memmove(tmp.coeffs,coeffs,ctypes.sizeof(coeffs))
        tmp.length=len(F)
        evaluate(values,ctypes.byref(tmp),points,len(points))
        return values
    finally:
        clear(ctypes.byref(tmp))

def term(j):
    return (4*j-3)*(4*j-1),18*j*j

def truncation(N,K,plain=False):
    steps=K-1
    u=S=c=1
    offset=0
    if steps and not plain:
        width=math.isqrt(steps)
        count=steps//width
        points=(W*count)(*(j*width for j in range(count)))
        polys=block(N,1,1+width)
        A,B,C=(evaluations(F,N,points) for F in polys)
        for a,b,g in zip(A,B,C):
            S=(b*u+g*S)%N
            u=a*u%N
            c=g*c%N
        offset=count*width
    for j in range(offset+1,K):
        a,g=term(j)
        S=(a*u+g*S)%N
        u=a*u%N
        c=g*c%N
    factor=math.gcd(N,c)
    return (S*pow(c,-1,N)%N if factor==1 else None),factor

def universal(N):
    early=math.gcd(N,6)
    if early>1:
        return None,early,early%4==3
    K=math.isqrt(N)//2+1
    T,factor=truncation(N,K)
    return T,factor,(factor%4==3 if factor>1 else T==0)

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description='Evaluate the CM truncation T0 test on a distinct odd semiprime N=1 mod 4.')
    parser.add_argument('N',nargs='?',type=int)
    args=parser.parse_args()
    if args.N is not None:
        if args.N < 5 or args.N % 4 != 1 or args.N.bit_length() > 8*ctypes.sizeof(W):
            parser.error("N must be a positive native-word semiprime congruent to 1 modulo 4")
        T,g,inert=universal(args.N)
        print(json.dumps(dict(N=args.N,cm_sum=T,denominator_gcd=g,both_3_mod_4=inert)))
        raise SystemExit
    out=[]
    for size in (32,48,64):
        for line in (HERE/f'rows{size}.jsonl').read_text().splitlines():
            row=json.loads(line)
            if row['index']>=4:
                continue
            N=row['N']
            L=1<<(size//2-1)
            start=time.perf_counter()
            T,factor=truncation(N,L)
            elapsed=time.perf_counter()-start
            assert factor==1 and math.gcd(T,N)==(N if row['c']==3 else 1)
            for s,aa in ((row['p'],row['ap']),(row['q'],row['aq'])):
                expected=0 if aa is None else (-1)**(((s-1)//4)%2)*2*aa*pow(3,(s-1)//2,s)%s
                assert T%s==expected
            if size==32:
                assert truncation(N,L,plain=True)==(T,1)
            out.append(dict(bits=size,N=N,T=T,gcd=math.gcd(T,N),seconds=elapsed))
        print('Verified CM sums at',size,flush=True)
    z1=sum(math.comb(4*k,2*k)*math.comb(2*k,k)*pow(64,-k,21) for k in range(21))%21
    assert z1==1
    (HERE/'hypergeom.json').write_text(json.dumps(dict(rows=out,z1_N21=z1),indent=2)+'\n')
