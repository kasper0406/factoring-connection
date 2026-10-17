"""Exact Round-2 diagnostics. Candidate functions receive N and public bounds only.
Harness factor labels occur only in main(). No external packages, floats, or factoring library.
Run: python3 -B round2_experiments.py > round2_results.json
"""
from math import gcd, isqrt, comb
from itertools import product
from collections import Counter
from fractions import Fraction
import json
import random
import time
import tracemalloc


def primes(b):
    return [x for x in range(2, b+1) if all(x%d for d in range(2,isqrt(x)+1))]


def trim(a):
    while a and a[-1]==0: a.pop()
    return a


def divexact(a,b):
    a=a[:]; out=[0]*max(0,len(a)-len(b)+1)
    while len(a)>=len(b):
        i=len(a)-len(b); c=a[-1]//b[-1]; out[i]=c
        for j,v in enumerate(b): a[i+j]-=c*v
        trim(a)
    assert not a
    return trim(out)


def cyclotomic(m,cache={}):
    if m not in cache:
        a=[-1]+[0]*(m-1)+[1]
        for d in range(1,m):
            if m%d==0: a=divexact(a,cyclotomic(d))
        cache[m]=a
    return cache[m]


class Field:
    def __init__(self,p,mod):
        self.p=p; self.mod=mod; self.k=len(mod)-1; self.Q=p**self.k
    def digits(self,x):
        r=[]
        for _ in range(self.k): r.append(x%self.p); x//=self.p
        return r
    def encode(self,a): return sum((v%self.p)*self.p**i for i,v in enumerate(a))
    def add(self,x,y): return self.encode([a+b for a,b in zip(self.digits(x),self.digits(y))])
    def neg(self,x): return self.encode([-a for a in self.digits(x)])
    def mul(self,x,y):
        a=self.digits(x); b=self.digits(y); c=[0]*(2*self.k-1)
        for i,v in enumerate(a):
            for j,w in enumerate(b): c[i+j]+=v*w
        for i in range(len(c)-1,self.k-1,-1):
            for j in range(self.k): c[i-self.k+j]-=c[i]*self.mod[j]
        return self.encode(c[:self.k])
    def power(self,x,e):
        r=1
        while e:
            if e&1: r=self.mul(r,x)
            x=self.mul(x,x); e//=2
        return r
    def roots(self,q):
        return sum(self.add(self.add(self.power(self.add(x,1),q),self.neg(self.power(x,q))),self.neg(1))==0 for x in range(self.Q))
    def logs(self):
        for g in range(1,self.Q):
            logs={}; x=1
            for j in range(self.Q-1): logs[x]=j; x=self.mul(x,g)
            if len(logs)==self.Q-1: return logs
        raise AssertionError('no generator')


def jacobi_count(field,q):
    # Exact arithmetic in Z[z]/Phi_(Q-1)(z); z is a primitive root of unity.
    m=field.Q-1; logs=field.logs(); pairs=[]
    for x in range(1,field.Q):
        y=field.add(1,field.neg(x))
        if y: pairs.append((logs[x],logs[y]))
    total=[0]*m
    for a in range(m):
        for b in range(m):
            A=Counter((q*(a*x+b*y))%m for x,y in pairs)
            B=Counter((-a*x-b*y)%m for x,y in pairs)
            for i,v in A.items():
                for j,w in B.items(): total[(i+j)%m]+=v*w
    phi=cyclotomic(m)
    while len(total)>=len(phi):
        c=total[-1]; d=len(total)-len(phi)
        for j,v in enumerate(phi): total[d+j]-=c*v
        trim(total)
    expected=(field.roots(q)-2)*m*m
    assert trim(total[:])==trim([expected]), (field.Q,q,total,expected)
    return field.roots(q)


def p_rem(a,b,p):
    a=trim([x%p for x in a]); inv=pow(b[-1],-1,p)
    while len(a)>=len(b):
        d=len(a)-len(b); c=a[-1]*inv%p
        for j,v in enumerate(b): a[d+j]=(a[d+j]-c*v)%p
        trim(a)
    return a


def p_gcd(a,b,p):
    while b: a,b=b,p_rem(a,b,p)
    return a


def quotient_distribution(p,q,r):
    N=p*q; f=trim([0]+[comb(N,j)%p for j in range(1,N)]+[0])
    rank=Counter(); vanishes=0
    for cs in product(range(p),repeat=r):
        h=list(cs)+[1]; rem=p_rem(f,h,p)
        vanishes+=not rem
        rank[r-(len(p_gcd(h,rem,p))-1)]+=1
    return {'p':p,'q':q,'r':r,'zero_remainders':vanishes,'rank_counts':dict(sorted(rank.items())),'total':p**r}


def bounds(N):
    lo=isqrt(N//2)
    while 2*lo*lo<N: lo+=1
    return lo,isqrt(N)


def carry_frontier(N,kmax):
    L,U=bounds(N); QL=(N+U-1)//U; QU=N//L
    states=[(1,1)]; out=[]
    for k in range(1,kmax+1):
        if k>1:
            M=1<<k; states=[(a+t*(M//2),pow(a+t*(M//2),-1,M)*(N%M)%M) for a,b in states for t in (0,1)]
        M=1<<k; good=[]
        for a,b in states:
            pmin=L+(a-L)%M; pmax=U-(U-a)%M
            qmin=QL+(b-QL)%M; qmax=QU-(QU-b)%M
            if pmin<=pmax and qmin<=qmax and pmin*qmin<=N<=pmax*qmax:
                good.append((a,b))
        states=good
        out.append({'k':k,'bare_low_bit_states':1<<(k-1),'after_integer_box_lookahead':len(states)})
    return out


def masks(N,ps):
    return [(l,{s for s in range(l) if (s*s-4*N)%l in {x*x%l for x in range(l)}}) for l in ps if N%l]


def possible(A,B,l,S):
    return any(A+(s-A)%l<=B for s in S)


def interval_mask_tree(N,ps,stop=1):
    A=isqrt(4*N)
    if A*A<4*N: A+=1
    B=isqrt(9*N//2)
    ms=masks(N,ps); stack=[(A,B)]; visited=0; leaves=[]; forced=0
    while stack:
        a,b=stack.pop(); visited+=1
        if any(not possible(a,b,l,S) for l,S in ms): continue
        if b-a+1>=max(ps): forced+=1
        if b-a+1<=stop: leaves.append((a,b)); continue
        c=(a+b)//2; stack.extend([(a,c),(c+1,b)])
    direct=[s for s in range(A,B+1) if all(s%l in S for l,S in ms)]
    assert leaves==[] or sorted(a for a,b in leaves)==direct
    return {'N':N,'interval':[A,B],'width':B-A+1,'max_aux_prime':max(ps),'visited_nodes':visited,'nodes_automatically_passing_each_separate_mask':forced,'survivors':len(direct),'survivor_values':direct[:20]}


def weighted_collision(N,ps,amax=4):
    L,U=bounds(N); found=None; candidates=0; raw=0; gcds=0
    for a in range(1,amax+1):
        for b in range(1,amax+1):
            if gcd(a,b)!=1: continue
            D=isqrt(4*a*b*N)
            if D*D<4*a*b*N: D+=1
            # Convex function a*x+b*N/x attains maximum at interval endpoints.
            H=max((a*L*L+b*N)//L,(a*U*U+b*N)//U)
            ms=masks(a*b*N,ps)
            for s in range(D,H+1):
                raw+=1
                if any(s%l not in S for l,S in ms): continue
                candidates+=1; j=s-D
                for exponent in [a+b*N-D,a*N+b-D]:
                    g=gcd(pow(2,j,N)-pow(2,exponent,N),N); gcds+=1
                    if 1<g<N and found is None: found=g
    return {'N':N,'a_max':amax,'raw_integer_offsets_scanned':raw,'mask_survivors':candidates,'gcd_calls':gcds,'first_factor':found,'note':'scanning implementation; no fast compressed successor claim'}


def cf_rows(N,K,T,B,schedule=None):
    ps=primes(B); rows=[]; maxbits=0; checked=0; zeros=0
    for k,steps in (schedule if schedule is not None else [(i,T) for i in range(1,K+1)]):
        D=k*N; a0=isqrt(D)
        if a0*a0==D: continue
        m=0; d=1; a=a0; hm2,hm1=0,1; vm2,vm1=1,0
        for j in range(steps):
            u=a*hm1+hm2; v=a*vm1+vm2
            hm2,hm1=hm1,u; vm2,vm1=vm1,v
            Q=u*u-D*v*v; checked+=1
            maxbits=max(maxbits,u.bit_length(),v.bit_length(),abs(Q).bit_length())
            if not Q: zeros+=1
            else:
                rem=abs(Q); es=[int(Q<0)]
                for ell in ps:
                    e=0
                    while rem%ell==0: e+=1; rem//=ell
                    es.append(e)
                if rem==1: rows.append((u%N,es,k,j,Q))
            m=d*a-m; d=(D-m*m)//d; a=(a0+m)//d
    # Remove repetitions: primitive information is the modular square identity.
    unique={ (u,tuple(es)): (u,es,k,j,Q) for u,es,k,j,Q in rows }
    basis={}; dependencies=0; factor=None
    selected=list(unique.values())
    for i,(u,es,k,j,Q) in enumerate(selected):
        vec=sum((e%2)<<t for t,e in enumerate(es)); witness=1<<i
        while vec:
            b=vec.bit_length()-1
            if b not in basis: basis[b]=(vec,witness); break
            v0,w0=basis[b]; vec^=v0; witness^=w0
        if not vec:
            dependencies+=1; X=1; totals=[0]*(len(ps)+1)
            for h,(uh,eh,*_) in enumerate(selected):
                if witness>>h&1:
                    X=X*uh%N; totals=[x+y for x,y in zip(totals,eh)]
            assert all(e%2==0 for e in totals)
            Y=1
            for ell,e in zip(ps,totals[1:]): Y=Y*pow(ell,e//2,N)%N
            assert (X*X-Y*Y)%N==0
            for g in (gcd(X-Y,N),gcd(X+Y,N)):
                if 1<g<N and factor is None: factor=g
    return {'N':N,'K':K,'T':T,'B':B,'convergents':checked,'smooth_rows':len(rows),'unique_rows':len(unique),'parity_rank':len(basis),'dependencies':dependencies,'factor':factor,'max_integer_bits':maxbits,'schedule':schedule,'note':'fixed multiplier control, no adaptive-source advantage claimed'}


def adaptive_cf(N,K=8,pilot=4,total=256,B=31):
    scores=[]
    for k in range(1,K+1):
        r=cf_rows(N,1,pilot,B,[(k,pilot)])
        scores.append((r['unique_rows'],r['parity_rank'],-k))
    chosen=[-x[2] for x in sorted(scores,reverse=True)[:2]]
    extra=(total-K*pilot)//2
    schedule=[(k,pilot+(extra if k in chosen else 0)) for k in range(1,K+1)]
    out=cf_rows(N,K,pilot,B,schedule)
    out['pilot_convergents_extra']=K*pilot
    out['total_convergents_including_recomputed_pilots']=out['convergents']+K*pilot
    out['chosen_multipliers']=chosen
    out['note']='pilot smoothness/rank score; all pilot and recomputation costs counted'
    return out


def generated_controls(seed=104729):
    rng=random.Random(seed); result=[]
    def prime_candidate(bits):
        while True:
            x=rng.randrange(1<<(bits-1),1<<bits)|1
            if all(x%d for d in range(3,isqrt(x)+1,2)): return x
    for bits in [8,12,16,24,32]:
        for rep in range(4):
            p=prime_candidate(bits); q=prime_candidate(bits)
            while q==p: q=prime_candidate(bits)
            N=p*q
            row={'factor_bits':bits,'rep':rep,'N':N,'diagnostic_factors':sorted([p,q])}
            for name,fn in [('fixed',lambda:cf_rows(N,8,36,31)),('adaptive',lambda:adaptive_cf(N))]:
                tracemalloc.start(); start=time.perf_counter(); r=fn()
                elapsed=time.perf_counter()-start; _,peak=tracemalloc.get_traced_memory(); tracemalloc.stop()
                if r['factor'] is not None: assert 1<r['factor']<N and N%r['factor']==0
                r['seconds']=elapsed; r['peak_traced_bytes']=peak; row[name]=r
            result.append(row)
    return {'seed':seed,'generation':'independent uniform odd candidates of stated bit length, rejection by exact trial-division primality; distinct factors','B':31,'budget_per_method':288,'instances':result}


def shared_wheel_layout(N,ps,cap=16384):
    D=isqrt(4*N)
    if D*D<4*N: D+=1
    W=isqrt(9*N//2)-D+1
    states=[0]; M=1; rows=[]; construction=0; used=[]
    for ell,S in masks(N,ps):
        offsets=sorted((x-D)%ell for x in S)
        if len(states)*len(offsets)>cap: break
        oldM=M; inv=pow(oldM,-1,ell)
        states=[t+oldM*((a-t)*inv%ell) for t in states for a in offsets]
        M*=ell; construction+=len(states); used.append((ell,S))
        lengths=[(W-1-t)//M+1 for t in states if t<W]
        H=max(lengths,default=0); total=sum(lengths)
        baby=min(H,isqrt(total)+1) if H else 0
        giant=sum((h+baby-1)//baby for h in lengths) if baby else 0
        # Exact set-coverage check; diagnostic scans are not charged as a fast algorithm.
        represented={t+M*h for t in states if t<W for h in range((W-1-t)//M+1)}
        expected={j for j in range(W) if all((D+j)%l in A for l,A in used)}
        assert represented==expected
        rows.append({'M':M,'R_explicit':len(states),'R_nonempty':len(lengths),'offsets':total,'baby_entries':baby,'giant_entries':giant,'cumulative_CRT_states_created':construction,'layout_entries_plus_explicit_class_setup':baby+giant+len(states)})
    return {'N':N,'W':W,'status':'exact list layout and coverage only; fast cross-collision extraction not implemented','stages':rows}


def main():
    fields=[Field(p,[0,1]) for p in [3,5,7,11,17,31]]+[Field(3,[1,0,1]),Field(5,[3,0,1]),Field(2,[1,1,0,1])]
    count=0; jc=[]
    for F in fields:
        for q in [3,5,7,11,13,17,37,41]:
            if q==F.p: continue
            R=jacobi_count(F,q); jc.append([F.Q,q,R]); count+=1
    qs=[]
    for p,q in [(5,7),(7,11),(11,13)]:
        for r in [1,2]: qs.append(quotient_distribution(p,q,r)); qs.append(quotient_distribution(q,p,r))
    diagnostic=[(101,193),(1009,1481),(10007,19001)]
    ps=[3,5,7,11,13,17,19,23,29,31]
    out={'status':'exact diagnostics, no asymptotic claims','jacobi_exact_checks':count,'jacobi_counts':jc,'quotient_distributions':qs,'mask_tree':[],'carry_frontiers':[],'weighted_collisions':[],'cf_control':[],'shared_wheel_layout':[]}
    for p,q in diagnostic:
        N=p*q; row=interval_mask_tree(N,ps)
        assert row['interval'][0]<=p+q<=row['interval'][1]
        assert all((p+q)%ell in S for ell,S in masks(N,ps))
        out['mask_tree'].append(row)
        out['carry_frontiers'].append({'N':N,'rows':carry_frontier(N,min(10,N.bit_length()//4+1))})
        out['weighted_collisions'].append(weighted_collision(N,ps))
        out['shared_wheel_layout'].append(shared_wheel_layout(N,ps))
        out['cf_control'].append(cf_rows(N,8,32,31))
    # Exact scalar counterexample to gcd(q,Q-1)-only prediction.
    out['cf_holdout']=generated_controls()
    out['same_gcd_counterexample']={'field':17,'q19_roots':Field(17,[0,1]).roots(19),'q41_roots':Field(17,[0,1]).roots(41),'both_gcd':1}
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
