"""Phase-1 exact diagnostics. Candidate probes receive N and public parameters only."""
from math import comb, gcd, isqrt, prod
from random import Random
from fractions import Fraction
import json, time, resource
from pathlib import Path

MET = {"ring_multiplications":0, "gcd_calls":0, "inversions":0,
       "largest_stored_coefficient_bits":0, "largest_quotient_degree":0}
def gd(a,b):
    MET["gcd_calls"] += 1
    return gcd(a,b)

def mul(a,b,N,r,c):
    out=[0]*r
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            MET["ring_multiplications"]+=1
            v=x*y % N
            if i+j>=r:
                v=v*c%N
                MET["ring_multiplications"]+=1
            out[(i+j)%r]=(out[(i+j)%r]+v)%N
    return out

def power(a,k,N,r,c):
    out=[1]+[0]*(r-1)
    while k:
        if k&1: out=mul(out,a,N,r,c)
        k//=2
        if k: a=mul(a,a,N,r,c)
    return out

def sketch(N,r,c):
    MET["largest_quotient_degree"]=max(MET["largest_quotient_degree"],r)
    MET["largest_stored_coefficient_bits"]=max(MET["largest_stored_coefficient_bits"],N.bit_length())
    f=power([1,1]+[0]*(r-2),N,N,r,c)
    f[0]=(f[0]-1)%N
    f[N%r]=(f[N%r]-pow(c,N//r,N))%N
    return f

def matrix(f,N,r,c):
    return [[(f[(i-j)%r]*(c if i<j else 1))%N
             for j in range(r)] for i in range(r)]

def unit_probe(N,A):
    """Full pivot search. No factorization, reference prime, or field library."""
    A=[row[:] for row in A]
    r=len(A)
    for k in range(r):
        pivot=None
        for i in range(k,r):
            for j in range(k,r):
                d=gd(A[i][j],N)
                if 1<d<N: return {"status":"factor","factor":d,"pivots":k}
                if d==1 and pivot is None: pivot=(i,j)
        if pivot is None: return {"status":"singular","pivots":k}
        i,j=pivot
        A[k],A[i]=A[i],A[k]
        for row in A: row[k],row[j]=row[j],row[k]
        inv=pow(A[k][k],-1,N); MET["inversions"]+=1
        for i in range(k+1,r):
            scale=A[i][k]*inv%N; MET["ring_multiplications"]+=1
            for j in range(k+1,r):
                A[i][j]=(A[i][j]-scale*A[k][j])%N
                MET["ring_multiplications"]+=1
            A[i][k]=0
    return {"status":"full","pivots":r}

def universal_mask(N,r,c):
    # W=X(X+1)(X^2+X+1) = X^4+2X^3+2X^2+X.
    a=[0]*r
    for k,v in [(1,1),(2,2),(3,2),(4,1)]:
        a[k%r]=(a[k%r]+v*pow(c,k//r,N))%N
    return a

def coefficient_probe(N,f):
    for a in f:
        d=gd(a,N)
        if 1<d<N:return d
    return None

def candidate(N,r,c):
    """Operational repaired candidate: only N,r,c; universal masks are baseline."""
    before=dict(MET)
    for v in (6,r,c):
        d=gd(v,N)
        if 1<d<N:return {"stage":"parameter","factor":d}
    w=unit_probe(N,matrix(universal_mask(N,r,c),N,r,c))
    if w["status"]!="full":
        return {"stage":"universal","result":w}
    f=sketch(N,r,c)
    coeff=coefficient_probe(N,f)
    ans=unit_probe(N,matrix(f,N,r,c))
    return {"stage":"power","coeff_factor":coeff,"result":ans,
            "f":f,"metrics":{k:MET[k]-before[k] for k in ("ring_multiplications","gcd_calls","inversions")}}

def reference_prime(a):
    return a>1 and all(a%d for d in range(2,isqrt(a)+1))
def random_prime(bits,rng):
    while True:
        a=rng.randrange(1<<(bits-1),1<<bits)|1
        if reference_prime(a):return a
def field_rank(A,p):
    A=[[a%p for a in row] for row in A]
    r=len(A);k=0
    for j in range(r):
        i=next((i for i in range(k,r) if A[i][j]),None)
        if i is None:continue
        A[k],A[i]=A[i],A[k];inv=pow(A[k][j],-1,p)
        for i in range(k+1,r):
            t=A[i][j]*inv%p
            for v in range(j,r):A[i][v]=(A[i][v]-t*A[k][v])%p
        k+=1
    return k

def seed_checks():
    count=0
    for N in range(2,1001):
        p=next((d for d in range(2,isqrt(N)+1) if N%d==0),N)
        a=[comb(N,k)%N for k in range(1,N)]
        if p==N:
            assert not any(a)
            continue
        count+=1
        assert not any(a[:p-1]) and a[p-1]==N//p
        assert all(k*v%N==0 and (not v or 1<gcd(k,N)<N) for k,v in enumerate(a,1))
        factorial=1
        for K in range(1,isqrt(N)+2):
            assert (not any(a[:K-1]))==(K<=p)==(gcd(factorial,N)==1)
            factorial=factorial*K%N
    assert all((pow(a+1,561,561)-1-pow(a,561,561))%561==0 for a in range(561))
    return {"composites_through_1000":count,"Carmichael_scalar_checks":561}

def section_checks():
    rows=[]
    for M in (2,4,8,16):
        c=2;r=M//2
        den=[comb(M,i)*(-1)**i for i in range(M+1)]
        den[M]-=c
        num=[0]*M
        for j in range(M-r):
            num[r+j]=comb(M-1-r,j)*(-1)**j
        # Exact integers, not reduced mod an input N.
        s=[sum(comb(t,k)*c**((k-r)//M)
               for k in range(r,t+1,M)) if t>=r else 0
           for t in range(3*M)]
        for t in range(3*M):
            val=sum(den[i]*s[t-i] for i in range(min(t,M)+1))
            assert val==(num[t] if t<M else 0)
        rows.append({"M":M,"r":r,"c":c,"minimal_linear_order":M,
                     "checked_terms":len(s),"quotient_coordinates":M})
    return rows

def derivative_checks():
    rows=[]
    for p,q in [(3,5),(5,7),(7,13),(11,17),(17,41)]:
        N=p*q;d=gcd(p-1,q-1)
        hits=sum(1<gcd((pow(a+1,N-1,N)-pow(a,N-1,N))%N,N)<N for a in range(N))
        expected=(d-1)*(p+q-2*(d-1))
        assert hits==expected
        roots=[]
        for ell in (p,q):
            # Exact polynomial division by 2X+1 via integer coefficient recurrence.
            coeff=[comb(N-1,k) for k in range(N-1)]
            rem=0
            qq=[]
            # H coefficients ascending. b_0=a_0, b_k=a_k-2b_(k-1).
            prev=0
            for k in range(len(coeff)-1):
                prev=coeff[k]-2*prev
                qq.append(prev)
            assert coeff[-1]==2*prev
            count=sum(sum(v*pow(a,k,ell) for k,v in enumerate(qq))%ell==0 for a in range(ell))
            assert count==d-2
            roots.append(count)
        rows.append({"N":N,"p":p,"q":q,"d":d,"split_count":hits,
                     "probability":str(Fraction(hits,N)),"stripped_root_counts":roots})
    return rows

def torus_checks():
    rows=[]
    for p,q,D in [(5,7,2),(7,11,2),(11,13,2),(10007,10037,2)]:
        assert reference_prime(p) and reference_prime(q)
        N=p*q
        eps=[1 if pow(D,(ell-1)//2,ell)==1 else -1 for ell in (p,q)]
        A=p-eps[0];B=q-eps[1];d=gcd(A,B);m=N-eps[0]*eps[1]
        assert gcd(m,A)==gcd(m,B)==d
        if max(p,q)<100:
            for ell,order in ((p,A),(q,B)):
                group=[(a,b) for a in range(ell) for b in range(ell) if (a*a-D*b*b)%ell==1]
                assert len(group)==order
                def mp(x,y):
                    return ((x[0]*y[0]+D*x[1]*y[1])%ell,
                            (x[0]*y[1]+x[1]*y[0])%ell)
                def pw(x,k):
                    ans=(1,0)
                    while k:
                        if k&1:ans=mp(ans,x)
                        x=mp(x,x);k//=2
                    return ans
                assert sum(pw(x,m)==(1,0) for x in group)==d
        prob=Fraction(d,A)+Fraction(d,B)-Fraction(2*d*d,A*B)
        rows.append({"N":N,"D":D,"orders":[A,B],"d":d,"split_probability":str(prob)})
    return rows

def sieve_checks():
    rows=[]
    primes=[3,5,7,11,13,17,19,23,29,31,37,41,43,47]
    for p,q in [(101,113),(107,193),(179,211),(251,401),(1009,1481)]:
        assert reference_prime(p) and reference_prime(q)
        N=p*q;lo=isqrt(4*N)
        if lo*lo<4*N:lo+=1
        hi=3*isqrt(N)+3
        survivors=list(range(lo,hi+1));period=1;classes=1
        for k,ell in enumerate(primes,1):
            squares={a*a%ell for a in range(ell)}
            mask={s for s in range(ell) if (s*s-4*N)%ell in squares}
            chi=1 if N%ell in squares else -1
            assert len(mask)==(ell+chi)//2
            period*=ell;classes*=len(mask)
            survivors=[s for s in survivors if s%ell in mask]
            assert p+q in survivors
            if k in (4,8,14):
                true=[s for s in survivors if isqrt(s*s-4*N)**2==s*s-4*N]
                rows.append({"N":N,"k":k,"interval":[lo,hi],"period_bits":period.bit_length(),
                             "explicit_CRT_classes":classes,"survivors":len(survivors),
                             "square_discriminants":true,"planted_sum":p+q})
    return rows

def rank_experiments():
    records=[];summaries=[]
    generation=Random(20261017);public=Random(43019)
    for bits in (8,12,16):
        bucket=[]
        for instance in range(12):
            p=random_prime(bits,generation);q=random_prime(bits,generation)
            while q==p:q=random_prime(bits,generation)
            N=p*q
            # Reproducible public sketch distribution; no diagnostic factors used.
            for r in (2,3,4,5,7,8):
                for repeat in range(8):
                    c=public.randrange(2,N)
                    ans=candidate(N,r,c)
                    row={"N":N,"diagnostic_p":p,"diagnostic_q":q,
                         "r":r,"c":c,"answer":ans}
                    if ans["stage"]=="power":
                        mat=matrix(ans["f"],N,r,c)
                        ranks=[field_rank(mat,p),field_rank(mat,q)]
                        row["diagnostic_ranks"]=ranks
                        if ranks[0]!=ranks[1]:assert ans["result"]["status"]=="factor"
                        if ans["result"]["status"]=="factor":
                            assert N%ans["result"]["factor"]==0
                    records.append(row);bucket.append(row)
        power_rows=[x for x in bucket if x["answer"]["stage"]=="power"]
        summaries.append({"factor_bits":bits,"instances":12,"probes":len(bucket),
                          "parameter_or_universal":len(bucket)-len(power_rows),
                          "power_probes":len(power_rows),
                          "power_splits":sum(x["answer"]["result"]["status"]=="factor" for x in power_rows),
                          "rank_gaps":sum(x["diagnostic_ranks"][0]!=x["diagnostic_ranks"][1] for x in power_rows),
                          "matrix_only_splits":sum(x["answer"]["result"]["status"]=="factor" and not x["answer"]["coeff_factor"] for x in power_rows)})
    return records,summaries

def tiny_repair():
    rows=[]
    for N,r,c in [(35,2,2),(77,2,2),(143,2,3)]:
        f=sketch(N,r,c)
        norm=(f[0]*f[0]-c*f[1]*f[1])%N
        rows.append({"N":N,"r":r,"c":c,"f":f,"norm_gcd":gcd(norm,N),
                     "universal_resultant_gcd":gcd(c*c+c+1,N)})
    # After excluding roots of W, N35 never has a rank defect.
    for r in range(2,9):
        for c in range(35):
            if gcd(r*c,35)!=1:continue
            w=unit_probe(35,matrix(universal_mask(35,r,c),35,r,c))
            if w["status"]=="full":
                mat=matrix(sketch(35,r,c),35,r,c)
                assert field_rank(mat,5)==field_rank(mat,7)==r
    return rows

def discriminant_checks():
    rows=[]
    for K in range(2,9):
        d=prod(j**(2*(K-j)) for j in range(1,K))
        rows.append({"K":K,"D_mod_35":d%35,"gcd":gcd(d,35)})
        if K>=3:
            prev=prod(j**(2*(K-1-j)) for j in range(1,K-1))
            nxt=prod(j**(2*(K+1-j)) for j in range(1,K+1))
            assert nxt*prev==K*K*d*d
    return rows


def random_quotient_pit_checks():
    def mm(a,b,N,h):
        e=len(h);v=[0]*(2*e-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):v[i+j]=(v[i+j]+x*y)%N
        for k in range(2*e-2,e-1,-1):
            for j,t in enumerate(h):v[k-e+j]=(v[k-e+j]-v[k]*t)%N
        return v[:e]
    def pw(a,k,N,h):
        e=len(h);v=[1]+[0]*(e-1)
        while k:
            if k&1:v=mm(v,a,N,h)
            a=mm(a,a,N,h);k//=2
        return v
    rng=Random(17391);rows=[]
    for N,kind,D in [(561,"F",561),(17,"F",17),(9,"nilpotent",3),(25,"nilpotent",5)]:
        e=D.bit_length()+2;hits=0
        for trial in range(128):
            h=[rng.randrange(N) for _ in range(e)]
            X=[0,1]+[0]*(e-2)
            if kind=="F":
                a=pw([1,1]+[0]*(e-2),N,N,h);b=pw(X,N,N,h)
                a=[(x-y)%N for x,y in zip(a,b)];a[0]=(a[0]-1)%N
            else:
                # 3(X^3-X) mod9 or5(X^5-X) mod25: zero on ALL scalar residues.
                p=isqrt(N);a=pw(X,p,N,h)
                a=[p*(x-y)%N for x,y in zip(a,X)]
                assert all(p*(pow(t,p,N)-t)%N==0 for t in range(N))
            hits+=any(a)
        if N==17:assert hits==0
        else:assert hits>0
        rows.append({"N":N,"kind":kind,"degree_bound":D,"quotient_degree":e,
                     "trials":128,"nonzero_remainders":hits})
    return rows

def random_gcd_baseline(records):
    rng=Random(19073);rows=[]
    for bits in (8,12,16):
        bucket=[x for x in records if x["diagnostic_p"].bit_length()==bits and x["answer"]["stage"]=="power"]
        expected=0.0;success=0;calls=0
        for x in bucket:
            N=x["N"];p=x["diagnostic_p"];q=x["diagnostic_q"]
            budget=x["answer"]["metrics"]["gcd_calls"];calls+=budget
            per=(p+q-2)/(p*q);expected+=1-(1-per)**budget
            hit=False
            for _ in range(budget):
                d=gcd(rng.randrange(N),N)
                hit |= 1<d<N
            success+=hit
        rows.append({"factor_bits":bits,"trials":len(bucket),"gcd_calls":calls,
                     "successful_random_gcd_trials":success,
                     "independent_random_gcd_expected_trials":expected})
    return rows

def robustness_checks():
    checked=0
    for N in range(2,201):
        for r in (2,3,4):
            for c in (2,3,5):
                a=candidate(N,r,c)
                d=a.get("factor") or a.get("result",{}).get("factor")
                if d:assert 1<d<N and N%d==0,(N,r,c,a)
                checked+=1
    return {"all_inputs_2_through_200":checked,
            "includes":"primes, prime powers, unbalanced, multiple prime factors"}



def exact_quadratic_rank_distributions():
    rows=[]
    for p,q in [(101,113),(107,193),(199,241),(251,401),(1009,1481)]:
        N=p*q;dist=[]
        for ell in (p,q):
            ct=[0,0,0]
            for c in range(ell):
                if c*(1-c)*(c*c+c+1)%ell==0:continue
                f=power([1,1],N,ell,2,c)
                f[0]=(f[0]-1)%ell
                f[1]=(f[1]-pow(c,N//2,ell))%ell
                a,b=f;det=(a*a-c*b*b)%ell
                rank=2 if det else (1 if a or b else 0)
                ct[rank]+=1
            dist.append(ct)
        a,b=dist
        prob=1-sum(Fraction(a[i],sum(a))*Fraction(b[i],sum(b)) for i in range(3))
        rows.append({"N":N,"diagnostic_primes":[p,q],"rank_counts":dist,
                     "exact_conditional_rank_gap_probability":str(prob)})
    return rows


def main():
    started=time.perf_counter()
    out={"seed":20261017,"public_sketch_seed":43019,
         "seed_checks":seed_checks(),"section_checks":section_checks(),
         "derivative_checks":derivative_checks(),"torus_checks":torus_checks(),
         "sieve_checks":sieve_checks(),"discriminant_checks":discriminant_checks(),
         "tiny_rank_repair":tiny_repair()}
    records,summaries=rank_experiments()
    out["rank_summary"]=summaries;out["rank_raw"]=records
    out["random_quotient_pit"]=random_quotient_pit_checks()
    out["matched_random_gcd_baseline"]=random_gcd_baseline(records)
    out["robustness"]=robustness_checks()
    out["exact_quadratic_rank_distributions"]=exact_quadratic_rank_distributions()
    out["metrics"]=MET
    out["wall_seconds"]=time.perf_counter()-started
    out["peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    Path("phase1_results.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:v for k,v in out.items() if k not in ("rank_raw",)},indent=2))
if __name__=="__main__":main()

