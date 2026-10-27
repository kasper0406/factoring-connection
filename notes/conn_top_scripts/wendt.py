# codex's generalized Wendt resultant D_N(c)=Res(X^N-1,(X+c)^N-1) mod N; when does gcd split N?
# prediction (ours): r | D_s(c) needs mu_s - mu_s to contain c in F_{r^k}, k = ord_s(r): only when k is small (M/A event).
import flint
from math import gcd
from sympy import primerange, n_order
def res_mod_prime(N,c,r):
    ctx=flint.fmpz_mod_poly_ctx(r)
    f=ctx([r-1]+[0]*(N-1)+[1]); g=ctx([c,1])**N - ctx([1])
    return int(f.resultant(g))
primes=list(primerange(11,80))
tot=0; split=0; predicted=0; agree=0
for i,p in enumerate(primes):
    for q in primes[i+1:]:
        N=p*q
        for c in (1,2,3):
            a=res_mod_prime(N,c,p)==0; b=res_mod_prime(N,c,q)==0
            tot+=1
            if a!=b: split+=1
            kp=n_order(q,p); kq=n_order(p,q)   # ord_p(q): mu_p in F_{q^kp}... careful
            # p | D_q(c) involves mu_q in F_{p^k}, k=ord_q(p); q | D_p(c) involves mu_p in F_{q^k}, k=ord_p(q)
            small = min(n_order(p,q), n_order(q,p))<=3
            if a!=b: predicted+= small
print("pairs*c",tot,"asymmetric splits",split,"of which min order<=3:",predicted)
