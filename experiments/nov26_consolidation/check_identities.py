"""Fresh small exact checks; not original scratch data or complexity benchmarks."""
from math import gcd, isqrt
from itertools import combinations
from random import Random
import json


def jacobi(a, n):
    a %= n
    s = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                s = -s
        a, n = n, a
        if a % 4 == n % 4 == 3:
            s = -s
        a %= n
    return s if n == 1 else 0


def primes(n):
    return [p for p in range(3, n, 2) if all(p % d for d in range(3, isqrt(p)+1, 2))]


def forms(D):
    for a in range(1, isqrt(D//3)+1):
        for b in range(-a, a+1):
            if (b*b+D) % (4*a):
                continue
            d = (b*b+D)//(4*a)
            if d >= a and gcd(gcd(a,b),d) == 1:
                if b < 0 and (a == d or a == -b):
                    continue
                yield a,b,d


def count(form, B, divisor=1):
    """Nonzero vectors modulo simultaneous negation, by complete-square bounds."""
    a,b,d = form
    D = 4*a*d-b*b
    total = 0
    for y in range(isqrt(4*a*B//D)+1):
        R = isqrt(4*a*B-D*y*y)
        lo = -((R+b*y)//(2*a))
        hi = (R-b*y)//(2*a)
        if y == 0:
            lo = max(lo, 1)
        for x in range(lo, hi+1):
            v = a*x*x+b*x*y+d*y*y
            total += int(v > 0 and v % divisor == 0)
    return total


def first(lo, hi, pred):
    assert pred(hi)
    while hi-lo > 1:
        mid = (lo+hi)//2
        if pred(mid):
            hi = mid
        else:
            lo = mid
    return hi


def decode(D, oracle):
    hi = 1
    while oracle(hi) % 4 == 0:
        hi *= 2
    a = first(0, hi, lambda t: oracle(t) % 4 != 0)
    hi = (D+a*a)//(4*a)
    d = first(a-1, hi, lambda t: (oracle(t)-isqrt(t//a)) % 4 != 0)
    b = isqrt(4*a*d-D)
    assert b*b == 4*a*d-D
    return a,b,d


def extract(N, ell, c, j, u, v):
    H = ell*ell*j+c*u*v
    disc = H*H-4*N*u*u*v*v
    if disc < 0:
        return None
    z = isqrt(disc)
    if z*z != disc:
        return None
    for w in (H-z,H+z):
        k,rem = divmod(w,2*v*v)
        if rem == 0 and 1 < k < N and N % k == 0:
            return k
    return None


def run():
    out = dict(quarter=0, reflection=0, half_class=0, central_class=0,
               slopes=0, grid_blocks=0, forms=0, canonical_hits=0,
               ellipse_cases=0, ellipse_queries=0)
    P = primes(100)
    for p,q in combinations(P,2):
        N = p*q
        chi = [jacobi(x,N) for x in range(N)]
        K = [0]
        for v in chi[1:]:
            K.append(K[-1]+(v == -1))
        if sorted([p%8,q%8]) == [1,3]:
            assert K[N//4] == (p-1)*(q-1)//8
            out['quarter'] += 1
        if N%4 == 3:
            for y in range(isqrt(N)+1):
                assert (y+K[y]+K[N-1-y])%2 == (y//p)%2
                out['reflection'] += 1
        if N%24 == 19:
            h = len(list(forms(24*N)))
            assert h == 4*sum(chi[5*N//24+1:7*N//24+1])
            out['central_class'] += 1
        if N%8 == 1:
            h = len(list(forms(8*N)))
            assert h == 4*sum(chi[:N//8+1])-2*sum(chi[:N//4+1])
            out['half_class'] += 1
    for k in range(3,13):
        m = 1 << k
        for n in (1,3,9,31,59):
            assert len({n*pow(u,-2,m)%m for u in range(1,m,2)}) == m//8
            out['slopes'] += 1
    rng = Random(20261126)
    for _ in range(700):
        p,q = rng.sample(P,2); N=p*q
        x,y,z,T = [rng.choice([i for i in range(1,N) if gcd(i,N)==1]) for _ in range(4)]
        u,v,L = rng.randrange(1,20),rng.randrange(1,20),rng.randrange(1,10)
        def g(u,v): return pow(x,u*u,N)*pow(y,v*v,N)*pow(z,u*v,N)%N
        A=g(u,v); B=pow(y,2*v,N)*pow(z,u,N)%N
        f=h=1
        for i in range(L):
            w=A*pow(B,i,N)*pow(y,i*i,N)%N
            assert w == g(u,v+i)
            f=f*(1-w*pow(T,-1,N))%N
            h=h*(1-g(u,v+i)*pow(T,-1,N))%N
        assert f == h
        out['grid_blocks'] += 1
    for D in range(5,2001):
        for a,b,d in forms(D):
            assert decode(D,lambda B:count((a,b,d),B)%4) == (a,abs(b),d)
            out['forms'] += 1
    # Same Delta and squared class, different canonical minima; mixed CRT root is diagnostic only.
    for N,p,q,c,A_expected in [(11009,101,109,204,1),(28829,127,227,336,2)]:
        ell=11; D=(4*N-c*c)//(ell*ell)
        assert 4*N-c*c == D*ell*ell and D == 20
        root=c*pow(ell,-1,N)%N
        a_root=root%p; b_root=-root%q
        t=(a_root+p*((b_root-a_root)*pow(p,-1,q)%q))%N
        if t%2 != D%2: t+=N
        assert (t*t+D)%(4*N)==0
        mixed=(N,t,(t*t+D)//(4*N))
        Q0=(1,D%2,(D+(D%2)**2)//4)
        for B in range(1,21):
            assert count(Q0,N*B,N) == 2*count(Q0,B)+2*count(mixed,B)
            out['ellipse_queries'] += 1
        local=root%p
        if local%2 != D%2: local+=p
        hidden=(p,local,(local*local+D)//(4*p))
        assert next(B for B in range(1,20) if count(hidden,B)) == A_expected
        assert next(B for B in range(1,20) if count(mixed,B)) == 1
        def oracle(B): return (count(Q0,N*B,N)-2*count(Q0,B))//2%4
        assert decode(D,oracle)==(1,0,5)
        out['ellipse_cases'] += 1
        for u in range(1,25):
            for v in range(1,25):
                if gcd(u,v)==1 and (2*p*v-c*u)%ell==0:
                    j=(p*v*v+q*u*u-c*u*v)//ell**2
                    assert (p*v*v+q*u*u-c*u*v)%ell**2==0
                    assert extract(N,ell,c,j,u,v)==p
                    if 4*j*j<D: assert j==A_expected
                    out['canonical_hits'] += 1
    return out


if __name__ == '__main__':
    print(json.dumps(run(),indent=2))
