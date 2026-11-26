import ctypes
import json
import math
import time
from pathlib import Path
from flint import fmpq, fmpz, fmpz_poly, nmod_poly
from experiment import native, W, HERE

class NMod(ctypes.Structure):
    _fields_ = [('n', W), ('ninv', W), ('norm', W)]

class Poly(ctypes.Structure):
    _fields_ = [('coeffs', ctypes.POINTER(W)), ('alloc', ctypes.c_long),
                ('length', ctypes.c_long), ('mod', NMod)]

ptr = ctypes.POINTER(Poly)
vec = ctypes.POINTER(W)
init = native("nmod_poly_init", None, ptr, W)
clear = native("nmod_poly_clear", None, ptr)
product_roots = native("nmod_poly_product_roots_nmod_vec", None,ptr,vec,ctypes.c_long)
evaluate = native("nmod_poly_evaluate_nmod_vec_fast", None,vec,ptr,vec,ctypes.c_long)

def progression(N,bits):
    L = 1 << (bits//2-1)
    first = L+(3-L)%4
    length = (2*L-1-first)//4+1
    residue = block_product(N,first,length)
    if bits == 32:
        z = 1
        for i in range(first,2*L,4):
            z = z*i % N
        assert z == residue
    return residue

def block_product(N,first,length):
    if not length:
        return 1
    width = math.isqrt(length)
    blocks, tail = divmod(length, width)
    # F(X) = product (X+4j), j = 0,...,width-1.
    roots = (W*width)(*((-4*j)%N for j in range(width)))
    points = (W*blocks)(*(first+4*j*width for j in range(blocks)))
    values = (W*blocks)()
    F = Poly()
    init(ctypes.byref(F), N)
    try:
        product_roots(ctypes.byref(F),roots,width)
        evaluate(values,ctypes.byref(F),points,blocks)
    finally:
        clear(ctypes.byref(F))
    residue = 1
    for value in values:
        residue = residue*value % N
    for i in range(blocks*width,length):
        residue = residue*(first+4*i) % N
    return residue

def controls():
    result = []
    for bits in (32,48,64):
        for line in (HERE/f'rows{bits}.jsonl').read_text().splitlines():
            row = json.loads(line)
            if row['index'] >= 4:
                continue
            N = row['N']
            start = time.perf_counter()
            z = progression(N,bits)
            elapsed = time.perf_counter()-start
            g = math.gcd(z,N)
            assert g == (N if row['c'] == 3 else 1)
            result.append(dict(bits=bits,N=N,c=row['c'],residue=z,g=g,seconds=elapsed))
    return result

def truncated_power(poly,k,prec):
    out = fmpz_poly([1])
    while k:
        if k & 1:
            out = (out*poly).truncate(prec)
        poly = (poly*poly).truncate(prec)
        k //= 2
    return out

def sigma(n,k):
    return sum(d**k for d in range(1,n+1) if n%d==0)

def pure_e4(k):
    m = k//4
    h = m//3
    E4 = fmpz_poly([1]+[240*sigma(n,3) for n in range(1,h+1)])
    E6 = fmpz_poly([1]+[-504*sigma(n,5) for n in range(1,h+1)])
    delta = ((E4**3-E6**2)//1728).truncate(h+1)
    multiplier = -2*k / fmpq.bernoulli(k)
    series = [fmpq(1)]+[multiplier*sigma(n,k-1) for n in range(1,h+1)]
    coeff = 0
    for j in range(h+1):
        base = (truncated_power(E4,m-3*j,h+1)*truncated_power(delta,j,h+1)).truncate(h+1)
        c = series[j]
        assert base[j] == 1
        coeff += c/1728**j
        for n in range(j,h+1):
            series[n] -= c*base[n]
    assert not any(series)
    return int(coeff.p),int(coeff.q)

def cyclotomic_gcd(N):
    # q-Pascal modulo X^N-1. The auxiliary prime exceeds 2m.
    m = N//4
    ell = 1009
    one = nmod_poly([1],ell)
    H = one.left_shift(N)-one
    row = [one]+[nmod_poly([],ell) for _ in range(m)]
    for j in range(1,2*m+1):
        for k in range(min(j,m),0,-1):
            row[k] = (row[k]+row[k-1].left_shift(j-k))%H
    return row[m].gcd(H//(nmod_poly([0,1],ell)-one)).degree()

def relation(N,s):
    m = N//4
    # A solution value is a scalar times its most recent free coordinate.
    coordinate, scale = 0,1
    for j in range(1,2*m+1):
        assert (2*(2*m+1-j) + 2*j-1)%s == 0
        if j%s == 0:
            assert scale == 0
            coordinate,scale = coordinate+1,1
        else:
            scale = scale*(2*j-1)*pow(2*j,-1,s)%s
        if j == m:
            middle_zero = scale == 0
    return middle_zero, coordinate>0 and scale!=0

def small():
    primes = [s for s in range(3,134,2) if fmpz(s).is_prime()]
    rows=[]
    for p in primes:
        for q in primes:
            N=p*q
            if p>=q or p%4 != q%4 or N>=400:
                continue
            m = N//4
            num,den = pure_e4(N-1)
            degree = cyclotomic_gcd(N)
            expected = sum(d-1 for d in (p,q) if 2*(m%d)>=d)
            assert degree == expected == ((p+q-2) if p%4==3 else 0)
            for s in (p,q):
                assert relation(N,s) == (s%4==3,True)
            rows.append(dict(N=N,p=p,q=q,num=str(num),den=str(den),gnum=math.gcd(num,N),
                             gden=math.gcd(den,N),cyclo_degree=degree))
    return sorted(rows,key=lambda r:r['N'])

if __name__ == '__main__':
    data = dict(small=small(), control=controls())
    (HERE/'structure.json').write_text(json.dumps(data,indent=2)+'\n')
    print('Small pairs',len(data['small']))
    print('Denominator nonunits',[(r['N'],r['gden']) for r in data['small'] if r['gden']!=1])
    print('Defined mismatches',[(r['N'],r['gnum']) for r in data['small'] if r['gden']==1 and ((r['gnum']==r['N'])!=(r['p']%4==3))])
    print('First small cases',[(r['N'],r['num'],r['den']) for r in data['small'][:3]])
    print('Control max seconds',max(r['seconds'] for r in data['control']))
