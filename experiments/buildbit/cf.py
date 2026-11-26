import ctypes as ct
import json
import sys
from math import gcd, isqrt, prod
from pathlib import Path

MULT = (1,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59)
SIZE = 16
lib = ct.CDLL(str(Path(__file__).resolve().with_name('basis.so')))
Word = ct.c_uint64
lib.select_multipliers.argtypes = (ct.POINTER(Word), ct.POINTER(ct.c_int))
lib.square_columns.argtypes = (ct.POINTER(Word), ct.c_int, ct.POINTER(Word), ct.c_int)
lib.square_columns.restype = ct.c_int

def columns(den):
    n = len(den)
    words = 5
    while True:
        data = (Word * (words*n))()
        if lib.square_columns((Word*n)(*den),n,data,words) >= 0:
            break
        words *= 2
    buf = memoryview(data).cast('B')
    return [int.from_bytes(buf[i*8*words:(i+1)*8*words],sys.byteorder) for i in range(n)]

def trajectory(N,k,steps=SIZE):
    D = N*k
    a0 = isqrt(D)
    a, P, Q = a0, 0, 1
    older, last = 0, 1
    norms, roots = [], []
    for step in range(steps):
        h = (a*last + older) % N
        P = a*Q-P
        Q = (D-P*P)//Q
        norms.append(Q)
        roots.append(h)
        older, last = last,h
        a = (a0+P)//Q
    return norms, roots

def dependencies(row, start=0):
    pivots = {}
    for i,v in enumerate(row,start):
        support = 1 << i
        while v:
            pivot = v.bit_length()
            old = pivots.get(pivot)
            if old is None:
                pivots[pivot] = v,support
                break
            v ^= old[0]
            support ^= old[1]
        if not v:
            yield i,support

def evaluate(N):
    # All selection reads only Q1 from unselected trajectories.
    paths = [trajectory(N,k) for k in MULT]
    candidates = (Word*(SIZE*len(paths)))(*(d for norm,_ in paths for d in norm))
    order = (ct.c_int*12)()
    lib.select_multipliers(candidates,order)
    den = [d for k in order for d in paths[k][0]]
    roots = [h for k in order for h in paths[k][1]]
    flags = set()
    splits = {}
    def log_split(g,via):
        if 1 < g < N:
            splits.setdefault(via,g)
    for k in range(6,len(paths)):
        log_split(gcd(paths[k][0][0],N),'candidate')
    # Check each available modular square root and all selected denominators.
    for j,(v,h) in enumerate(zip(den, roots)):
        assert (h*h - (v if j%2 else -v)) % N == 0
        log_split(gcd(v,N),'denominator')
    mask = sum(1 << j for j in range(0,len(den),2))
    row = columns(den)
    individual = []
    for i,k in enumerate(order):
        if any((mask&s).bit_count()%2 for _,s in dependencies(row[SIZE*i:SIZE*(i+1)],SIZE*i)):
            individual.append(MULT[k])
    witnesses = []
    first_root = None
    for j,s in dependencies(row):
        odd = (mask&s).bit_count()%2
        if odd:
            flags.add('full')
            if j < 6*SIZE:
                flags.add('first')
        indices = []
        while s:
            a = s & -s
            indices.append(a.bit_length()-1)
            s ^= a
        qprod = prod(den[i] for i in indices)
        v = isqrt(qprod)
        assert v*v == qprod
        g = gcd(v,N)
        log_split(g,'relation_denominator')
        if g != 1:
            continue
        h = prod(roots[i] for i in indices) * pow(v,-1,N) % N
        assert h*h % N == (N-1 if odd else 1)
        if odd:
            if first_root is None:
                first_root = h
                witnesses.append({'root': h, 'terms': [(MULT[order[i//SIZE]],i%SIZE+1) for i in indices]})
            else:
                log_split(gcd(h-first_root,N),'two_negative_roots')
                log_split(gcd(h+first_root,N),'two_negative_roots')
        else:
            log_split(gcd(h-1,N),'positive_root')
            log_split(gcd(h+1,N),'positive_root')
    return dict(N=N,first='first' in flags,full='full' in flags,individual=individual,splits=splits,
                order=[MULT[k] for k in order],witness=witnesses)

def run(n):
    path = Path(f'cf{n}.tsv')
    with Path(f'cf{n}_out.jsonl').open('w') as out:
        for line in path.read_text().splitlines():
            N,p,q,label = map(int,line.split())
            row = evaluate(N)
            assert not row['full'] or label or row['splits']
            assert all(z in (p,q) for z in row['splits'].values())
            row['label'] = label
            out.write(json.dumps(row,separators=(',',':'))+'\n')

if __name__ == '__main__':
    run(int(sys.argv[1]))
