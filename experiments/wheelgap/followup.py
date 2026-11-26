import ctypes as ct
import json
import math
import numpy as np
from collections import defaultdict
from pathlib import Path
from time import process_time
from benchmark import signature, ODD

lib=ct.CDLL(str(Path('followup.dylib').resolve()))
lib.first_integer.argtypes=[ct.POINTER(ct.c_int),ct.c_int,ct.c_uint32,ct.c_uint64]
lib.first_integer.restype=ct.c_uint64
lib.first_rational.argtypes=[ct.POINTER(ct.c_int),ct.c_int,ct.c_uint32,ct.POINTER(ct.c_int),ct.POINTER(ct.c_int)]
lib.first_rational.restype=None

def collision():
    P=ODD[:31]
    N=2**22
    sieve=np.ones(N+1,dtype=bool)
    sieve[:2]=False
    for p in range(2,math.isqrt(N)+1):
        if sieve[p]: sieve[p*p::p]=False
    ns=np.flatnonzero(sieve & (np.arange(N+1)>131))
    signs=np.zeros(len(ns),dtype=np.uint32)
    for i,p in enumerate(P):
        qr=np.ones(p,dtype=np.uint32)
        qr[(np.arange(1,p)**2)%p]=0
        signs|=qr[ns%p]<<i
    for with8 in (False,True):
        seen={}
        for x,s in zip(ns.tolist(),signs.tolist()):
            key=(x%8 if with8 else 0,s)
            if key in seen:
                y=seen[key]
                assert signature(x,P)==signature(y,P)
                print('prime collision',with8,y,x,'mod8',y%8,x%8,flush=True)
                break
            seen[key]=x
        else: print('no prime collision',with8,flush=True)

def run():
    collision()
    data=defaultdict(dict)
    for line in open('trials.jsonl'):
        r=json.loads(line)
        data[r['bits'],r['typ'],r['trial']][r['method']]=r
    with open('exact.jsonl','w',buffering=1) as f:
        for (bits,typ,j),methods in data.items():
            r=methods['lll']
            wheel=r['wheel']
            k=len(wheel)
            p=(ct.c_int*k)(*wheel)
            a,b=ct.c_int(),ct.c_int()
            start=process_time()
            lib.first_rational(p,k,r['t'],ct.byref(a),ct.byref(b))
            rat_time=process_time()-start
            assert a.value and signature(a.value*b.value,wheel)==r['t']
            start=process_time()
            h=lib.first_integer(p,k,r['t'],min((1<<50),a.value*b.value,r['h']))
            int_time=process_time()-start
            assert h and signature(h,wheel)==r['t']
            scan,xor=methods['scan'],methods['xor']
            assert (scan['h']==h if scan['h'] else scan['prefix']<h)
            assert (max(xor['a'],xor['b'])==a.value if xor['h'] else xor['prefix']<a.value)
            f.write(json.dumps(dict(bits=bits,typ=typ,trial=j,h=h,a=a.value,b=b.value,int_secs=int_time,rat_secs=rat_time))+'\n')
            if j==31 and typ=='S': print('completed',bits,flush=True)
if __name__=='__main__':run()
