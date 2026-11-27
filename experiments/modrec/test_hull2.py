import random, sys, time
from math import isqrt
from hull_count import count_box, Region
rng = random.Random(2)
def colsum(b,c,T,X):
    K = Region(b,c,T,0,X)
    s = 0
    for x in range(0, min(X, isqrt(T))+1):
        s += min(K.top(x), X) + 1
    return s
bad = 0
for (b,c) in [(0,1),(1,6)]:
    for e in range(12, 37, 2):
        T = rng.getrandbits(e)
        for X in [isqrt(T)//3, isqrt(T//(1+b+c)), isqrt(T)+5]:
            st = {}
            if e <= 30:
                ok = count_box(b,c,T,X,stats=st) == colsum(b,c,T,X)
                bad += not ok
            else:
                count_box(b,c,T,X,stats=st)
print('bad', bad)
