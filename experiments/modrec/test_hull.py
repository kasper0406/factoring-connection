import random, sys
from hull_count import count_box, brute_box
rng = random.Random(1)
bad = 0; n = 0
for (b, c) in [(0,1),(1,6),(1,1),(0,5),(1,2)]:
    for _ in range(300):
        T = rng.randrange(0, 4000)
        X = rng.randrange(0, 80)
        n += 1
        if count_box(b,c,T,X) != brute_box(b,c,T,X):
            bad += 1; print('BAD', b,c,T,X, count_box(b,c,T,X), brute_box(b,c,T,X)); 
            if bad > 5: sys.exit()
print('cases', n, 'bad', bad)
