"""Test public integer-representative features on held-out moduli."""

import json
from math import gcd, isqrt, log, sqrt
from random import Random

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from adaptive_exponent_screen import fixtures
from check_quarter_fourier import quarter


TIERS = (36, 64)
TRAIN_MODULI = 16
TEST_MODULI = 8
SAMPLES = 4096
SEED = 261402


def jacobi(a, n):
    a %= n
    result = 1
    while a:
        parity = (a & -a).bit_length()-1
        a >>= parity
        if parity % 2 and n % 8 in (3,5):
            result = -result
        if a % 4 == n % 4 == 3:
            result = -result
        a, n = n % a, a
    return result if n == 1 else 0


def public_features(n, c):
    inverse = pow(c,-1,n)
    c2 = c*c % n
    c3 = c*c2 % n
    values = (c, inverse, c2, c3, pow(c,c,n), pow(c,n-1,n), pow(c,65_537,n))
    features = []
    for v in values:
        features.extend((v/n, min(v,n-v)/n, v.bit_count()/n.bit_length(),
                         isqrt(v)/isqrt(n), quarter(n,v),
                         jacobi(v-1,n), jacobi(v+1,n)))
        features.extend((v >> i) & 1 for i in range(6))
    features.extend(((c*inverse-1)//n/n, c*c//n/n, n%c/n,
                     n.bit_length(), n % 16, jacobi(c,n)))
    x,y = n,c
    previous, denominator = 0,1
    cutoff, last = isqrt(n),1
    quotients = []
    while y:
        quotient, remainder = divmod(x,y)
        quotients.append(quotient)
        x,y = y,remainder
        previous,denominator = denominator,quotient*denominator+previous
        if denominator <= cutoff:
            last = denominator
    features.extend((len(quotients), last/cutoff, c*last%n/n))
    features.extend(log(1+a) for a in (quotients+[0]*6)[:6])
    return features


def gather(primes, rng):
    rows, labels, planted, moduli = [], [], [], []
    for p,q in primes:
        n = p*q
        hidden = 1+p*(-2*pow(p,-1,q) % q)
        moduli.append(n)
        for _ in range(SAMPLES):
            c = rng.randrange(1,n)
            while gcd(c,n) != 1:
                c = rng.randrange(1,n)
            label = quarter(n,c*hidden) > 0
            rows.append(public_features(n,c))
            labels.append(label)
            planted.append(int(label if rng.random()<.55 else not label))
    return np.array(rows,dtype=np.float32), np.array(labels), np.array(planted), moduli


def main():
    fixture_rng, sample_rng = Random(SEED), Random(SEED+1)
    results = {}
    with threadpool_limits(limits=2):
        for bits in TIERS:
            pairs = list(fixtures(bits,fixture_rng,TRAIN_MODULI+TEST_MODULI))
            train, train_y, train_signal, train_moduli = gather(pairs[:TRAIN_MODULI],sample_rng)
            test, test_y, test_signal, test_moduli = gather(pairs[TRAIN_MODULI:],sample_rng)
            record = {"features":train.shape[1], "train_moduli":train_moduli,
                      "test_moduli":test_moduli, "correct":{}}
            for name, build in (
                ("logistic",lambda:make_pipeline(StandardScaler(),LogisticRegression(max_iter=200))),
                ("boosted",lambda:HistGradientBoostingClassifier(max_iter=90,max_leaf_nodes=23,
                                                                    early_stopping=False,random_state=SEED)),
            ):
                for plant in (False,True):
                    model = build()
                    features = np.column_stack((train,train_signal)) if plant else train
                    heldout = np.column_stack((test,test_signal)) if plant else test
                    model.fit(features,train_y)
                    correct = model.predict(heldout)==test_y
                    record["correct"][name+("_planted" if plant else "")] = (
                        correct.reshape(TEST_MODULI,SAMPLES).sum(axis=1).tolist()
                    )
            record["correct"]["constant"] = test_y.reshape(TEST_MODULI,SAMPLES).sum(axis=1).tolist()
            results[bits] = record
    radius = sqrt(log(2*len(TIERS)*5/.05)/(2*TEST_MODULI*SAMPLES))
    print(json.dumps({"seed":SEED,"samples":SAMPLES,"radius":radius,"tiers":results},indent=2))


if __name__ == "__main__":
    main()
