# Three elementary checks for the strategy review

These are deductions made during the read-only review, not changes to the paper.

## 1. Prime presieving in the arithmetic large-sieve denominator

Let `P = 2^(n/2)`, `s = |V(H)|`, and fix a prime q in the theorem's q-range and one band k.
Assume the scale separations in Theorem D′. Let U consist of **prime** p in [P/2,P]
uncovered by this band. Let `Q = floor(sqrt(P))` and `Z = floor(Q/(32Y))`.
For all sufficiently large n, `1 < Z < Y`, `log Z` is comparable to n, and every
band prime is less than Q. The following avoided classes are available:

- modulo each band prime ell, the s distinct classes Aq/B;
- modulo every prime r <= Z, the class 0, since every p in U is greater than r.

The arithmetic large sieve gives

    |U| <= (P/2 + 1 + Q^2)/J,
    J = sum_{d<=Q} mu(d)^2 product_{r|d} omega(r)/(r-omega(r)).

Retain the terms `d = ell*a`, with ell in the chosen band and squarefree `a <= Z`.
They satisfy d <= Q, have no band prime dividing a, and are all distinct as (ell,a)
varies. Thus

    J >= [sum_ell s/(ell-s)] [sum_{a<=Z} mu(a)^2/phi(a)]
      >= c lambda log Z
      >= c' lambda n,

where `lambda = s sum_ell 1/ell`. The elementary estimate in the second line follows
already from `1/phi(a) >= 1/a` and the positive density of squarefree integers.
Consequently `|U| << P/(lambda*n)`, uniformly in q and the choice of band.

There are O(P/n) prime q, and the bounded-ratio population has size comparable to
P^2/n^2. Choosing the clean band as a function of q therefore gives the coverage
exceptional fraction

    O(1/lambda) = O(n/H^2).

The other exceptional terms in D′ are unchanged. The original n^2/H^2 coverage bound
is valid but loses a logarithm by counting all integer p and retaining only prime
terms in J. This refinement is an application of the standard arithmetic large sieve,
not a new sieve theorem. It does not prove a higher inverse power of lambda.

## 2. The paper's algorithmic DCH includes a counterexample

The conjecture in paper/main.tex includes all sufficiently large N coprime to 6,
all integer candidates p in [sqrt(N/2),sqrt(N)], and excludes A=B.

Take N=t^2 with t a large prime greater than 3 and candidate p=t. Then

    Bp^2 - AN = (B-A)t^2.

For the conjecture's fixed polylogarithmic H, asymptotically every band prime ell
satisfies `H < ell < t`. Since `0 < |B-A| < H`, none can divide this expression.
This disproves that statement as written. It does not affect D′, whose success
population consists of products of distinct primes. Excluding square inputs repairs
this immediate counterexample; it does not prove the repaired conjecture.

## 3. Ratio-dependent work and the Pollard–Strassen envelope

Write R=q/p=N^theta. The proposed localized-weight regime has work

    N^(1/6) R^(1/6) Delta^(2/3),

where Delta is the local ratio width and equals H^2 in the stated design. With
subpolynomial Delta its exponent is `(1+theta)/6`. Pollard–Strassen up to the smaller
factor costs `sqrt(p) = N^((1-theta)/4)` up to logarithms. Their lower envelope is

    min((1+theta)/6, (1-theta)/4),

whose maximum on [0,1/3] is 1/5, attained at theta=1/5. The interval on which both
exponents exceed 1/6 is `0 < theta < 1/3`.

This verifies an optimization of the proposed costs, not a theorem that every
almost-all factoring algorithm must pay this cost. Harvey–Hittmeir supplies the
rigorous 1/5 fallback if the proposed ratio-extension success proof is unavailable.
