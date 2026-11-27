# The unit-coefficient mean extracts an integer

Use the literal all-scale squared quartic transform at zero offset in
[the quadratic chirp calculation](quadratic-chirp.md). Average its value
uniformly over unit outer coefficients A modulo M=N², for odd N>1.
The exact mean is the number H of ordered collisions of the quartic modulo N.

Indeed D=(4A)⁻¹ permutes the units. In the local product formula, each
nonconstant term contains a primitive additive character on at least one
prime-power component of M. Its average on that component's units vanishes:
group units modulo p^(2e) by their residue modulo p^(2e-1) and sum the p
equally spaced phases of the lifts. CRT makes these averages independent.
Only the product of the local collision counts survives. This averaging is
over the outer coefficient; the inner transform still includes nonunit scales.

For distinct odd primes N=pq the constant itself splits N. Write
a=5+χ_p(-1), b=5+χ_q(-1), so a,b belong to {4,6}. Then
H=(3p-a)(3q-b). For each of four guesses form S=(9N+ab-H)/3,
discarding nonintegral values. For the correct guess S=bp+aq, and
S²-4abN=(bp-aq)². Taking the integer square root and trying both signs
in (S±√(S²-4abN))/(2b) recovers p. Verify every candidate by division.
An absolute error strictly below 1/2 in the mean suffices after rounding.
Conversely the factors give the exact mean directly.

There is a quantitative cost to uniform sampling even on balanced inputs.
Let X,Y be the independent primitive phases for the two components.
Each has mean zero and modulus one. The centered transform is
4C_q X + 4C_p Y + 16XY. Orthogonality gives exact complex variance
16(C_p²+C_q²)+256. The ordinary average of m independent exact samples
has mean squared absolute error equal to this quantity divided by m.
Thus constant RMS accuracy costs Θ(N) independent samples when p and q
are comparable. This is a bound for that estimator; it leaves open
structured summation and other uses of individual samples.

For arbitrary odd N the mean identity still holds. If N has a repeated
prime factor p then p divides H by the local formula, so gcd(N,H) can
already split some such inputs. It may equal N. The four-guess decoder
above requires exactly two distinct primes of exponent one.

## Deterministic decoding on all odd inputs using power queries

Let N have k distinct prime divisors. Write e for the exponent of p in N
and χ for the Legendre symbol of −1 modulo p. Substituting ej for e in
the branch's prime-power collision formula gives the polynomial
$$
F(j)=H(N^j)/N^{j-1}
=\prod_{p^e\parallel N}p^{e-1}\big[5e(p-1)j-(2p+\chi)\big].
$$
The equality is for positive integer j. The right side defines the polynomial
at other arguments. Its degree is k and all its roots are rational. Interpolate
from j=1 through floor(log₂ N)+1. The moduli have quadratic bit length
in the bit length of N, and the counts are bounded by the squares of those
moduli. Exact interpolation and deterministic rational polynomial factorization
(by the LLL polynomial factoring algorithm) take polynomial bit time.

For each distinct rational root r try each exponent e up to floor(log₂ N)
and both signs χ. Skip 5er=2. Test the integer candidates
p=(5er+χ)/(5er−2) for divisibility into N, retaining divisors at least 3.
The true root satisfies 5er−2=(2+χ)/(p−1)>0 and supplies the prime p.
Coincident roots cause no ambiguity because all exponent/sign guesses are
tested. Strip the candidates in increasing order with multiplicity. Any
spurious composite divisor has already lost every one of its prime factors
when reached. For even inputs first strip the powers of two.

Thus a polynomial-time exact mean evaluator on arbitrary odd moduli is
deterministically polynomial-time Turing equivalent to factoring. Factoring
computes the mean via the local product formula. The forward reduction also
works with absolute errors strictly below N^(j−1)/2 in H(N^j), rounding
after division. It relies on availability at nonsquarefree powers and does
not give an algorithm for their means from unfactored input. The reference
decoder uses SymPy for exact rational roots. [Implementation](quartic_collision_powers.py)

## A small-output evaluator and the first informative binary digit

Consider an oracle returning H(m) modulo a supplied prime ℓ ≥ 7 dividing m
exactly once. For an odd target n, select auxiliary primes ℓ not dividing n
and take m=nℓ. CRT on the pairs of coordinates gives H(m)=H(n)H(ℓ).
The multiplier H(ℓ)=3ℓ−5−χ_ℓ(−1) is −4 or −6 modulo ℓ, hence invertible.
Recover H(n) modulo each ℓ and stop when the product exceeds n², the number
of ordered coordinate pairs. This determines H(n) exactly.

The product of primes up to x has logarithm Θ(x) by Chebyshev's bounds.
Removing 2,3,5 and the primes dividing n costs at most log(30n) from that
logarithm. Thus the largest required prime is O(log(2n)), and simple sieving
or trial division, oracle queries, and incremental CRT all take polynomial
bit time outside the oracle. The oracle returns O(log log(3n)) bits per
call. Its input includes the auxiliary prime, so identifying it is free.
Combining with the power decoder gives a deterministic all-input reduction.
Conversely the local formula computes these residues in polynomial bit time
given a factorization. This also proves equivalence for H(m) mod m, which
can answer each small-prime query by reduction. The supplied divisor ℓ is
coprime to the target; the remaining cofactor can still be hard to evaluate.

At the fixed binary modulus 4 the answer actually is public: H(n)=1 mod 4
for every odd n. Substitute p=1 or 3 mod 4 into the local prime-power formula,
then multiply. At modulus 8, for squarefree odd n, each p=1,3 mod 8 gives
local count 5, and each p=5,7 mod 8 gives 1. Thus with b=(H(n)−1)/4 mod 2,
one has μ(n)=(-1)^b (−2|n), where the last factor is the public Jacobi symbol.
Indeed b counts the primes with (−2|p)=+1 modulo two, and the Jacobi symbol
counts the complementary primes by sign. This gives a precise interpretation
of one fixed bit. The CRT reduction uses varying primes and establishes no
factoring algorithm from that fixed bit alone.

## A zero test suffices, even without the nonzero residues

Let b be the bit length of odd n ≥ 3. Choose primes ℓ > max(b,5) coprime
to n. For j=1,...,ℓ ask only whether ℓ divides H(n^j). If the input format
requires ℓ to divide the queried modulus once, query n^j ℓ instead: H(ℓ)
is a unit modulo ℓ. The factor n^(j−1) in the power formula is also a unit.
The resulting set of positive answers, read modulo ℓ, is exactly the union
of the roots of the local linear factors of F. A factor that vanishes
identically makes the whole set positive; skip that prime.

Write s for the number of distinct *rational* roots. For a non-skipped prime,
a primitive integral linear factor reduces to either a linear polynomial or
a nonzero constant. Consequently the observed number of distinct roots is
at most s. If exactly s appear, all the rational roots have denominators
invertible modulo ℓ and remain distinct; hence the monic polynomial on
the observed set is the reduction of the same monic squarefree rational
polynomial at every such prime. Multiplicities need never be inferred.

Here are explicit bounds for discovering s without knowing it in advance.
Put L=5bn. The unreduced local root numerators and denominators have absolute
value ≤ L. Outside primes dividing n, a local denominator, or a nonzero
cross determinant of two distinct roots, exactly s roots appear.
The product of all exceptional primes divides the absolute product of those
nonzero integers, which is at most
B=n L^b (2L²)^(b(b−1)/2).
The numerator and denominator of each reduced monic coefficient are bounded
by U=(2L)^b, by expanding the product of the primitive linear factors.
Keep primes attaining the largest observed root count, resetting their CRT
accumulator whenever this count increases. Once the retained prime product
exceeds max(B,2U²), the observed count must be s. Coefficientwise rational
reconstruction is then unique. The old rational-root decoder recovers n's
prime factors, including exponents. Removing all primes below max(b,5) and
all exceptional primes costs only polynomially many bits of prime product.
Chebyshev's bounds therefore give a largest required prime polynomial in b.
All ℓ exponents can be queried with moduli of polynomial bit length.

For squarefree n the native test now has a concrete equivalent formulation.
For each prime p dividing n, with χ determined by p mod 4, the local condition
is (5j−2)p = 5j+χ modulo ℓ. When 5j=2 modulo ℓ there are no hits since
2+χ is a unit. Otherwise the two signs specify at most two residue classes
modulo 4ℓ by CRT. The question is whether n has a prime divisor in their
union. Testing known prime factors in these classes answers the question;
testing n's own residue generally does not. A polynomial-time uniform native
decision procedure for these power queries is still missing.
[One-bit reference decoder](quartic_collision_zero.py)

## Unfactored reference evaluation and an infinite family of misses

[Native zero test](quartic_collision_native.py) searches for any prime
component satisfying the local linear congruence. On a composite part m
of the original odd input n, try every possible exponent e in n. Combine
the local condition modulo ℓ with the sign condition p mod 4. Discard
zero and incompatible residues. For a resulting unit class a modulo 4ℓ,
search up to min(√m, n^(1/e)) in class a. If the exponent permits p>√m,
also search up to √m in class m/a. If a witnessing prime divides m,
either p or m/p occurs in the searched classes. A nonunit gcd gives a
proper split because the candidates are between 2 and √m. A miss excludes
every witness in m. Recurse on both parts of any split; at a prime leaf
measure its valuation in the original n, including when the split was
noncoprime. This computes the bit without a supplied factorization.
The reference uses trial division to certify probable primes beyond the
deterministic range of its primality library.

There is a worst-case obstruction to this particular scan. Fix ℓ=7 and j=1,
and take primes X<p<q<2X both equal to 1 mod 28. The prime number theorem
in arithmetic progressions supplies such pairs for all sufficiently large X.
Below √(pq) the only possible nonunit is p. Its residue 1 never satisfies
the local condition for any exponent: the local expression there is −3 mod 7
with the actual positive character. Since pq is also 1 mod 28, residue 1
cannot be the complementary class of any target either. Yet the e=1,
positive-character target 9 mod 28 is scanned up to √(pq). The reference
therefore takes Ω(√n) gcds on these no instances, including all empty probes.
A faster way of rejecting whole progressions remains a distinct problem.

### Batching the arithmetic progressions

With python-flint 0.9 installed, the native evaluator now screens a progression
of K ≥ 64 candidates in blocks of w=floor(√K); pass batched=False to reproduce
the individual scan. Construct P(X)=∏(X+i·stride), 0≤i<w, by balanced
polynomial multiplication over Z/mZ. Multipoint evaluation at every full
block start computes its product. A unit gcd clears the block, a proper gcd
splits m, and a gcd equal to m requires at most w individual gcds to find
a proper split. Every candidate is strictly between 1 and m, so that split
exists. The tail has fewer than w terms. The polynomial is monic and
multipoint remainders have monic divisors, including over a composite ring.

With fast polynomial arithmetic, this uses Õ(√K) modular arithmetic work
per progression including misses (plus bit costs polynomial in log m).
On the n=1,000,406,036,309 fixture the answer was false with 850 gcds,
774 linear factors constructed and 775 full-block values evaluated,
versus 143,301 individual gcds. Polynomial arithmetic costs additional time;
these counts do not measure wall-clock speedup. The paired-prime family
above still forces Ω(n^(1/4)) full-block values for fixed ℓ=7, since there
is no split and at least one progression has Θ(√n) candidates. This is a
barrier for the blocked progression traversal, not a general lower bound
for the bit. Without python-flint the evaluator uses the original scan.

### Where the fixed zero bit is easy on random integers

For each fixed prime ℓ≥7 and fixed j with 5j≠2 mod ℓ, there is a
polynomial-bit-time, one-sided certificate finder that succeeds on a
density-one subset of odd n coprime to ℓ. The implementation checks primes
up to the bit length of n and verifies the full valuation of any divisor.
It returns true only if a local factor is zero. Otherwise the exact evaluator
continues with the progression search.

For e=1 the target consists of s=1 or 2 invertible classes modulo 4ℓ.
At most one numerator 5j±1 vanishes, so s≥1. Each prime p in those classes
with p dividing n exactly once certifies the bit. On X≤n<2X take
Y=c log X with a sufficiently small fixed c so Y is below the bit length
and 2ℓ times the product of p² for target primes p≤Y is at most √X
(Chebyshev's bound). On a whole period CRT makes the fraction of eligible
integers missing all these witnesses exactly ∏(1−1/p+1/p²).
Truncating the final period costs O_ℓ(X^(-1/2)) in the fraction.
Mertens' theorem for fixed arithmetic progressions bounds the miss
probability by O_ℓ((log log X)^(-s/φ(4ℓ))). This bounds fallback frequency,
not expected running time. It is a fixed-query result and gives no uniform
density bound for primes growing with input length.

For ℓ=7, j=1 the exponent is only 1/6. In independent groups of 2,000
uniform odd b-bit integers coprime to 7 (Python Random seed 2701, widths
64, 256, 1024 in order), the implementation certified 92, 164, and 226
positive answers respectively. It gets no small-prime certificate on large
balanced semiprimes. At the exceptional j=2/5 mod ℓ every squarefree
input gives false. Higher valuations can change the bit: ℓ=7, j=6 and
n=13² gives true. A uniform declaration of false there is unsound.

Even supplementing the small-prime check by the product's residue modulo
any *fixed* auxiliary integer D cannot decide the ℓ=7, j=1 bit on all
balanced semiprimes. Set L=lcm(28,D). For negative instances choose two
distinct primes in (X,2X) both 1 mod L. For positives choose a unit a mod L
lifting 9 mod 28 (avoid the other prime divisors of L by CRT), then choose
primes in (X,2X) of residues a and a⁻¹ mod L. The prime number theorem in
fixed progressions supplies both constructions at arbitrarily large X.
Both products are 1 mod L and have no divisor up to their bit length;
only the positive pair contains a target prime. This obstruction concerns
fixed modular summaries, and does not exclude adaptive arithmetic on n.

### Selecting a prime guaranteed to have a positive exponent

Let b be the bit length of odd N>1. At any prime ℓ>max(5,b) coprime to N,
the coefficient 5e is a unit for every exponent e in N. The local linear
factor admits a root exactly when p≠1 mod ℓ: in the other case it is the
nonzero constant −2−χ. Thus every j misses iff *every* p dividing N is
1 mod ℓ. A public sufficient test for at least one hit is N mod ℓ outside
{0,1}. The least prime satisfying these conditions is O(b): the excluded
prime product above max(5,b) divides N(N−1), whereas the log prime product
up to Cb exceeds log(N(N−1)) plus the log prime product up to max(5,b)
for a sufficiently large absolute C. Trial division of these small
candidate primes is polynomial in b. This selects ℓ, leaving the hitting j
unknown. The decoder now skips residues 0 and 1, sacrificing at most
log(N(N−1)) of available log prime product.

The magnitude O(b) is worst-case sharp even when ℓ need only be ≥7.
Let M be the product of primes up to y. Choose distinct comparable primes
p,q equal to 1 modulo M and of size at most M^O(1). To justify comparability
without assuming a prime-gap conjecture, choose a prime r between y² and
2y². Modulo Mr there are r−1 unit classes reducing to 1 mod M. Linnik's
least-prime theorem gives distinct primes in those classes, all ≤C(Mr)^L
for absolute constants C,L. There are only O(log(Mr))=O(y) dyadic size
intervals, so two primes lie in the same interval. For N=pq one has b=Θ(y)
and every eligible ℓ≤y misses at every exponent. For a finite example take
M=97# (the product of all primes through 97), p=134M+1 and q=135M+1.
Here N has 256 bits, the first hitting prime is 101, and the public rule
selects 257, with j=132 and 224 giving positive answers modulo 257.

One can also shrink the full decoder's worst-case prime budget. Let K
bound the number k of distinct primes; the implementation takes the largest
K with 3·5···(2K+1)≤N. A local numerator and denominator both have magnitude
≤h_p=5ep. Since e≤p^(e−1), their height product over p is ≤D=5^K N.
For the distinct rational roots the coefficient bound becomes U=2^K D.
The product of all nonzero cross determinants is at most
2^(K(K−1)/2) D^(K−1), by using 2 h_p h_q for each pair. Including local
denominators and N gives B=N 2^(K(K−1)/2) D^K. These replace the coarser
bounds above. Chebyshev's bounds now put the largest auxiliary prime in the
complete reconstruction at O(Kb), hence O(b²/log(2b)) since K=O(b/log(2b)).
This bounds the sizes of oracle queries; their native evaluation cost remains.

### A single noisy hitting exponent per prime suffices

For odd n>1 let K be the public upper bound on the number of distinct
prime factors used above and b its bit length. Select primes ℓ>max(b,5)
with n mod ℓ outside {0,1}. Ask for only one exponent per prime. Fix
a positive integer h. For the selected primes with product M suppose the
product of those on which the response actually hits is at least M^(1/h).
The decoder need not recognize which responses hit. The hidden prime causing
a hit may change on every query.

Use the integral polynomial F obtained by multiplying the local linear
factors 5e(p−1)X−(2p+χ), omitting the positive scalar in the power formula.
It has degree k≤K and Euclidean coefficient norm at most its coefficient
1-norm, which is at most 10^K n by the product-of-heights bound. At each
correct response F vanishes modulo the selected prime.

Combine all answers by CRT into x mod M. Take d=2Kh and the lattice of
integer coefficient vectors (g₀,...,g_d) satisfying G(x)=0 mod M. The
congruence is surjective via its constant term, so the determinant is M.
The usual LLL determinant bound (δ=3/4) supplies a nonzero polynomial
with Euclidean norm at most 2^(d/2) M^(1/(d+1)), allowing slack.
The implementation verifies this bound by an integer inequality.

Let R be the resultant of F and G. Every successful prime divides R since
the two polynomials have a common finite root in that field (equivalently
evaluate the integral resultant Bezout identity at x). Hence if R is nonzero,
|R|≥M^(1/h). Hadamard's bound on the Sylvester matrix gives
|R|≤(10^K n)^d · 2^(dK/2) · M^(K/(d+1)).
This is strictly smaller than M^(1/h) once
log₂ M > 2hd(b+5K+2). Thus R=0. Since F splits over the rationals, G
has a rational root belonging to a true prime component. Trying the public
exponents and signs as before gives a verified proper divisor for composite
n. Recurse on the two smaller inputs; the prediction guarantee is required
on each call. Powers of two can be removed directly.

The product cutoff has polynomial bit length in b and h; Chebyshev's bounds
and the n(n−1) exclusion give polynomially sized auxiliary primes. This is a
deterministic reduction under a log-weighted agreement promise on the actual
query set. It gives no efficient native predictor. In particular, the fact
that some exponent always hits does not supply the required rate for a
uniform guess as the prime range grows. [Implementation](quartic_collision_hit.py)

### Optimal prediction from the public product residue

Consider a squarefree balanced semiprime, and condition on its eligible unit
residue v modulo 4ℓ. In the finite model choose the first prime's unit class
uniformly and set the second class to v divided by the first. The best possible
hit probability of a single prediction j is **exactly 2/(ℓ−1)** for every
prime ℓ≥7 and every eligible v. This probability concerns the distribution
of hidden classes, not a guarantee on any particular factored input.

To see the upper bound, the proposed slope t=5j determines at most two
target classes S. A hit requires the first class to belong to S or v/S.
This union has at most four elements among 2(ℓ−1) units. To attain four,
exclude t=2,1,−1. The two remaining target labels modulo ℓ are
a=(t+1)/(t−2), c=(t−1)/(t−2), with respective signs +1 and −1 modulo 4.
For v=1 mod 4 a collision between S and v/S requires a²=v or c²=v
modulo ℓ; for v=3 mod 4 it requires ac=v. Each equation excludes at most
two more slopes. The polynomials after clearing denominators cannot vanish
identically: evaluate at t=2, where the numerators are nonzero mod ℓ.
At most seven slopes in total are excluded, so for ℓ≥11 some t in 0,...,7
attains the upper bound. Direct enumeration handles ℓ=7. This constructs
a predictor with constant many modular inversions. [Predictor and exhaustive
finite-model audit](quartic_collision_local_prediction.py)

There is also a typical-input obstruction at the growing prime sizes used
by the LLL decoder. Choose distinct primes independently without replacement
in [X,1.1X], with X a power of two so the product has fixed bit length b.
For each ℓ≤Y=poly(b), allow a prediction depending arbitrarily on the product
residue modulo a predetermined m_ℓ≤poly(b), as well as on b and ℓ. In
particular the predictor may use a different modulus for each prime.
Under Siegel–Walfisz, uniformly for these moduli, each prime's class modulo
lcm(m_ℓ,4ℓ) is uniform among units up to relative error o(1) (take the
logarithmic error exponent larger than the fixed modulus-growth exponent).
Conditioning on the product class still leaves the first class asymptotically
uniform. Excluding identical primes changes probabilities negligibly.
Thus the probability of an eligible hit is at most (2+o(1))/(ℓ−1).
The expected total log weight of hits through Y is O(log Y) by the prime
harmonic sum, even if only a subset of those hits counts.

For the actual decoder with agreement parameter h≤poly(b), its endpoint
is bounded by Y=poly(b), since excluded primes divide n(n−1) or are ≤b.
Its required correct weight is at least log(M)/h=Ω(hKb), from the stated
product cutoff. Markov's inequality makes the fraction of these balanced
inputs meeting that promise O(log Y/(hKb)). This rules out a uniform
agreement guarantee for this predictor family, even with arbitrary lookup
tables. Predictors using the full integer or an accumulated transcript of
unbounded joint modulus are outside the argument. The bound concerns the
sufficient agreement promise, not all ways LLL might succeed incidentally.

### An unfactored rejection test for algebraic predictions

Suppose a full-input predictor exposes an explicit nonzero G in Z[X], and
each proposed j is a root of G mod ℓ. Its choice among roots can depend on
N, ℓ and the entire transcript. A list of public templates can be combined
by multiplication, including a linear factor for any constant fallback.
For odd composite N one can either obtain a verified split directly or
upper-bound the *total product* of distinct hitting primes over the entire run:
$$
P_{hit} \le (10^K N)^{\deg G}\,\|G\|_1^K .
$$
Here K is the public upper bound on the number of prime factors from the
one-hit decoder. Extract G's rational roots, try the existing exponent/sign
decoder, and test the candidates by division. In the absence of a split, F
and G are coprime over Q: any common root belongs to a prime component of N
and gives that prime upon the correct public exponent/sign guess.
Every hit then divides their nonzero integral resultant. Hadamard's bound
and the earlier norm bound on F give the displayed inequality without
knowing the factors. If the right side is below the required product
threshold, the prediction strategy can be rejected before evaluating a
single collision bit. This test needs an explicit degree/height bound;
large templates and predictors without such a relation remain open.

For a concrete full-input predictor, set t to the greatest odd integer
at most √N, and propose the slope for the estimated prime t (exponent one),
using j=0 when its denominator is noninvertible. The public template is
G(X)=X[5(t−1)X−(2t+χ_t)], with χ_t determined by t mod 4. For distinct
odd primes N=pq there is a sharper factor-supplied audit: define
D_r=(2+χ_t)(r−1)−(2+χ_r)(t−1). Every hitting prime divides
(t−1)D_p D_q. This follows by cross multiplication or the denominator
fallback, using ℓ≥7. If D_r=0, trying the two sign choices directly
recovers r from t. A numerically close proxy with nonzero discrepancies
has only the prime divisors of those discrepancies available for agreement.

Using the one-hit decoder's h=1 query schedule, the proxy t=101 split
10403 directly. For N=1,500,011,500,021 it hit 2 of 427 primes, with
log₂ agreement product 18.32 against total weight 4313.39. Before seeing
any answers the public resultant test bounds its agreement product to
fewer than 408 log₂ units on that input. Neither the certificate nor
the direct split requires the audit factorization.
[Public test and fixtures](quartic_collision_template.py)

### Additive precision and multiplicative amplification

Write A(m)=H(m)/m for odd m. A uniform evaluator of A to absolute error
one suffices for deterministic polynomial-time factoring. Choose a product
s of distinct primes p≥7 coprime to n until A(s)>2n. Each local ratio
3−(5+χ)/p exceeds 2. Only O(log n) primes are needed, even counting the
skipped divisors; log s=O(log n log log(3n)). For every power query in
the exact decoder one has n A(sn^j)/A(s)=H(n^j)/n^(j−1), an integer.
Scale the approximate answer and round; its error is strictly below 1/2.
One padding s serves all powers of n. The promise concerns absolute error
in A; multiplication by A(s) does not improve relative precision.
[Decoder](quartic_collision_amplify.py)

More generally let T_c(m)=exp(c log m / log log m) for sufficiently large m,
using natural logarithms and fixed c≥0. Consider evaluating A(m) with
absolute error at most T_c(m). Its leading-constant boundary is as follows.

* For c≤log 3, such an all-input evaluator suffices to factor in polynomial
  time. To recover the exact H(u), take y=(log u)² and let s contain the
  primes from 7 through y that do not divide u. Write m=su. One has
  log A(s) ≥ π(y)log 3 − O(log u + log log y): omitted divisors cost
  at most O(log u), and the local corrections to log 3 are O(1/p).
  Also log m ≤ θ(y)+log u. By the classical prime number theorem with
  remainder, θ(y)=y+o(y/log²y) and
  π(y)=y/log y+y/log²y+o(y/log²y). Hence even at c=log 3,
  log A(s)−log T_c(m) ≥ (log 3+o(1))y/log²y−O(log u).
  This exceeds log(2u) for large u. Scaling the answer by u/A(s)
  and rounding gives H(u). The padded input has O(log²u) bits.
  For a uniform implementation, grow y by doubling until the strict
  inequality holds with a fixed extra margin, using bounded-precision
  enclosures; using c=log 3 suffices throughout. Apply separately to
  the power queries. The reverse reduction uses the local formula given
  the factorization.
* For fixed c>log 3, returning zero meets the tolerance at all sufficiently
  large odd inputs. To see this put L=log m and z=L/(log L)^3.
  Each local A(p^e) is ≤5e−2≤3^e, directly from the local formula.
  For p≤z bound each by 5 log₂ m; their total logarithm is at most
  z log(5 log₂ m)=o(L/log L). For p>z the sum of exponents is ≤L/log z.
  Thus log A(m)≤(log 3+o(1))L/log L uniformly. The displayed elementary
  upper bound can be checked before returning zero; on the finitely many
  remaining inputs compute the count exactly. This is polynomial bit time
  for each fixed c with effective comparison (one can use a fixed rational
  constant strictly between log 3 and c).

This boundary locates where a coarse additive estimate becomes uninformative.
It supplies no polynomial native evaluator at unit accuracy or at c=log 3,
and it uses the freedom to query moduli with growing numbers of components.

### Relative accuracy recovers the number of components

Let A(m)=H(m)/m. On *any* odd n>1, two calls to an evaluator of A with
10% relative error determine k=ω(n) exactly. Put b=bitlength(n), J=8b.
In A(n^j), the local factor for p^e exactly dividing n is
5e(1−1/p)(j−r), where r=(2p+χ)/(5e(p−1)) lies strictly between 0 and 1.
Thus A(n^(2J))/A(n^J) equals 2^k times a product of
1+r/(2(J−r)). This correction lies between 1 and exp(k/(2(J−1))).
Independent adversarial ±10% errors in the two values change their ratio
by a factor between 9/11 and 11/9. The resulting ratio divided by 2^k
lies strictly between 1/√2 and √2, so rounding the base-two logarithm
gives k. Queried moduli have O(b²) bits. This consequence does not by
itself give a split.

In the reverse direction, on a promised squarefree odd n of known component
count k, there is a deterministic unfactored relative evaluator for every
A(n^j). With integer inverse error R≥1, trial-divide up to 2k(R+1).
For each prime found use its exact local factor. For each of the other
components use the common leading value B=5j−2. The omitted factor in the
true/estimated ratio is 1−(5j+χ)/(Bp)≥1−2/p. The product of these
corrections is at least one minus their summed losses, hence at least
R/(R+1). The estimate lies between the exact value and (1+1/R) times
that value. The bit complexity is polynomial in b, log(2j), and R
when n and j are supplied separately. For semiprimes the count k=2 is
part of the promise; no factor oracle is used. [Both algorithms](quartic_collision_relative.py)

On balanced distinct odd semiprimes p<q<2p, the same argument gives an
absolute-error bound for the public value (5j−2)² without any trial
division: its excess over A(n^j) is at most 4(5j−2)²/p. Thus for any fixed
polynomial exponent budget in b, even all the unpadded normalized power
queries have inverse-polynomial absolute accuracy on sufficiently large
balanced inputs using this public expression alone. This explains why the
extra coprime components in the preceding additive reduction matter.

There is a sharp input-length charge for that specific amplification.
To achieve A(s)>2n using a known coprime multiplier s, the uniform
maximal-order bound above forces log s=Ω(log n log log n): first log s
is Ω(log n), then apply log A(s)=O(log s/log log s). The distinct-prime
construction has log s=O(log n log log n), attaining this order. This
is a limit on the multiplicative amplification route.

### Repeated factors at relative precision

A 10%-relative oracle for A on all odd inputs computes the Möbius function
exactly on every positive integer with at most two queries. It also estimates
the product D=∏e of the prime exponents. Strip primes up to 16b, including 2,
tracking their counts and exponents; call the residual m. If m>1 use
J=16b and estimate A(m^J), A(m^(2J)). Recover k=ω(m) from their ratio as
above. Divide the first estimate by (5J)^k. Its true normalized value is
∏e · ∏(1−1/p)(1−r/J), with products over the residual primes and 0<r<1.
Because k≤b and p>16b, the combined correction is between 1−2b/(16b)=7/8
and 1. Multiplying back the stripped exponent product yields an estimate
between (63/80)D and (11/10)D, including the 10% measurement error.
The cutoff 7/5 separates D=1 from D≥2. In the former case the answer
is (−1)^ω(n), in the latter it is zero. On cube-free inputs D=2^v where
v counts primes of exponent two, so geometric rounding gives v exactly.
This recovers the complete exponent histogram on that promise, without
identifying the primes. [Decoder and profile-aware evaluator](quartic_collision_relative.py)

Given an exponent histogram on any odd input, trial division up to
2k(R+1) gives an upper estimate within relative error 1/R for A(n^j),
even with unbounded exponents. Replace the contribution of each remaining
p^e by 5ej−2. The discarded local relative loss is at most 2/p for
all ej≥1, so the earlier product bound applies unchanged.

Conversely, for any fixed exponent bound E, the histogram follows from
E+1 power queries with sufficiently small constant relative error depending
on E. To see identifiability, use the E adjacent log differences at
j=1,...,E. The leading models for a component of exponent e give entries
L_{je}=log((j+1−s_e)/(j−s_e)), s_e=2/(5e).
This matrix is invertible. Its determinant is the integral over
t_j in [j,j+1] of det(1/(t_j−s_e)), by row multilinearity.
The Cauchy determinant formula gives a fixed nonzero sign in the interior
since the t_j are strictly ordered and the s_e are distinct.

For an explicit stability argument choose a positive constant η with
4η||L^(-1)||∞<1/2. Strip primes up to 4b/η and record their valuations.
Replacing the remaining local factors by 5ej−2 incurs total log error
at most η (use −log(1−2/p)≤4/p for p≥5). Require each answer's log error
to be at most η. The observed log-difference vector then has infinity error
at most 4η from L times the unknown integer histogram. Solve and round
each coordinate, then add back the stripped histogram. A strict margin
allows numerical error using polynomially many bits; E and the matrix are
fixed. This proof gives no uniform inverse-polynomial precision bound
for growing E. Counting exponents and recovering the prime values are
separate computational objectives here.

### Growing exponents: sufficient precision and indistinguishable profiles

Let b be the bit length of the odd input and E a public bound on every prime
exponent. Use differences at the separated power pairs (2j,2j+1), j=1,...,E.
For the ideal log-difference matrix L this gives an explicit bound
$$
 \|L^{-1}\|_\infty\le D_E
 = E(2E+1)^{2E-1}(3E^2)^{E-1}.
$$
Indeed write each row as the integral of C_{je}=1/(t_j−s_e) over the
unit interval [2j,2j+1], where s_e=2/(5e). Each entry of the inverse of
the averaged matrix is the average of the corresponding entry of C inverse
with weight proportional to |det C|. This follows by integrating the
cofactor and determinant separately; the deleted row has integration
interval of length one and the determinant has constant sign.
In the explicit Cauchy inverse each entry has 2E−1 numerator differences
bounded by 2E+1. The t differences are at least one; the E−1 differences
involving distinct s are each at least 1/(3E²). Summing the entry bound
over a row gives D_E.

Set η=1/(16D_E) in the preceding stripping and rounding argument. For
E=O(log b/log log b), D_E is polynomial in b, so inverse-polynomial
relative precision, 2E queries, and polynomial total bit work outside
the oracle suffice. The query moduli have O(Eb) bits. The log matrix
and vector can be approximated rationally: entry errors below
1/(64 E b D_E) perturb the linear system applied to the true histogram
by less than 1/(64D_E), with the same strict rounding margin. These
precision bounds give a uniform algorithm given E.

There is a limitation to inference based only on the relative values of
all positive powers. Take d large, a=d ceil(log₂ d), and give one input
binomial(d,i) prime components of exponent a+i for each even i in [0,d];
give the other those components for odd i. Each has 2^(d−1) such primes
and identical total exponent (the alternating sums of 1 and i vanish).
Use the same primes in each input, and add one further common prime of
exponent one. Both numbers then have perfect-power exponent one and
different histograms.

These inputs may even have identical bit length. Let S be their common
total exponent, L=2^d, and choose all primes between 2^L and
(1+1/(4S))2^L. The prime number theorem with classical remainder gives
enough primes for large d: the relative interval width exceeds 1/L²
eventually. Both bit lengths are LS+1 since (1+1/(4S))^S<2.
Thus log b=Θ(d) and E=Θ(log b log log b).

For each power j≥1 the difference between the two idealized log means
is, up to sign, a dth finite difference of log(x−2/(5j)) at x=a.
Integrating its dth derivative over a unit d-cube bounds it in absolute
value by (d−1)!/(a−2/5)^d. The finite-prime corrections to both log means
total at most 8(2^(d−1)+1)/2^L, uniformly in j, by the local relative-loss
bound. The sum is exp(−Ω(d log log d)), smaller than b^(−C) for every
fixed C as d grows. Thus even the entire sequence of inverse-polynomially
accurate relative power means (accuracy in the original bit length, or in
query bit length for polynomially bounded powers) can have mutually
consistent answers for two different profiles. This is a measurement
limitation; the public integers are different, and algorithms using their
other arithmetic properties are not excluded. The intermediate exponent
range between the sufficient regime and this obstruction remains unresolved.
[Numerical checks](quartic_collision_precision.py)

### Localization by the collision difference

For a uniform ordered quartic collision modulo an odd composite n, taking
gcd(x−y,n) returns a proper divisor with probability at least 3/8. A pair
is unnecessary: sampling its difference modulo n suffices.

At p^e exactly dividing n let a be the probability that p divides x−y.
The diagonal mod p has p−3 smooth points and three transverse crossings
(at 0,1,−1 along the diagonal). The branch's lift counts above give
a=[p+3(e−1)(p−1)]/[(5e−2)p−5e−χ]. Direct comparison gives
1/3≤a≤3/5 for odd primes and e≥1 (for the upper bound use 4p≥15+3χ).
Under CRT the events are independent. For two distinct components,
different indicator values imply a proper gcd. Minimizing
a(1−a')+a'(1−a) on that rectangle gives 4/9. For a single p^e, e≥2,
subtract the p^e exact diagonal pairs from the divisible pairs. The
proper-gcd probability becomes 3(e−1)(p−1)/[(5e−2)p−5e−χ]≥3/8,
approaching 3/8 for e=2 and growing p. More generally the exact split
probability equals 1−n/H(n)−∏(1−a), by subtracting the n diagonal pairs
and unit differences.

Let an evaluator return relative estimates of the count with difference
in a specified half-open interval of [0,n). Start at the whole interval;
bisect, query both child weights, and choose a child proportionally to its
estimate. Zero mass must yield zero. With ε=1/(16b) and b=bitlength(n),
the approximate Bernoulli transition differs from the true one in total
variation by at most 2ε/(1−ε); this follows by bounding the approximate
child probability between the true probability times (1−ε)/(1+ε) and
its reciprocal factor. Coupling the at most b successive choices costs
at most 2bε/(1−ε), so the split probability is still above 1/5 on every
odd composite. Verify the gcd and repeat independently on failures, strip
powers of two and use polynomial-time primality testing at the leaves.
A polynomial-bit-time interval evaluator at this relative precision
therefore yields randomized expected polynomial-time factoring. This is
a conditional localization result; evaluating the whole mean relatively
does not answer the child queries. [Reference sampler](quartic_collision_sampling.py)

Even a singleton interval already measures local quadratic characters.
For a unit difference d with 4−d² also a unit, substituting t=x+y gives
t(t²−(4−d²))=0 modulo n. At each p^e the roots are simple and number
2+χ_p(4−d²). Thus the singleton weight is 3 to the power of the number
of distinct prime components on which 4−d² is a square. For arbitrary
unit a, put d=2(1−av²)/(1+av²); then 4−d²=a(4v/(1+av²))².
After removing primes ≤5, try v=1,...,5, checking gcds with n of
1±av². A proper gcd splits, and at most four v are forbidden even modulo
one prime, so a unit choice exists. On a two-prime input the singleton
weight is 9, 3 or 1 according to the two local characters. Relative
error 1/4 distinguishes them, including quadratic residuosity on
Jacobi-positive inputs. This last reduction on its own does not show
that singleton queries suffice to factor.

### Terminal intervals of polynomial width already suffice

The localized relative evaluator can omit every interval shorter than
W=max(1, floor(n^(1/28))/8 rounded down), with integer division applied
before the maximum. On the promise n=pq for distinct primes it can omit
widths below max(1, floor(n^(1/6))/8 rounded down). Bisect while length
is at least 2W. Every requested child has length at least W and terminal
intervals have lengths in [W,2W). Here is a uniform constant-probability
argument for the general cutoff.

Let V=min(e,v_p(x−y)) at a component p^e. Up to the common scale p^(e−1),
the collision masses at V=1,...,e−1 are each 3(p−1), and the mass at
V=e is p. Indeed the three diagonal crossings each contribute (p−1)
at any strict intermediate valuation, by the transverse local coordinates
and the condition p^e | uv. The exact diagonal has p^e pairs. With
D=(5e−2)p−5e−χ, the following three bounds result:
P(V=0)≥2/5 as above,
E[V/e]=[p+3(e−1)(p−1)/2]/D≥3/10,
and P(V≥ceil(e/2))=[p+3 floor(e/2)(p−1)]/D≥3/10.

For 2≤k≤6 distinct primes, require the component with the smallest
e log p to have V=0, and all others to have V≥ceil(e/2). CRT independence
gives probability at least (2/5)(3/10)^5; the resulting proper gcd is
at least n^(1/4). For k≥7 again exclude the smallest component. The
expected log gcd from the others divided by log n is at least (3/10)(6/7)
=9/35. This random variable lies in [0,1], so its chance of reaching
1/4 is at least (9/35−1/4)/(1−1/4). Independently the excluded component
has V=0 with probability ≥2/5. For a prime power e≥2, the intermediate
valuations from ceil(e/4) through e−1 include at least (e−1)/2 values,
hence have mass at least half the earlier 3/8 split bound. For a distinct
semiprime the larger prime divides the difference and the smaller does
not with probability at least (1/3)(2/5); this gives a proper gcd ≥√n.

For an interval starting at a, let X be its length minus one. In dimension
s=8 (general) or s=4 (semiprime), use coefficient vectors of n and
(x+a)x^i, 0≤i≤s−2, scaling the degree-t column by X^t. Their determinant
is n X^(s(s−1)/2). Standard LLL gives a first vector of 1-norm at most
√s 2^((s−1)/4) n^(1/s) X^((s−1)/2). The stated W gives strict bound
n^(2/s) whenever W>1; shorter trivial intervals can be checked directly.
If d=a+r has proper gcd g≥n^(2/s), the decoded polynomial at r is divisible
by g and smaller than g in absolute value. It vanishes over the integers.
Recover its integer roots in [0,X] and verify their gcds. This arithmetic
uses polynomial bit time. The reference checks the norm certificate itself.

With relative error ε=1/(16b) every positive-mass child has approximate
transition probability at least (1−ε)/(1+ε) times the exact one. A terminal
interval's path has at most b steps, so its probability is at least
(1−2ε)^b≥7/8 times its exact probability. This retains a constant chance
of a sufficiently large proper gcd even for the smallest bound above.
Fresh trials plus primality testing and recursion suffice for randomized
expected polynomial-time factoring under the coarse interval oracle.
This changes the resolution requirement; computing those interval weights
from unfactored input remains open. [Lattice extractor](quartic_collision_coarse_split.py)

### Deterministic constant-precision searches

For the original quartic assume N=pq with p<q<2p odd primes. In the integer
band √N<d≤√(2N), the only nonunit is q: p lies below the band and 2p lies
above. At a unit difference the local weight is ≤3 at each prime (also
when 4−d² vanishes), so the weight is ≤9. At d=q the q component allows
all q choices and the p component has at least one. Thus any interval of
length w≤floor(√N)/64 has mass ≤9w on a miss and ≥q>64w on a hit. Even
adversarial 50% relative error separates these cases at threshold 16w.
Cover the band by O(1) equally sized integer blocks, find the positive
block, and bisect deterministically down to the preceding sixth-root
lattice cutoff. For large N all queried widths are Ω(N^(1/6)); the small
cases can be handled directly. Only O(log N) count calls are needed. This
gives a deterministic reduction on this promise without drawing collisions.

A related all-input target comes from the subfamily x²=y² mod N of the
quartic collisions. For odd N, using d=x−y and t=x+y gives the bijective
condition dt=0 mod N, so the restricted weight at d is exactly gcd(d,N).
The repo previously considered the *complete* Pillai sum. For its prefix
G_N(T)=∑_{1≤d≤T} gcd(d,N), a 1/8-relative estimate suffices at each T.
If p is the least prime factor, G_N(T)=T for T<p. Otherwise let
k=floor(T/p)≥1. Then G_N(T)≥T+k(p−1), and T≤(k+1)p−1 implies
k(p−1)/T≥(p−1)/(2p−1)≥2/5. Thus even the bottom of the positive estimate
range exceeds 1.2T, and the top of the negative range is only 1.125T.
Binary search through floor(√N) finds the least factor or certifies primality.
Split off twos first and iterate. Only broad prefixes are needed and the
relative tolerance is fixed for every N. In the balanced promise, e.g.,
endpoints near the detected transition have magnitude Θ(√N).

The restricted family removes the extra quartic component; the all-input
claim does not presume an oracle for the original unfiltered quartic can
compute the restricted count. Evaluating these prefixes natively is open.
At T=p on a balanced semiprime there is exactly one exceptional summand.
Uniform point sampling needs order p draws just to encounter it with
constant probability; the large mean gap alone does not provide an
efficient estimator.
[Decoders](collision_prefix_gap.py) · [Repository's complete-sum result](../../notes/nativeinv_results.md)

### Factored-input approximation and the noise threshold of the prefix test

Given the prime factorization of N, every gcd prefix can be approximated in
deterministic time polynomial in log N and inverse relative error R≥1,
without listing all divisors. Gauss's totient identity gives
G(T)=∑_{d|N} φ(d) floor(T/d)=∑_{i=1}^T F(floor(T/i)), where
F(x)=∑_{d|N,d≤x} φ(d). Build a nondecreasing step function approximating F.
Start with unit mass at 1. For each new prime power p^e, replicate each
retained jump m at d as a jump φ(p^a)m at p^a d for 0≤a≤e. The uncompressed
update is a positive sum of scaled copies of the preceding function, so it
preserves its pointwise approximation factor. In sorted order retain a new
cumulative mass only when it exceeds (1+1/(2kR)) times the last retained
mass, where k is the number of distinct primes. The discarded cumulative
mass at every position is within that factor. After k stages the approximation
lies between F/(1+1/R) and F, since exp(1/(2R))≤1+1/R. The empty
factorization gives unit mass directly. At most O(1+kR log N) jumps survive
at each stage: the first mass is one and total mass is at most N. Temporary
replication by e+1 adds only a polynomial factor. Sum each retained jump m
times floor(T/d) to evaluate any prefix; pointwise guarantees on the cumulative
function pass through the positive sum over i. On the product of the first
30 primes (155 bits), 29,508 final jumps represented 1,073,741,824 divisors
at R=8. The whole-prefix estimate was 98.84% of the exact value. [Implementation](gcd_prefix_compression.py)

The other direction tolerates every fixed adversarial relative error ε<1/3
using only endpoints T≤√N. Write δ=2−(1+ε)/(1−ε)>0. Trial-divide through
1+floor(4/δ) or √N, whichever is smaller. Grow endpoints by at most a
factor 1+δ/4 (round the increment down, at least one, and cap at floor√N).
For the first endpoint T≥p, with p the least prime factor, T≤(1+δ/4)p
and 1/p<δ/4. Thus G(T)/T≥1+(p−1)/T>2−δ. The lower permitted estimate
strictly exceeds (1+ε)T, while all permitted estimates at endpoints below
p are at most (1+ε)T. Binary-search the first crossing bracket; every point
in that bracket at or above p has the same margin. Absence of a crossing
certifies primality. This costs O((log N)/δ) probes and gives deterministic
polynomial-time factoring for fixed ε, or inverse-polynomial slack.
[Decoder](collision_prefix_gap.py)

At exactly ε=1/3 there is a public simulator on the full queried endpoint
range for prime inputs and balanced semiprimes p<q<2p. Return 4T/3.
For T≤√N a balanced pair has at most the one nonunit d=p, so T≤G(T)<2T
(and equality to T on prime inputs). The simulator satisfies the tolerance
throughout. Consequently this restricted promise oracle at that accuracy
conveys no extra information beyond N; this is not a lower bound on
factoring or on oracles allowing other input moduli and larger endpoints.

### Native one-third evaluator on every semiprime prefix

The preceding simulator extends to every semiprime N and every nonnegative
endpoint T, including multiple full periods. Even N and perfect squares
give the factors immediately under this promise, so take distinct odd
primes p<q. For 0≤t<N the exact sum is
G(t)=t+(p−1)floor(t/p)+(q−1)floor(t/q).
For t≤√N the q term vanishes, and t≤G(t)≤2t. Return 4t/3.
For √N<t<N, at least one multiple of p occurs. If k=floor(t/p),
then k(p−1)≥(t−1)/2 because t≤(k+1)p−1 and (k−1)(p−2)≥0.
Also G(t)≤3t−2: if both prime terms occur each is ≤t−1, and otherwise
the same upper bound holds. Hence (3t−1)/2≤G(t)≤3t−2, and the public
value 2t−1 is within one third throughout this range.
One whole period has mass (2p−1)(2q−1) between 3N and 4N since p≥3,q≥5.
For T=hN+t return 4hN plus the tail estimate. Relative error is preserved
under addition of positive quantities. No balance promise or supplied factors
are needed. The complexity is polynomial in the input and endpoint lengths.

On balanced distinct semiprimes, a short lattice search improves the guarantee
slightly below one third throughout the earlier decoder's range t≤√N.
Let L=max(1,floor(N^(1/6))/128 rounded down). If 2t²≤N there is no factor
in the prefix. Otherwise look for p in [t−L+1,t], clipped to integers above
√(N/2). The implied reciprocal interval for q has length at most 2L+1,
since the absolute derivative of N/x on this range is at most 2. Its
candidate q≥√N is covered by the certified dimension-four lattice radius
above (for L=1 the image of a single integer has at most one integer
candidate). A split gives the exact sum. On a miss, either no p≤t is
possible, or t≤G(t)≤U=2t−L−1, since the prefix contains at most one
multiple of p. Return 2tU/(t+U); the error is at most (U−t)/(U+t),
which improves on one third by at least 2(L+1)/(9√N). Thus polynomial
time achieves a uniform margin Ω(N^(−1/3)) on these short prefixes.
This margin is exponentially small in bit length, whereas the decoder
requires polynomially bounded inverse slack for its geometric probing
bound. [Implementation](gcd_prefix_semiprime.py)

### Power queries at the one-third boundary

Given a balanced odd semiprime N=pq with p<q<2p, **all** G_{N^j}(T)
have a deterministic 1/3-relative estimate using time polynomial in the
bit lengths of N, j and T. The base N and exponent j are supplied; T may
exceed N^j. A division by N^j or enumeration to j is unnecessary.

Set B=floor(√(N/2)), U=floor(√(2N))+1, so B≤p,q≤U≤3B for B≥3.
The divisor expansion has one summand φ(d)floor(T/d) for each d=p^a q^b,
0≤a,b≤j. Each is at most T and its degree s=a+b puts d between B^s and
U^s. If U^s≤T/32 the contribution is at least
cT with c=(1−1/B)^2(31/32). Let m be the largest degree ≤2j with
B^m≤T, and v the analogous degree for U^v≤T/32. Let C(s) count exponent
pairs of degree at most s. Then C(v)cT and C(m)T sandwich the answer.
The successive shell sizes divided by s+1 are nonincreasing (equal to
one through j, then decreasing). Thus C(m)/C(v) is at most
(m+1)(m+2)/((v+1)(v+2)) for 0≤v≤m.

Here are explicit checks that the evaluator's lower and upper bounds are
within a factor two for B≥32·3^19. For 3≤m≤19 one has U^(m−1)≤T/32,
so their ratio is at most (m+2)/(mc)<2. For m≥20 take v'=floor(9m/10).
Then U^v'≤T/32 since B^(m/10)≥32·3^(9m/10) under this base bound.
The shell ratio with v' is at most (m+1)(m+2)/[(0.9m)(0.9m+1)],
again less than 2c. For m=2,j≥2, degree one contributes at least 3cT
including degree zero. Below N there are at most four possible terms;
at or above N the known pq term contributes at least (1−1/B)^2 T/2
in addition to 3cT (or all six terms already meet the T/32 test).
The bound ratio is below two in either case. For T<B² only the
base-semiprime prefix is needed. Handle j=1 by the preceding estimator.
For the finitely bounded smaller bases the program first tries the same
certificate and uses the exact divisor expansion only if it fails. Both
degree loops for the exact expansion are capped by log₂(T) for positive T.

The harmonic mean of the certified upper and lower bounds gives the
claimed relative error. This shows that querying powers alone at the
one-third tolerance need not reveal any information beyond a balanced
semiprime base already in hand. It does not apply to unrelated or
multiplicatively padded query bases. At any fixed error strictly below
one third the earlier least-factor reduction works already at j=1.
[Estimator](gcd_prefix_semiprime.py)

### Coprime padding at exactly the one-third boundary

Given a balanced odd semiprime base N, positive j, and the prime
factorization of s coprime to N, there is a deterministic one-third
relative simulator for every G_{sN^j}(T), polynomial in the binary input
lengths (including log s). Positive convolution gives the exact expression
∑_{h|s} φ(h) G_{N^j}(floor(T/h)). Direct enumeration only handles padding
with polynomially many divisors. An additional argument prevents divisor
compression from losing the exact one-third tolerance.

Put A=floor(√(N/2))+1. For any z let [l(z),u(z)] be the public bounds
on the base-power prefix. Below A it equals z. In the band A≤z<2A use
the two semiprime bounds at √N, since z<p²: u≤2l and u≤3z. Above the
band one has u≤1.9l. Indeed, below B² both p and q have occurred so
[2z−1,3z] suffices. At larger cutoffs, the constants in the preceding
degree proof give a ratio below 1.9 for all large bases: at degree 2 the
worst ratio is at most 6/(3c+(1−1/B)^2/2), at degrees 3 through 19 at
most 5/(3c), and above that the displayed shell ratio divided by c.
When j=1 exclude q² from the candidate list from the outset. Use exact
evaluation for bounded bases when the sharper certificate fails.
In every case l≥z.

Write F and F' for the exact and compressed cumulative totient masses
of s, so F'≤F≤(1+δ)F'. For any nonnegative nonincreasing kernel on h,
the weighted sum on compressed jumps is between the true sum divided by
1+δ and the true sum. This follows by summation by parts. Apply it to
the actual prefix kernel; its public lower and upper bounds need not
themselves be monotone.

If T<2A, only h=1 can have its reduced endpoint in the band. Keep that
term separately. Other endpoints have exact kernel floor(T/h).
Deleting the common unit jump leaves cumulative functions F−1,F'−1
within factor 1+4δ (δ≤1/8): if F>1 then integer F'≥2. Thus the other
terms have a ratio bound strictly below two, and so does their sum with
the separately bounded term, or at worst two if the others vanish.

If T≥2A, call an h high when floor(T/h)≥2A, and in-band when its endpoint
is in [A,2A). Give each h elementary weight φ(h)floor(T/h). Every in-band
divisor h>1 has a high parent h/r obtained by removing one prime factor,
of at least the same weight: φ(h)≤rφ(h/r), and the parent endpoint is at
least r times as large. Each parent has at most k=ω(s) possible children.
Consequently the true total elementary weight above or in the band is at
most (k+1) times its high part. Both truncated elementary kernels are
nonincreasing, so the compressed in-band elementary weight is at most
(k+1)(1+δ) times the compressed high weight. If L denotes the compressed
sum of local lower bounds, its outside-band portion is therefore at least
L/(1+6(k+1)): the in-band lower bound is at most three times its elementary
weight and the high lower bound dominates its elementary weight. For the
compressed sum U of upper bounds this implies
U≤[2−0.1/(1+6(k+1))]L. The true answer is in [L,(1+δ)U]. Taking
δ≤1/(256(k+1)) makes the enclosure ratio at most two. Return its harmonic
mean. The implementation refines only if the current certificate requires it.
[Simulator](gcd_prefix_semiprime.py)

This permits any polynomial-length sequence of such queries, including
adaptive choices of j, s and T, when the simulator is given the base,
power and padding decomposition. It does not cover unrelated composite
bases or an unknown padding factorization.

### Second moments tolerate square-root-scale multiplicative distortion

Change the observable to S_N(T)=∑_{d=1}^T gcd(d,N)^2. This counts triples
(d,t,u) with dt=du=0 modulo N and 1≤d≤T. Both equations use the same d.
For every fixed 0<c<1, two-sided multiplicative factor
F(T)=max(1,c√T) on all short prefixes T≤√N suffices for deterministic
polynomial-time factoring. Here the permitted range is [S/F, FS], rather
than additive relative-error notation.

Strip 2 and trial-divide small integers through a cutoff of order
1/c²+1/(1−c). Probe endpoints t growing by factor at most 1/c, rounding
down and capping at floor√N. Below the least prime p, S=t. At the first
probe t≥p and throughout its crossing bracket on the positive side,
t<p/c because the preceding negative probe was below p. Now S≥p²>c²t².
Even the lowest permitted answer exceeds tF(t); below p no permitted
answer does. Binary search the bracket. This charges polynomially many
probes if 1/(1−c) and 1/c are polynomial in bit length. In particular
for c=1/2 the allowable multiplicative factor increases without bound.

There is a parallel consequence for the quartic on distinct odd semiprimes:
let w_N(d) be its collision weight at difference d. Query ∑w_N(d)^2
over the prefix. For d<p, w_N(d)≤9 and at p it is ≥p. The identical
test works with factor max(1,c√T/9) and negative ceiling 81T; trial-divide
up to order 81/c² first. This counts pairs of quartic collisions with
matching difference. No reduction from the original first-moment estimate
to this second-moment oracle is assumed.

At critical factor √T, the gcd second moment has a native short-prefix
simulator on every distinct semiprime p<q. If gcd(T,N)>1 the endpoint gives
p and the exact sum. Otherwise return floor(√(T³)). For T≤√N only p
can contribute, and S=T+(p²−1)floor(T/p). When p<T and p does not divide
T, p floor(T/p)≤T−1 and p≤T−1, giving S≤T²−T+1. For T≥3 this leaves
at least one unit of rounding slack after division by √T. When no p
occurs, S=T directly. This establishes the boundary on this promise;
it does not rule out different oracle input families.

With a supplied factorization the earlier compressed divisor table works
by replacing φ with the Jordan totient J₂(p^a)=p^(2a)−p^(2a−2), for a≥1
and J₂(1)=1. The total mass bound becomes N², keeping preprocessing and
query work polynomial. Thus even though accuracy demands on the new
unfactored evaluator are weak, evaluating it with the subcritical guarantee
is factoring-equivalent (use, for example, a tighter-than-required
approximation when F>1, and direct summation for the bounded small cutoffs
where F is one or close to one).
[Second-moment decoders](squared_collision_prefix.py) ·
[Compressed factored evaluator](gcd_prefix_compression.py)
