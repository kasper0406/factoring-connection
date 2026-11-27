# Recovering mixed roots from partial oracle answers

A univariate small-root lattice can finish the decoder once the oracle has confined a root to an interval of length $N^{1/2-\varepsilon}$ for any fixed $\varepsilon>0$. The same holds for an affine progression with that many terms and step coprime to $N$. Thus the nonadaptive quarter-window reduction needs only $(1/2+\varepsilon)\log_2 N+O_\varepsilon(1)$ scaled disk queries, plus its reference. The implemented five-dimensional lattice uses a wider safety margin and needs $(3/5)\log_2 N+O(1)$ scaled queries. [Quarter-window reduction](nonadaptive-quarter.md)

## Explicit lattice

For a monic quadratic $f(x)\bmod N$ and an integer bound $X$, use the coefficient vectors after substituting $x\mapsto Xx$ in
$$N^{m-i} f(x)^i,\qquad xN^{m-i}f(x)^i\quad (0\le i<m),\qquad f(x)^m.$$
Each polynomial is divisible by $N^m$ at a root of $f\bmod N$. The basis is triangular by degree and has determinant $N^{m(m+1)}X^{m(2m+1)}$. LLL and Cauchy–Schwarz give a nonzero combination whose absolute value for $|x|\le X$ is at most
$$\sqrt{2m+1}\,2^{m/2} N^{m(m+1)/(2m+1)}X^m.$$
When this is strictly less than $N^m$, every bounded modular root is an ordinary integer root of the combination. Factor that integer polynomial and verify each candidate against the original congruence. Fixed $m$ permits $X\le c_m N^{m/(2m+1)}$ for a positive constant $c_m$. The code uses $m=2$ and $X=\lfloor N^{2/5}/4\rfloor$, with direct tests for intervals of at most two integers.

For an interval or progression $u=r+ht$, substitute into $u^2-1$ and divide by the unit $h^2\bmod N$ to make the quadratic monic. In the quarter-window reduction, partial binary digits give an interval for the lower mixed root. Stop generating oracle queries once the known interval has at most $X+1$ integers. Decode by the lattice and take gcds of the roots minus 1. For the affine CRT decoder, the same argument permits stopping the auxiliary prime list once its product exceeds the window length divided by $X+1$. If the trial retained a singleton, the CRT residues then leave an arithmetic progression short enough for the lattice. Extra candidates can always be verified. [Affine CRT decoder](affine-parity.md)

## Square-root cutoff for approximate quotients

The same determinant estimate gives an explicit variable-dimension version of the [translated-interval result](single-quadratic-count.md#short-translated-intervals-via-univariate-small-roots). Set $m$ to the bit length of $N$. Taking the $m$-th root in the strict shortness condition above requires
$$X (2m+1)^{1/(2m)} 2^{1/2} N^{1/(4m+2)} < \sqrt N.$$
The three factors multiplying $X$ have product less than four: $2m+1\le 3^m$ and $N<2^m$. Thus $X=\lfloor\sqrt N/4\rfloor$ suffices when positive. The rank is linear in the input bit length, and the entries have polynomial bit length after reducing the input quadratic coefficients modulo $N$. Exact LLL and integer polynomial factorization therefore still take polynomial bit time.

For an upper-half root, a quotient interval of width $W$ restricts the root to an integer interval of span at most $W$. Divide it into blocks of $X+1$ consecutive integers, translate the monic quadratic in each block, and verify its integer roots. For every fixed $c$, this gives a polynomial-time finish even for $W=O(\sqrt N(\log N)^c)$. Small $N$ can be scanned directly. The quotient routine implements this blocking and uses the first m≥2 satisfying (3X)^(2m+1)≤N^m for each block span X≥2. This suffices because √2 (2m+1)^(1/(2m))<3. The prefix decoders retain their five-row cutoff. This estimate changes the accuracy required of a quotient predictor; finding one with sufficient coverage remains open.

## Why short windows alone do not complete the reduction

A uniform random translate retains each mixed root with probability $w/N$, even for a fixed unit multiplier. Thus with $w=O(\sqrt N)$ the probability of retaining any mixed root is at most $2^{\omega(N)}w/N=N^{-1/2+o(1)}$. Here $\omega(N)=O(\log N/\log\log N)$ follows from the product of the smallest distinct primes. This bounds the probability of successful isolation before the lattice is applied. The lattice gives a native polynomial-time evaluator for the contents of already narrow affine windows, while the reduction needs wide windows or a source of partial coordinates.

## Checks

The code matched exhaustive root enumeration on 300 bounded progressions, including a repeated-prime input. It also split three inputs end to end after querying the shortened quarter-window transcripts. The 19-bit input 471,311 took 13 scaled disk queries and one reference. The code uses SymPy for lattice reduction and integer polynomial roots. Run `python bounded_quadratic_roots.py`.

## When divisibility at every lift forces the standard shifts

Let N>1 be odd, let y be a unit square, and put f(T)=(A+T)²−y. Fix m≥1 and a nonzero H in ℤ[T]. Write s=min(m,1+deg H). If s is at most the smallest prime factor of N, the following are equivalent:

* N^m divides H(t) for every integer t satisfying f(t)=0 mod N, including every representative of each solution class.
* H belongs to the ideal (N,f)^m in ℤ[T].

Ideal membership immediately gives divisibility. For the converse, take an exact prime-power factor P=p^e of N. Write y=h² mod P. The two roots a,b=−A±h mod P differ by a unit since p is odd and h is a unit. For every integer u, CRT supplies t congruent to a+Pu modulo P^m and to a root modulo N/P. Hence H(a+Pu)=0 mod P^m. Expand H(a+Z) with integer coefficients c_j. Terms j≥m vanish modulo P^m on substitution Z=Pu. Evaluate the remaining polynomial in u at 0 through s−1. Its Vandermonde determinant is a unit modulo P^m since s≤p. Thus P^m divides P^j c_j for j<s, and any coefficients with s≤j<m are zero. This gives membership in (P,T−a)^m. The same holds at b.

The two local ideals (P,T−a) and (P,T−b) are comaximal because a−b is a unit modulo P. Their intersection is (P,f): modulo P the product of the two monic linear factors equals f, and comaximal ideals intersect in their product. Thus the intersection of their m-th powers is (P,f)^m. Across exact prime powers the ideals (P,f) are again comaximal with intersection (N,f) by CRT modulo N. Taking products of powers proves the claim.

In particular if H has polynomial degree in the input bit length, trial division through deg H+1 either finds a prime factor or certifies the interpolation condition for every m. This uses the degree of the expanded polynomial; a short arithmetic circuit can represent a polynomial of much higher degree.

Consequently H/N^m lies in ℤ[T,f/N]. For this promise, even a different proof of divisibility at every lift gives the same class of normalized polynomials. The [capacity obstruction](https://arxiv.org/abs/1605.08065) then applies to a uniform real-interval smallness certificate. Divisibility proved only on bounded representatives or on roots satisfying an extra CM condition falls outside the equivalence.

### Two opposite root classes already suffice

Keep the same hypotheses, choose a unit square root r of y, and write a=−A+r and b=−A−r. It even suffices to have N^m divide H at the 2s integers a+jN and b+jN for 0≤j<s. Each progression may start at any representative of its class. Indeed, expanding H(a+jN) in j modulo N^m gives a polynomial of degree less than s. Its Vandermonde matrix at 0 through s−1 is invertible modulo N^m. Thus the coefficient of (T−a)^i in H is divisible by N^(m−i) for i<m, giving H in (N,T−a)^m. Repeat at b. Since a−b is a unit modulo N, the two ideals are comaximal, and modulo N their product is generated by f. Their m-th powers therefore intersect in (N,f)^m. The argument uses no factorization, even in the proof.

Thus selecting one global sign pair and retaining every lift supplies no extra polynomial divisibility of this kind. Clipping the representatives to a shorter window can invalidate the interpolation premise. For a one-sided condition there is a direct first check at m=1: reduce H modulo the monic f and N to uT+v. If H vanishes at a promised root modulo N and this remainder is nonzero, gcd(u,N) either gives a proper factor or u is invertible and −v/u recovers the shifted root. The case u=0 mod N would force v=0 mod N. The construction of such an asymmetric relation is still open.

### Extracting a one-sided relation at arbitrary multiplicity

Use the same odd modulus and unit-square hypotheses. Let a be an unknown root of f modulo N, and suppose N^m divides H(a+jN) for s consecutive integers j, with s=min(m,1+deg H) at most the smallest prime factor of N. Given explicit H outside (N,f)^m, one can extract a root or a proper factor without knowing a or any of the lifts.

By repeated monic division write H=Σ L_i(T) f(T)^i with deg L_i≤1. Ignore zero coefficients and let k be the minimum of i plus the largest exponent v such that N^v divides both coefficients of L_i. Membership in (N,f)^m is equivalent to divisibility of L_i by N^(m−i) for every i<m: expansion of the generators N^(m−j)f^j times arbitrary integer polynomials gives necessity, and summing the terms gives sufficiency. Thus k<m. For each i≤k put B_i=L_i/N^(k−i), using zero for missing terms. At least one B_i is nonzero modulo N.

For t=a+jN let q_j=f(t)/N. Dividing H(t) by N^k and reducing modulo N gives Σ B_i(a) q_j^i=0 over i≤k. Modulo N the q_j form an arithmetic progression of step f′(a)=2(A+a), a unit. The polynomial in q has degree at most min(k, floor(deg H/2)), strictly less than s. Pairwise differences among the s consecutive evaluation points are units, so interpolation gives B_i(a)=0 modulo N for every i. Choose any B_i=uT+v nonzero modulo N. Then gcd(u,N) is either a proper factor or one, in which case −v/u recovers a mod N. Adding A gives a square root of y.

The decoder uses monic divisions, integer divisions and gcds, all polynomial in the explicit input size. The existence of a compact public construction of H remains the substantive requirement. Divisibility at a single clipped representative does not supply the interpolation equations.

### Evaluating a compressed auxiliary polynomial

For m bounded by a polynomial in log N, the same extractor accepts division-free arithmetic circuits over the integers without expanding H. Evaluate every gate in R=(ℤ/N^mℤ)[T]/(f^m). As f^m is monic of degree 2m, the canonical remainder uses at most 2m coefficients of O(m log N) bits. Monic reduction and ring arithmetic take polynomial bit time per gate. Lift the final remainder to an integer polynomial P and feed P to the explicit extractor. The difference H−P lies in (N^m,f^m), a subideal of (N,f)^m, so standard-ideal membership is unchanged. The difference also vanishes modulo N^m on every lift of every root of f mod N. Thus the one-sided promise passes to P. Using m lifts per sign as applicable and testing primes below m suffices independently of the expanded degree. This covers integer addition, subtraction and multiplication circuits with explicitly encoded constants.

A planted check at N=175, m=4, A=3 and y=4 squared T+A successively 128 times, reducing each gate. Call the resulting power U. For H=U c(T+A−2) f³ the decoder returned root 2 for c=1 and divisor 7 for c=7. For H=U f⁴ it returned no relation. The power has degree 2^128 before reduction. These examples exercise compressed evaluation only; the useful linear factor was supplied.

### Integer-valued auxiliary polynomials given by evaluations

Suppose a proposed auxiliary H is rational of known degree at most d,
integer-valued on the integers, and publicly evaluable at 0 through d
modulo N^m. This is a conditional interface for a CM proposal.
No evaluation oracle of this type with improved divisibility is assumed.
Newton interpolation writes H as a sum of binom(T,j) times the jth
forward difference at zero. Thus D=d! clears every denominator. If
N≤d, trial division costs polynomial time in d. Otherwise check the
integers through d against N: a nonunit gives a proper factor, or D
is a unit modulo N^m. In the latter case the differences computed
from d+1 samples reconstruct the ordinary monomial coefficients of
H modulo N^m using inverses of j! modulo N^m. Multiplying these by D
also gives the coefficients of the integer polynomial DH at the same
precision. All costs beyond the samples are polynomial in d,m,log N.

Divisibility of H on promised integer lifts is preserved on multiplying
by D. Since the standard ideal (N,f)^m contains N^m, membership of DH
is determined by those coefficient residues. Changing representatives
modulo N^m preserves all the required lift congruences as well. Hence
the previous one-sided extractor applies to a nonstandard DH under its
stated lift hypotheses, even if the public construction exposes only
integer values. This argument does not preserve an archimedean norm
bound on clearing denominators and supplies no obstruction to a
certificate valid only at one clipped representative. Neither a new
size bound nor a public CM evaluator follows from interpolation.

A concrete factorial-denominator control is
H(X)=prod_{j=0}^{m−1}(f(X)−jN)/m!, assuming gcd(m!,N)=1.
For integer X take z congruent to f(X)/N modulo m!: the numerator
is congruent to N^m times the falling factorial of z, hence divisible
by m!. At f(X)=Nt the value is N^m binom(t,m). Throughout the real
range −1<t<m the absolute binomial is less than one. For −1<t<0
bound each absolute factor by its counterpart in m!; on (j,j+1)
for 0≤j<m bound the product by (j+1)!(m−1−j)! with strict inequality.
The interior integer points are zeros. For f=X²−1 this certifies
|X|<sqrt(mN+1). It is exactly a packaged search of the m public
equations f=jN. Shifting the quotient range also leaves m such equations.
The construction and strict-size checks agreed at 3,666 integer inputs
over four moduli and admissible multiplicities up to eight.

More generally take any integer-valued rational polynomial H of exact
degree d≥1 with d!H in the standard ideal (N,f)^m for monic quadratic
f. Repeated monic division gives its f-adic expansion with coefficients
of degree below two, the i-th divisible by N^max(0,m−i). (Expand the
standard generators in this basis.) Consequently d! times the leading
coefficient of H is a nonzero multiple of C=N^max(0,m−floor(d/2)).
On an integer interval of endpoint span L≥d, a d-th finite difference
with step h=floor(L/d) has absolute value at least C h^d. A strict
bound |H|<N^m throughout the grid forces C h^d<2^d N^m, hence
h<2 sqrt(N) and L<d(2 sqrt(N)+1). Thus a polynomial-bit-degree
candidate in this class cannot improve the exponent 1/2 by a fixed
positive amount. This extends the existing integer-shift grid bound
to factorial denominators; it does not assume a bound between integers
or cover relations outside the standard ideal.

### Taking a field norm of standard CM shifts

Let K be a number field of degree r with integer ring O, let
I=(N,f) in Z[X] for the same monic quadratic f, and take a nonzero
H in K[X] of degree d≥1 with H(t) in O for every integer t. Suppose
d!H belongs to the extended ideal I^m O[X]. This hypothesis specifies
which divisibility argument is being tested. Define P as the determinant
of multiplication by H in K[X] on a rational basis of K. Then P has
degree rd (its leading coefficient is the nonzero norm of the leading
coefficient of H), and P(t) is a rational integer at every integer t.

In an integral basis, multiplication by d!H has every matrix entry in
I^m, by the integral structure constants for multiplication in O.
Its determinant is (d!)^r P and belongs to I^(rm). The multinomial
coefficient (rd)!/(d!)^r is integral, so (rd)!P belongs to the same
ideal. Thus if the product of the conjugate absolute values of H(t)
is strictly less than N^(rm) throughout the integer interval, the
preceding grid bound gives endpoint span less than rd(2 sqrt(N)+1).
In an imaginary quadratic field that product is simply |H(t)|².
For example, a Gaussian-integer shift H=U+iV yields P=U²+V².
This argument concerns the norm as a mathematical certificate and
requires no algorithm for finding an integral basis.

For a proposed asymmetric local divisor, the quantity to check instead
is v_p(Norm H(t)) = sum over primes P above p of
f(P/p) v_P(H(t)), with v_P normalized by v_P(P)=1. At a split prime
of a quadratic field, divisibility of order m at only one prime above
p guarantees only p^m in the norm. At an inert prime the same order
guarantees p^(2m). Any claimed improvement needs both this local
accounting and the product bound across all embeddings. The standard
ideal hypothesis above need not hold for such an asymmetric construction.

### Testing CM coordinate vanishing before building a lattice

For a public monic quadratic f=X²+bX+c and a Gaussian auxiliary
H=U+iV, reduce U and V modulo f and N. Write U+ZV=A(Z)X+B(Z)
in this quotient, where A,B are linear or constant in the indeterminate
Z. Compute the gcd of N with the three coefficients of

    R(Z) = B(Z)² − b A(Z) B(Z) + c A(Z)².

This is the resultant of f and U+ZV. For every prime divisor r of N,
r divides all three coefficients exactly when f,U,V share a root
over the algebraic closure of the residue field. Indeed R factors as
the product of U(t)+ZV(t) over both roots t of f, counted with
multiplicity. A product is the zero polynomial over a field exactly
when some factor is identically zero. The coefficients must be used
jointly to test this equivalence. The calculation works on public
residues and also on division-free circuits evaluated in the quadratic
quotient, without finding t.

At an inert Gaussian prime p, vanishing of H at a base-field root of
f forces both coordinates to vanish there. If at a split prime q no
root of f kills both coordinates, the gcd above splits N=pq even if
one conjugate of H vanishes at a root modulo q. One-sided vanishing at
a specified q-root alone does not exclude a simultaneous zero at the
other q-root. The no-common-root condition matters.

For example take N=91, f=X²−4, U=X−2+35 and V=X−2+7. Both vanish
at 2 modulo 7; at 2 modulo 13 only the conjugate U−5V vanishes,
with 5²=−1 modulo 13. Their difference 28 is nonzero modulo 13,
excluding any simultaneous zero there. The resultant is
(35+7Z)(31+3Z), whose coefficient gcd with 91 is 7. This example
uses known primes to plant the asymmetry. It checks the decoder and
does not construct an auxiliary for unfactored inputs.

**Which inputs can a fixed resultant construction separate?** Suppose
the monic quadratic and auxiliary coordinates are fixed integer
polynomials in the parameter n and the variable X, evaluated at n=N.
Monic remaindering and the displayed resultant commute with reduction
modulo N. Thus its three coefficients at N reduce to their values at
n=0. Write C for the nonnegative gcd of these three integer values,
allowing C=0. The decoder's answer on every N is exactly gcd(N,C).
If C=0 it returns N for every input; otherwise the exceptional primes
are precisely the prime divisors of this fixed integer. In the planted
N=91 example the three coefficients 1085,322,21 have gcd 7 outright.

In particular a fixed polynomial construction of this kind cannot make
the resultant detect all inert Gaussian primes and reject all split
ones. Factorization type can change for infinitely many primes, but the
resultant tests existence of a *common geometric root*, and that test has
only the exceptions encoded by C. For a menu shared at a size scale B,
a pair p,q of distinct primes in [B,2B] can split only if their
membership vectors among the supports of the nonzero C's differ.
Even adaptive selection within the menu cannot separate identical
vectors. Covering every pair requires all but at most one prime to
divide the product of those C's, requiring total explicit bit length
Ω(B) by the prime number theorem. This last storage estimate does
not bound the cost of representing a huge integer by a short circuit.
Formulas using digits of N, input-dependent exponents or divisions
with rounding are outside the fixed-polynomial specialization.

**Gaussian N-th powering control.** For N=rs, distinct odd primes,
sample a uniform unit z of the Gaussian algebra modulo N. Test each of
z^N=z and z^N=conjugate(z) using a gcd of N with both coordinate
differences. At r≡3 mod 4 the unit group is cyclic of size r²−1 and
conjugation is r-th powering. The numbers of local units passing the
two tests are respectively gcd(N−1,r²−1) and gcd(s−1,r²−1).
At r≡1 mod 4 the units are pairs of nonzero scalars and conjugation
swaps the pair. The counts are gcd(s−1,r−1)² and gcd(s²−1,r−1)
out of (r−1)²: the swap test determines the second coordinate from
the first and requires the first to be fixed by s²-th powering.
Full local enumeration verified both counts on all 72 ordered pairs of
distinct odd primes below 30.

For factors in [B,2B], the conjugation test hits either component with
probability O(1/B), by the two counts. For the identity test put
g=gcd(s−1,r−1). The local hit probability is at most g/(r−1)
in both cases: in the inert case gcd(N−1,r²−1)≤g(r+1), and in
the split case square the ratio. Thus proper splitting probability
at least δ implies (r−1)/g≤2/δ on some side. Since
g=gcd(N−1,r−1), the public multiple
(N−1)·lcm(1,...,ceil(2/δ)) then covers that side's scalar group
order, enabling the earlier repeated-squaring decoder. These are tests
of full coordinate equality on uniform units; they do not analyze
other sampling laws, partial-coordinate gcds, or additive identities.

**Additive Frobenius comparison.** For distinct odd primes p<q in
[B,2B], test (z+1)^N=z^N+1 by a joint coordinate gcd. Use a uniform
element of the full Gaussian algebra here, including nonunits. Write
h_r for the corresponding local scalar pass probability in F_r. For
r split the two coordinates are independent, giving h_r². For r inert
the local probability is at most (s−1)/r², where s is the other prime:
by injectivity of r-th powering the test in F_(r²) reduces to the zeros
of (X+1)^s−X^s−1, a nonzero polynomial of degree s−1 since r≠s.

A uniform scalar gcd already splits N with probability
P=h_p(1−h_q)+h_q(1−h_p). In fact P≥max(h_p,h_q)/3. To see why, at
the larger prime q the scalar polynomial has at most p−1 roots.
Away from 0,−1 also use x^p=x^(−k) with k=q−1−p≥1. Clearing
denominators gives x^k−(x+1)^k−x^k(x+1)^k=0, of degree 2k.
Including the two endpoints gives at most 2(q−p) roots. The minimum
of these two bounds is at most 2q/3. Thus h_q≤2/3; whichever local
rate is larger, the other is at most 2/3. As a check, over all pairs
in [B,2B) for B=50,200,800 the largest observed h_q values were
16/97, 46/349, and 148/1279 respectively (6,757 pairs in total).

The Gaussian split probability is at most the sum of its two local
rates, hence at most 6P+4/B. Thus inverse-polynomial success of this
Gaussian sampler on balanced semiprimes would already imply
inverse-polynomial success of the scalar additive test on the same
inputs for sufficiently large B. The bound gives no lower bound on P
on arbitrary inputs, and does not treat chosen z or individual-coordinate
gcds.

**A family where uniform additive tests fail.** For p<q in [B,2B]
write d=q−p and assume d+1<p. On F_p the exponent N acts as d+1,
so the additive defect is a nonzero polynomial of degree d. On F_q
the inverse-power bound above gives at most 2d roots including 0,−1.
Thus the scalar per-draw split probability is at most d/p+2d/q.
A public unit shift a in (x+a)^N−x^N−a^N gives the same rate by
substitution x=at. This covers a union of polynomially many prescribed
unit shifts on a uniform sample, or adaptive shifts with fresh uniform
samples. Exhaustive local counts at (1009,1163) and (10007,10771)
gave scalar splitting probabilities 0.0056704 and 0.0011279.

For all sufficiently large B there are Ω(B^(7/4)/(log B)^2) pairs in
this dyadic interval with B^(2/3)<d≤B^(3/4). Indeed partition into
O(B^(1/4)) intervals of length at most B^(3/4). By Cauchy–Schwarz
on the bucket populations and the prime number theorem, same-bucket
pairs number Ω(B^(7/4)/(log B)^2). At most O(B^(5/3)) pairs of
integers in the entire range have separation at most B^(2/3).
Subtracting leaves the claim. On these inputs the scalar test and
the preceding Gaussian full-coordinate test both have per-draw
probability O(B^(−1/4)). Ordinary multiplier-one Fermat search also
takes Ω(B^(1/3)) midpoint increments since (p+q)/2−sqrt(pq)
equals d²/[2(p+q+2sqrt(pq))]. This rules out inverse-polynomial
worst-case success for these uniform tests, without asserting that the
same inputs resist other factoring methods.

**Choosing x=floor(sqrt(N)).** This deterministic sample also misses
most inputs in a provable family. Fix 1/2<β<α<2/3, put D=B^α,
and take B≤p<q≤2B with B^β<d=q−p≤D. Set
k=(p+q)/2−floor(sqrt(pq)); the midpoint identity above gives
1≤k=O(1+D²/B). The sample reduces to v=d/2−k modulo p and
−1−w modulo q with w=d/2+k−1. Both v and w are positive for
large B (pointwise k=O(1+d²/B)). Reflection x↦−1−x preserves
zeros of the additive defect since N is odd.

For fixed d,k, a p-side hit requires p to divide the strictly positive
integer (v+1)^(d+1)−v^(d+1)−1. A q-side hit requires q to divide
w^(d−1)−(w+1)^(d−1)−[w(w+1)]^(d−1), which is strictly negative:
this is the cleared inverse-power relation, valid since 0<w<w+1<q.
Both integers have O(D log B) bits. Each has at most O(D) prime
divisors of size at least B, and either prime together with d fixes
the pair. Across all d,k this gives O(D²(1+D²/B)) pairs with a hit.
The bucket argument above gives Ω(BD/log² B) pairs in the specified
gap range (subtract O(B^(1+β)) shorter pairs). The hit fraction is
O(log² B (D/B+D³/B²)), tending to zero for α<2/3. Ordinary Fermat
search takes Ω(B^(2β−1)) steps on this same family.

Separately, exhaustive experiments on the *previous* band
B^(2/3)<d≤B^(3/4) compared this one fixed sample with four uniform
draws per pair. Counts are in ascending lexicographic pair order
with Python Random(7403) reset for each B. The last column counts
fixed-sample splits missed by both Fermat gcds on x and x+1.

| B | pairs | sqrt sample splits | uniform splits / draws | additional splits |
|---:|---:|---:|---:|---:|
| 128 | 41 | 1 | 12 / 164 | 1 |
| 512 | 383 | 1 | 18 / 1,532 | 0 |
| 2,048 | 4,019 | 5 | 53 / 16,076 | 3 |
| 8,192 | 39,064 | 5 | 108 / 156,256 | 3 |
| 32,768 | 373,899 | 20 | 248 / 1,495,596 | 17 |

The exact expected uniform split counts for one draw per pair were
2.1530 and 4.4554 at the first two scales, computed from exhaustive
local rates. The deterministic proof is specific to this sample and
gap range; the table does not estimate success on general semiprimes.

**Checking the exceptional hits against the six symmetries.** For
F(x)=(x+1)^N−x^N−1 with odd N, F(−1−x)=F(x) and
F(1/x)=F(x)/x^N. When x and x+1 are units, reflection and inversion
generate six fractional linear transforms: x, −1−x, 1/x, −1/(1+x),
−(1+x)/x and −x/(1+x). Every transformed sample has exactly the
same additive gcd with N; scanning this orbit gives no new probes.

A separate post-hoc check asks whether the *output* y=x^N at a hit
is itself one of these six transforms of the input modulo the factor.
Replaying all 373,899 pairs at B=32,768 reproduced the 17 additive
splits missed by both previous Fermat gcds. In all 17 cases no orbit
equality held modulo the recovered factor, and none of the six
cross-multiplied comparisons gave any proper gcd with N. For example
p=32,783, q=35,081 gives x=33,912, y=6,579 modulo p, and additive
gcd p. Thus these witnesses do include nontrivial preservation of
addition beyond the checked comparisons. The replay identifies no
public feature that predicts which input pairs give such a witness.

**Primitive-root inputs and the existing AKS controls.** For prime
ℓ coprime to odd N, let ζ be a primitive ℓ-th root and evaluate the
gcd with N of Res(Φ_ℓ,F), where F=(X+1)^N−X^N−1. This detects
vanishing at *any* primitive ℓ-th root modulo a hidden prime, unlike
the earlier coefficient gcd tests. The scalar is exactly the field
norm of u^N−σ_N(u), for u=1+ζ and σ_N(ζ)=ζ^N. For odd ℓ, u is
a cyclotomic unit since Φ_ℓ(−1)=1. At ℓ=3 it is a sixth root of
unity and the defect vanishes identically for N coprime to six.

At ℓ=5 write 1+ζ=−ζ³φ for a suitable choice of golden ratio
φ=(1+√5)/2 (choose ζ=exp(2πi/5)). For each prime divisor r of
squarefree N coprime to 10 the local resultant vanishes exactly when
φ^N=σ_N(φ) at some prime above r. The condition is invariant under
φ↦−1/φ since N is odd. Hence it holds in the entire reduced quadratic
algebra modulo r. For N≡±1 mod 5 this means φ^(N−1)=1; otherwise it
means φ^(N+1)=−1. With the Fibonacci sequence F_k, compute φ^k as
F_k φ+F_(k−1). A joint gcd of the two residual coefficients therefore
reproduces the resultant gcd for squarefree N. This identifies the
order-five calculation as a quadratic-unit powering test.

On the same 200 unplanted 32- and 48-bit-prime inputs as the
coefficient audit, we reused its cyclic_power(N,ℓ,1) defect. The
primitive-resultant gcd was N on every order-three case and 1 on
every order-five, seven, and eleven case. The order-five decoder
also agreed with modular Fibonacci fast doubling on all 200 inputs.
These checks cover shift one and these four orders only.

**A bound across primitive orders.** Let B≤p<q≤2B be primes,
d=q−p, and fix 1/2<β<α<1 with B^β<d≤D=B^α. For any fixed C
allow *all* odd orders 5≤l≤(log B)^C. Write u=1+ζ_l and σ_a(ζ_l)
=ζ_l^a for units a modulo l. For large B both primes are units
modulo l. Modulo any prime ideal over r=p or q, σ_r acts as r-th
powering. The original defect u^(pq)−σ_(pq)(u) is therefore the
r-th power of u^s−σ_s(u) in the residue field, with s the partner.
At p this latter expression equals u^d σ_p(u)−σ_q(u).
At q, multiplying it by u^d gives σ_q(u)−u^d σ_p(u).
Here u is nonzero at every residue-field primitive root since l is
odd and coprime to r. Consequently **both** local vanishing tests
are characterized by divisibility of the same integer

    A(l,d,p mod l) = Norm(u^d σ_p(u) − σ_q(u)).

Once d and p modulo l are fixed, q modulo l is fixed too. This
integer has O(l(d+1)) bits, since every conjugate of u has modulus
at most 2. It is nonzero throughout this gap range for large B:
in the embedding ζ_l=exp(2πi/l), |u|≥2cos(π/5)>1, every conjugate
has magnitude at least 2sin(π/(2l))≥2/l, and |σ_q(u)|≤2. Thus
|u|^d |σ_p(u)|>2. Each integer can have at most O(l(D+1)/log B)
prime divisors in the dyadic range. Either such prime and d determines
the pair, given the side. Summing over at most D gaps, l residue
classes, and all permitted l gives O(D² polylog B) pairs with even
one local hit. The prime bucket count above gives Ω(BD/log² B) pairs
in this band. The covered fraction is at most B^(α−1) polylog B and
tends to zero. This rules out a worst-case splitter based solely on
these shift-one primitive resultants, even with adaptive choice among
all the allowed orders. The identity for A is a counting device:
constructing it as written uses the unknown gap.

**Powering to any conjugate.** The preceding gap argument extends to
a polynomial-size shared menu of monic irreducible integer polynomials
f whose degrees k and coefficient bit lengths are polynomial in log B.
Exclude zero and roots of unity as characteristic-zero roots. For each
prime divisor r of N, allow a hit whenever some root x of f in the
algebraic closure of F_r has x^N again a root of f. This is more
permissive than specifying one particular conjugate.

For a gap d=q−p>2 let a_i denote the complex roots of f and form

    C(f,d) = product over all i,j,t of (a_i^d a_j − a_t).

It is a rational integer by Galois invariance and integrality. Every
factor is nonzero: the absolute logarithmic Weil heights of the three
conjugates are equal to the same positive h, while a_i^d=a_t/a_j
would imply d h≤2h. Positivity here uses Kronecker's theorem. Cauchy's
root bound gives log(2+|C|)=O(k³(d+1) log(2+max |coeff(f)|)).
No bound on the degree of the splitting field is needed.

The local root set is permuted by r-th powering, even if the reduction
has repeated roots. Thus x^N in the root set implies x^s in it for
the partner s=N/r. Both x^p and x^q are then roots and x^d x^p=x^q.
Reduce the product defining C at a prime of the splitting field above
r: all local roots occur among the reductions, so r divides C.
For a fixed d, each menu entry therefore covers O(d polylog B) pairs,
since knowing either divisor prime and its side determines the pair.
Sum over d≤D=B^α and compare with the Ω(BD/log² B) pairs in
B^β<d≤D, 1/2<β<α<1. Even the union of all menu hits covers at most
a B^(α−1) polylog B fraction. This includes shifted cyclotomic inputs
a+ζ_l when their defining polynomials meet the stated size and
non-torsion conditions and powering is tested against a conjugate.
It does not analyze input-specific polynomial construction outside a
shared menu, or tests comparing powers to nonconjugate targets.

### Powering from one algebraic root set into another

Allow distinct monic integer source and target polynomials f,g from
a shared menu of polylogarithmic size, degrees and coefficient bit
lengths at scale B. On N=pq test the gcd of N with Res(f,g(X^N)).
Equivalently, at each prime test existence of a source root whose
N-th power belongs to the target set over the algebraic closure.
The same vanishing coverage bound holds for proper gcds on the gap
band above, even for reducible polynomials.

First take a nonzero, non-torsion irreducible source component. Denote
its roots by a_i and the target roots by b_j. At the smaller prime a
hit forces divisibility of the product over i,k,j of a_i^d a_k−b_j.
At the larger prime use the product of a_i^d b_j−a_k. Indeed Frobenius
permutes each root set, so at the smaller prime x^p is a source root
and x^q a target root; the roles reverse at the larger prime. Both
products are rational integers. A triple can vanish in characteristic
zero for at most one positive gap d, since powers of a_i are distinct
and a_k is nonzero (a zero b_j never causes an exception).
The entire menu thus has only polylogarithmically many exceptional
gaps. Allowing every pair at those gaps costs O(B polylog B).
At other gaps the previous root-size bound gives O(D polylog B) bits
per product, hence O(D² polylog B) covered pairs in total. Bounds
remain uniform for source factors: their roots are among those of f.

A zero source root hits exactly when the prime divides g(0). For a
cyclotomic source component Φ_m, φ(m)≤deg f implies m≤2(deg f)²
by the prime-power product formula for φ. Hence for large B the primes
are coprime to m and N-th powering permutes the primitive m-th roots
at either prime. Existence is measured by the fixed integer Res(Φ_m,g).
If a fixed integer here is zero the full query hits both primes and
cannot split. Otherwise it has polynomially many bits in log B (bound
the product of g over roots of unity). It contributes at most
O(B polylog B) pairs across the menu. Combining the cases and dividing
by Ω(BD/log² B) gives a proper-split fraction at most
(D/B+1/D) polylog B. This is an aggregate root-existence result.
An algorithm using other information in the resultant or constructing
a different menu for each N is outside the count.

### One residue of the input-indexed cyclotomic polynomial

A benchmark outside the shared menu is V=Φ_N(1+N) mod N², taken in
[0,N²). For odd N with at least two distinct prime divisors,

    V = 1 + N φ(N)/2.

Indeed Φ_N(1)=1 for such N. Differentiating the reciprocal identity
Φ_N(X)=X^φ(N) Φ_N(1/X) at one gives 2Φ'_N(1)=φ(N)Φ_N(1).
Taylor reduction modulo N² proves the claim, with no wraparound since
0<φ(N)<N. At N=77 this yields V=2311 and φ(N)=60, hence p+q=18.
Explicit cyclotomic coefficients and modular Horner evaluation agreed
on all 143 odd N below 501 having at least two distinct prime factors.

On a distinct semiprime the one residue therefore gives deterministic
quadratic recovery. In general φ(N) supplies a multiple of every local
unit group exponent. Using random units and gcds along the successive
squarings of the odd part of φ(N) splits any modulus with at least two
distinct odd primes with probability at least one half per trial
(the first identity indices at two primes are independent and each
atom has probability at most one half). Recurse with the same exponent
or query again. Remove powers of two and handle odd prime powers by
primality and exact root tests. Conversely given the factors the
displayed residue is immediate on the stated domain. The
factorization used by standard divisor-product evaluation of Φ_N
would already suffice. The input-index notation alone does not give
an efficient circuit for the residue. This is a scalar equivalent to
totient evaluation on the stated domain, not a new factoring algorithm.

**Fixed base, small moduli.** For distinct odd primes p,q and N=pq,
let l be a prime not dividing N, M=4^l−1 and V=Φ_N(4) mod M.
V is a unit: a prime dividing both would divide 4^N−1 and 4^l−1
and hence 3, whereas Φ_N(4)=1 mod 3 since Φ_N(1)=1.
By the divisor-product identity, the public residue
3(4^N−1)/V mod M equals (4^a−1)(4^b−1) mod M for
a=p mod l, b=q mod l. No division by either nonunit 4^a−1 or
4^b−1 is required. Compute 4^N by modular powering.

This code uniquely determines the unordered pair of nonzero residues.
In cyclic base-four positions its signed coefficients are +1 at 0
and at a+b mod l, and −1 at each of a,b, counting repetitions.
For two codes the constant +1 cancels and each difference coefficient
has magnitude at most 3. Their sum is zero, so the represented
integer has absolute value strictly below 3(1+4+...+4^(l−1))=M
(the only extremal vectors have all digits of the same sign).
Equality modulo M now means integer zero, forcing every difference
digit zero by successive reduction modulo four. Since a,b are nonzero,
neither coincides with the sum position or the constant position;
the negative digits recover their multiset. A lookup over O(l²)
pairs suffices.

Accumulate a+b mod l by CRT until the prime product exceeds N.
This recovers p+q exactly and the quadratic gives both factors.
A prime dividing N encountered during selection already splits it.
The prime number theorem allows all l=O(log N), so the number of
queries and the binary length of every modulus are polynomial in
the input length. Direct divisor-product residues verified all 2,870
eligible triples from distinct odd prime factors below 80 and query
primes below 50. This is a deterministic semiprime reduction to
modular evaluation at base four. Generating these residues without
factor information is the unresolved step.

**One fixed-base query for general inputs.** Remove factors 2 and 3
first; the following identity holds for any n>2 coprime to 3.
Take L=3^k the least power of three exceeding n and M=L². Query
V=Φ_(Ln)(4) mod M. Put y=4^(L/3) mod M, z=y³ mod M, and take
Δ=(z−y) mod M. Then

    φ(n) = [2 ((V−1)/L) · (Δ/L)^(-1)] mod L.

Use least nonnegative residues before the exact divisions. Indeed
v_3(4^(3^(k−1))−1)=k by repeated cubing, and z−y has the
same valuation since y(y+1) is a unit. The cyclotomic identity for
3 not dividing n gives Φ_(Ln)(4)=Φ_n(4^L)/Φ_n(4^(L/3)).
Writing A=Φ_n(1), the constant A is a 3-adic unit (it is one or
a prime divisor of n). Reciprocity gives Φ'_n(1)/A=φ(n)/2.
Both Taylor arguments are 1 modulo L. In the ratio modulo L² all
quadratic terms vanish, leaving V=1+Δ φ(n)/2. Since φ(n)<L
this recovers the exact totient. Both index and modulus are O(n²)
in numerical size and therefore O(log n) bits. Given φ(n) use the
preceding randomized splitting reduction; this needs no additional
cyclotomic queries.

For n=77 the query index is 6237, M=6561, V=5590 and Δ=2592,
yielding 60. Modular Horner evaluation of the numerator and denominator
using explicit Φ_n coefficients checked 100 inputs, all 3≤n≤150
coprime to three plus 385 and 1225. This verifies the decoder using
reference polynomials. Conversely, on this query family the displayed
linear formula constructs V directly from φ(n). Thus the new oracle
interface has exactly the totient content even away from semiprimes;
constructing its residue in polynomial time from the index alone
remains open.

**Precision for the existing balanced lattice decoder.** The query
identity works unchanged with L=3^k smaller than n: it yields φ(n)
modulo L, hence s=p+q modulo L for n=pq. Solve (2q−s)²=s²−4n
modulo L. If the discriminant is zero, q lies in one class of stride
3^ceil(k/2). Otherwise its true valuation 2v<k is even, and the
unit square root lifts to two classes of stride 3^(k−v), since two
is invertible. This is the ternary version of the earlier supplied-sum
decoder. It makes explicit the precision loss when p and q are close
3-adically.

On the existing 100 unplanted 48-bit-prime controls, the first k
passing a uniform sufficient rank-nine lattice radius test ranged from
16 to 20 (median 17). Use the full larger-factor range from
floor(sqrt(n))+1 to floor(sqrt(2n)), its endpoint difference W, and
the conservative centered radius max(1, floor((floor(W/stride)+1)/2)).
Test the repository determinant inequality with dimension 9 and
multiplicity 4. Reference residues were computed with known p,q using
Φ_pq(t)=S_q(t^p)/S_q(t), S_q=1+...+t^(q−1), at the two ternary
Taylor arguments; both denominators are units modulo L². The decoded
sum was checked at every precision. This benchmark checked the
sufficient lattice certificate rather than rerunning LLL. The reference
computation uses the factors. At the median the oracle supplies 17
ternary digits of the totient (about 27 bits), while its raw residue
is taken modulo 3^34.

**The first ternary digit.** At L=3 the formula simplifies to
Φ_(3N)(4) ≡ 1+3φ(N) mod 9, since y=4 and z=64. For N=pq with
distinct odd primes greater than three, N≡2 mod 3 forces one factor
to be 1 mod 3, hence the residue is always 1. For N≡1 mod 3 the
factors are either both 1 or both 2 mod 3. The residue is 1 in the
first case and 4 in the second. This distinguishes exactly whether
−3 is a square modulo N: over either prime, completing the square
identifies a square root of −3 with a root of X²+X+1, which exists
precisely when the local multiplicative group has an element of order
three. CRT combines the local tests. The Jacobi symbol is +1 in both
ambiguous cases. Either decision oracle on this promise computes the
other; no extraction of a square root from the decision is asserted.

Exact integer cyclotomic quotients over all 190 pairs of primes from
5 through 79 gave 100 forced residues and 45 of each ambiguous type.
Exhaustive square testing agreed on all 90 ambiguous inputs. Compare
the [earlier dyadic sum-bit result](single-quadratic-count.md#the-first-ambiguous-sum-bit-as-quadratic-residuosity):
both low-precision formulations reach a fixed quadratic-residuosity
question. The ternary calculation identifies the specific target and
the input class on which the residue is informative.

**Fixed index at higher precision.** For N=pq with p≡1 and q≡2 mod 3,
the same V=Φ_(3N)(4) determines q mod 3^h from V mod 3^(h+1).
Write G(u)=(4^(2u)+4^u+1)/3. The cyclotomic quotient gives
V G(p)G(q)=G(1)G(N). Starting from b=2 mod 3, try the three lifts
of b modulo 3^(k+1), setting a=N/b in that ring. Keep the lift
satisfying the quotient identity modulo 3^(k+2).

For k≥1 put w=4^(3^k)−1, of ternary valuation k+1. Writing y=4^u,
the increment in G under u→u+3^k equals w(2y²+y)/3 modulo 3^(k+2).
The coefficient is 1+2u modulo 3. Thus the change on the a branch
vanishes at this precision and the three b lifts give distinct outputs,
since every G is 1 modulo 3. The same increment bound ensures that
the test is independent of representatives. Exact integer quotients
recovered q for all 100 mixed-residue pairs of primes from 5 through
79, testing every level up to 3^h>N. This refines the earlier oracle
reduction on this promise. Computing V from N alone remains open here.

**Two fixed bases cover all distinct-prime semiprimes coprime to 6.**
Let S=p+q. Use the same index 3N at both bases 4 and 16. For r=1,2
put G_r(t)=(4^(2rt)+4^(rt)+1)/3. As above the supplied cyclotomic
residue gives D_r=G_r(p)G_r(q) via G_r(N)G_r(1)/Φ_(3N)(4^r).
All divisions here are by ternary units after normalizing G.

For a candidate sum s, compute D_r(s) as the norm of G_r(T) in
T²−sT+N=0. More precisely, put U=(4^(rT)−1)/3. Its binomial expansion
has kth term 3^(k−1)/k! times the falling factorial of rT for k≥1.
These coefficients are ternary integral, vanish modulo 3 for k≥2,
and tend to zero at least as fast as valuation (k−1)/2. Thus U≡rT
modulo 3 as a restricted power series with integral coefficients.
Use G_r=1+3(U+U²). The norm is the product at the two formal roots,
so (D_r(s)−1)/3 reduces coefficientwise to
r s + r²(s²−2N) modulo 3 by trace. Quadratic reduction and the norm
formula preserve coefficientwise convergence in s.

Consequently F(s)=(D_1(s)−D_2(s))/3 equals −s+3H(s), where H is
an integral restricted power series. For any s,t, H(s)−H(t) is
divisible by s−t in the ternary integers. Hence F preserves ternary
distance, and exactly one of the three lifts matches at every digit,
including the initial digit. The two residues modulo 3^(h+1) recover
S modulo 3^h. Taking 3^h>N recovers S and then both factors by the
quadratic formula. Evaluation on a candidate uses O(h) series terms
in a rank-two algebra with polynomial bit cost. [Runnable check on
190 prime pairs](check_fixed_index_sum.py). This removes the mixed
residue restriction on the decoder. Producing its two inputs remains open.

**Can a congruence table supply the first missing digit?** For every
fixed input modulus M, mixed semiprimes with the same N mod M can
give different Φ_(3N)(4) mod 27. Take D=lcm(9,M). Choose q in
either class 2 or 5 mod 9 and impose q≡1 at each prime power in D
coprime to 3. Lift the ternary class to a unit modulo the 3-part of D.
For each resulting unit class choose primes q and p in that class
and its negative inverse, respectively, by Dirichlet's theorem.
Then N≡−1 mod D in both cases, p≡1 mod 3, and the fixed-index lift
at k=1 forces distinct outputs because q mod 9 differs and N mod 9
agrees. For a small numerical check, 623=7·89 and 299=13·23 both
equal 56 mod 81 and give outputs 19 and 10 mod 27 respectively.
This rules out a single fixed congruence table and any finite
collection of fixed input congruences (take their lcm). It leaves
algorithms using the whole input or input-dependent moduli open.

**Fixed polynomial root counts also miss the digit.** Fix any finite
list of univariate integer polynomials of positive degree and any
congruence modulus M. Even exact counts of their distinct roots modulo
N, together with N mod M, cannot determine Φ_(3N)(4) mod 27 on every
mixed semiprime. Here is a simultaneous collision argument.

Let D=lcm(M,9) and let K be the compositum of the splitting fields
and the D-th cyclotomic field. Pick a unit u mod D equal to 2 mod 9
and pick elements σ,τ of Gal(K/Q) restricting on the cyclotomic field
to u and −1/u respectively. These lifts exist by normality.
Chebotarev gives distinct good primes q,p with Frobenius conjugate to
σ,τ, and q',p' with Frobenius conjugate to their respective inverses.
Exclude primes dividing leading coefficients or discriminants of the
squarefree parts, and primes ramified in K. Each polynomial has the
same local root count at q and q', and at p and p'. A permutation and
its inverse fix precisely the same roots. CRT multiplies the local
counts. Both products N=pq and N'=p'q' are −1 mod D, but q and q'
are 2 and 5 mod 9 respectively. The k=1 lift therefore forces different
outputs. Chebotarev supplies infinitely many such pairs for each fixed
list. This argument covers existence tests as well. Polynomials chosen
from N, and algorithms using actual roots or other data, lie outside
the claim.

**Sharp conditional success rate.** Fix the polynomial list and M,
put D=lcm(M,9), and sample p,q uniformly and independently among
primes in [X,2X] in their respective mixed classes modulo 3. Condition
on pq≡−1 mod D. From the specified counts and N mod M, the best
possible limiting success probability for the digit is 2/3 when the
list may include one suitably chosen cubic.

Indeed, in the same compositum K, Chebotarev assigns each conjugacy
class weight equal to its size divided by the group order, also in
dyadic intervals as X grows. The independent pair has product weights.
Inverting both group elements preserves those weights and the
conditioning. The conditional q classes 2,5,8 mod 9 have equal mass
because each admissible unit u mod D pairs with exactly −1/u.
Inversion exchanges the first two while preserving every observable.
They require distinct digit outputs, so even a randomized rule can
succeed on at most half their combined mass. The third class gives
the upper bound 2/3 overall on this conditioned family.

For attainment, X³−3X+1 has roots ζ+ζ⁻¹ for primitive ninth roots ζ,
with inverse pairs yielding the same root. Away from 3 the polynomial
has three roots at primes congruent to ±1 mod 9 and zero otherwise.
On the conditioned family its root count modulo N is nine exactly
when q≡8 mod 9. Recognize that class and always guess 2 on the other
two. Direct local root enumeration at all 166 primes from 5 through
997 confirmed the rule on 2,294 mixed pairs with N≡−1 mod 9.
This comparison assumes access to the root count. Its sharp rate
concerns the fixed-test model and the specified conditional distribution.

**Extracting a root on the mixed promise.** The standard
root-extraction reduction applies even to an algorithm that, given N
and a random unit cube a, only returns an integer y satisfying
y³≡a mod q, where q≡2 mod 3 is hidden. Draw a uniform unit x and
supply a=x³ mod N, keeping x private. Compute gcd(y−x,N). Since
cubing is bijective modulo q, any correct y agrees with x there.
Conditioned on a, the three equally likely values of x have distinct
residues modulo p. For any correct y independent of x given a, at
most one of the three agrees modulo p. Hence the gcd gives q with
probability at least 2/3 conditional on correctness, or exactly 2/3
for an extractor that always returns a full root. Every gcd can be
verified. Exhaustive enumeration on the 100 mixed pairs from primes
5 through 79 checked all 48,440 unit-cube fibers and every choice of
full root. This identifies an already factoring-hard witness task on
this promise. It supplies no method to construct that witness.

**Overlap with public powering.** For q≡2 mod 3 the inverse of cubing
modulo q uses exponent (2q−1)/3. Substituting N for q gives an integer
e=(2N−1)/3 on the mixed promise. With a=x³ and x a uniform unit,
gcd(a^e−x,N)=gcd(x^(2(N−1))−1,N) exactly. Write
A=gcd(2(q−1),p−1)/(p−1) and B=gcd(2(p−1),q−1)/(q−1).
Cyclicity and CRT give proper split probability A+B−2AB. Exhausting
the units for all 100 mixed pairs with prime factors from 5 through
79 verified this formula. More generally any public positive integer
exponent e gives the same gcd as the direct powering test at 3e−1.
The cube substitution therefore contributes no additional splitting
power to this implementation. See [earlier exponent coverage](#fixed-polynomial-exponents-on-separated-balanced-primes).

### Evaluation from the prime factorization

Given the complete factorization of n, let R be its radical and
replace the base b by c=b^(n/R) modulo 3^h using repeated squaring.
The identity Φ_n(X)=Φ_R(X^(n/R)) absorbs all multiplicities.
For n=1 return b−1. Otherwise compute Φ_R(1+T) in
(Z/3^h Z)[T]/(T^h). Evaluation at T=c−1 is well-defined since 3
divides c−1. Start with any prime p of R other than 3 if available.
The kth coefficient of Φ_p(1+T) is binom(p,k+1). For each remaining
prime l, compute Φ_(ml)(X)=Φ_m(X^l)/Φ_m(X). The substitution
T↦(1+T)^l−1 preserves the truncation ideal. The denominator's
constant is p at the first division and 1 at every subsequent
division. Hence each division uses a unit; a lone prime 3 requires
no division. Ordinary Horner composition takes O(h³) coefficient
operations per distinct prime. The binomial coefficients used have
O(h log(n+1)) bits before reduction, and there are at most log₂ n
prime steps. This is polynomial in h and the bit lengths of n and b
even when the number of prime divisors grows.

The [implementation](factored_cyclotomic.py) agreed with independent
exact divisor recursion in 4,500 evaluations for indices 1 through
500. The factorization is an input. Combined with the reciprocal-base
totient reduction, this gives randomized polynomial-time equivalence
between general factoring and this near-one modular evaluation problem
(at polynomial precision). Strip powers of 3 before totient recovery.

### Varying coprime index multipliers at fixed precision

Fix h≥1. For distinct primes p,q above 3 and N=pq, consider all
queries Φ_(mN)(b) mod 3^h with gcd(m,N)=1 and b≡1 mod 3. For
each admissible m,b, the answer depends on the prime factors only via
their unordered residue classes modulo 3^h. This holds without a size
bound on m.

Write [u]_y=1+y+⋯+y^(u−1). For y≡1 mod 3 the geometric sum
[3^h]_y vanishes mod 3^h: factor it into the h terms
1+y^(3^j)+y^(2·3^j) for 0≤j<h, each divisible by 3.
Thus [u]_y mod 3^h is periodic in u with period 3^h.
The semiprime formula Φ_N(y)=[N]_y/([p]_y[q]_y) has unit
denominators. This proves the assertion for m=1. For coprime m
the cyclotomic divisor product gives
Φ_(mN)(b)=∏_(d|m) Φ_N(b^d)^μ(m/d). Every factor is a ternary
unit and every b^d is still 1 mod 3, giving the assertion for all m.

For example, 35=5·7 and 3599=59·61 have matching prime residues
modulo 27 and different residues modulo 81. Their outputs agreed on
all 1,179 combinations of admissible multipliers through 200 and
bases 1,4,…,25 modulo 27, checked with the factored evaluator.
The bound concerns the information supplied by these queries at
common parameters. Algorithms can also use N itself, and a multiplier
with a proper gcd with N already supplies a factor.

### Two digits at several primes

For n with at least two distinct prime factors, set t=φ(n)/2.
The derivative identity above gives Φ_n(1)=1 and Φ'_n(1)=t.
Consequently Φ_n(1+ℓ) mod ℓ² is exactly 1+ℓ(t mod ℓ) for each
prime ℓ, including 2 and any prime dividing n. These queries all
have index n. Read the second digit and combine by CRT, stopping
when the prime product reaches n/2 since 0≤t<n/2. Primality and
perfect-power checks handle the omitted input classes separately.
By the prime number theorem for the primorial, primes up to O(log n)
suffice and the number of queries is O(log n/log log n). Each
base and modulus therefore has O(log log n) bits. Even a uniform
evaluator polynomial in the *numerical* modulus and log n suffices
for a polynomial-time reduction here.

Exact divisor recursion recovered all 806 eligible totients for
3≤n≤1000 in 3,674 queries. For n=77, for example, the second digits
are 0,0,0,2 at primes 2,3,5,7, giving t=30. This is a change of
query parameters for the existing reciprocity reduction. Computing
the digit from n alone remains the missing operation.

**Wrong digits can be corrected on a semiprime.** Convert each
predicted t=φ(N)/2 mod ℓ into a predicted S=N+1−2t mod ℓ for
S=p+q. Combine these residues to r in [0,M), with M the product of
distinct queried primes. Let B be the product of the primes where
the resulting sum digit is wrong. Given S≤H, a sufficient recovery
condition is 2HB²<M. For comparable factors q≤2p one may use
H=3⌈√N⌉, and H=N works without balance.

Indeed M/B divides r−S, so (r−S)/M=k/B. Under the bound, its
distance from r/M is S/M<1/(2B²). Also r−S cannot be negative:
it is a multiple of M/B and S<M/B. The reduced fraction therefore
occurs among the continued-fraction convergents of r/M. For each
convergent a/b with b dividing M, test s=r−Ma/b by solving the
quadratic with sum s and product N and verifying the integer factors.
This uses polynomial time without knowing B or the error positions.
For N=1009·1013 and the nine odd primes through 29, all 68 tested
subsets of at most three shifted sum digits satisfying the bound
recovered S. This is a standard noisy-CRT decoder applied to the
present reduction. The guarantee requires one prediction per prime
with enough correct log-weight on the same N. A local list alone
does not meet it. [Previous CRT decoding analysis](../../notes/crtdecode_results.md).

### Evaluation at multiples of three

For n>1 and 3 dividing b, the cyclotomic divisor product reduces modulo
3^h to the product of (1−b^d)^μ(n/d) over divisors d<h. Every term
is a unit, and all omitted terms equal one. Scan integers d<h for
divisibility. This costs polynomial bit time with a Möbius oracle.
Conversely Φ_n(3)≡1−3μ(n) modulo 9, distinguishing all three Möbius
values. Thus on unrestricted indices this evaluation problem and
computing μ are polynomial-time Turing equivalent. The claim concerns
this class of bases.

On the distinct-semiprime promise the oracle can be removed. For each
observed divisor d the exponent is +1 if d=1 or n, and −1 otherwise.
At b=3, with p the smaller prime, the answer equals 1−3 through
precision h=p. At h=p+1 the inverse of 1−3^p introduces the first
change, of exact ternary valuation p. Hence on balanced semiprimes
the precision required to see a factor-dependent digit grows on the
order of √n. The [runnable check](check_nilpotent_base.py) verifies
both the truncation and the two sides of this precision threshold.

### Division gates and the quadratic norm

In the truncated ring R above, put J=(N,f)R and X=T+A. Since N^m=f^m=0, J^(2m−1)=0. Reduce a proposed denominator D modulo J to uX+v, with X²=y. Its norm is Δ=v²−yu². If gcd(Δ,N)=1, the base inverse is (v−uX)/Δ modulo J. Lift it in R by z←z(2−Dz). The error 1−Dz squares each time, so at most ceil(log₂(2m−1)) steps suffice. This extends the compact circuit evaluator to divisions by units at polynomial multiplicity.

If gcd(Δ,N) is a proper factor, return it. If Δ=0 mod N and u is a unit, −v/u is a square root of y mod N. If u has a proper gcd with N, return that factor. In the remaining case u=0 mod N, we have v²=0 mod N, so nonzero v mod N yields a proper gcd. Thus every nonunit denominator whose base residue is nonzero supplies either a root or a proper factor. A denominator reducing to zero lies in J and is nilpotent; even an exactly divisible numerator cannot generally be divided by it using only this truncated residue. Such a gate needs a cancellation identity or additional precision.

At N=175, y=4, m=4, the four denominators N+f+uX+v with (u,v)=(1,1),(1,3),(1,2),(0,0) respectively gave a verified lifted inverse, factor 5, root 173, and zero base residue. These are diagnostic inputs with chosen coefficients.

### A public power denominator

For N=pq with distinct odd primes, try D=(X+c)^(N−1)−1 in X²=y, with uniform c mod N and a fixed unit-square target. Repeated squaring takes O(log N) ring multiplications. For either prime ℓ=p,q take a root r of y modulo ℓ. The norm vanishes at ℓ exactly when c+r or c−r lies in the subgroup of (N−1)-st roots of unity in the multiplicative group. Its size is g=gcd(p−1,q−1) at both primes, since N−1 is congruent to q−1 modulo p−1 and conversely. Thus the hit probability at ℓ is the size of the union of two translates of that subgroup divided by ℓ, at most 2g/ℓ. The chance of any nonunit norm is at most 2g(1/p+1/q), uniformly in y.

There are infinitely many consecutive prime pairs with gap O(log p): the prime number theorem gives order x/log x primes in [x,2x] and the internal gaps sum to at most x. Along such a sequence g divides q−p and both primes are comparable. Even counting every nonunit as a success, the per-trial rate is O(log N/√N). This excludes a uniform inverse-polylogarithmic guarantee for uniform independent shifts and this exponent.

An exhaustive check over c with y=4 counted (unit norm, proper gcd, norm zero mod N) as (63,64,16) for N=11·13, (0,465,1426) for N=31·61, and (9603,784,16) for N=101·103. These counts concern the norm diagnostic; they do not use the known square root to split the modulus.

### Average coverage for separated balanced primes

Choose primes p in [x,2x] and q in [3x,4x] independently and uniformly. Averaging over these ranges, E[gcd(p−1,q−1)] is O(log³ x). Indeed gcd(a,b)≤Σ_{d|a,d|b} d. For each d≤2x there are at most x/d+1 integers in each interval congruent to 1 modulo d. Summing d(x/d+1)² gives O(x² log x), and the prime number theorem gives Θ(x²/log² x) prime pairs. Hence the mean hit probability of the power denominator is O(log³ N/√N) on this distribution, uniformly in the unit square chosen for each N. A union bound retains negligible mean coverage for any polynomial-in-bit-length number of uniform shifts.

There is also a bound for predetermined shifts at fixed N. Under a uniform unit-square target y, one may write y=r² with r uniform among units. For fixed c, each of c±r at ℓ=p,q hits the subgroup with probability at most g/(ℓ−1). Thus a list of L shifts chosen using N alone has success at most 2Lg(1/(p−1)+1/(q−1)), even when the most promising shift on that list is selected after inspecting y. This concerns the norm test; dependence on y that computes shifts outside such a list remains open.

### Target-dependent rational shifts

Fix N=pq and a shift c=P(y)/Q(y) with integer polynomials chosen from N alone. Suppose both degrees are at most d and Q is nonzero as a polynomial modulo each prime. For a uniform unit square y write y=r² with r a uniform unit. At ℓ=p or q, there are at most 2d possible nonzero r where Q(r²)=0. Count every pole as a potential success, including a gcd discovered while attempting division. For any of the g subgroup elements h and either sign, a hit at a regular input requires
P(r²)−hQ(r²)±rQ(r²)=0 modulo ℓ.
This polynomial is nonzero: its odd part is ±rQ(r²) while the remaining terms have even degree. Its degree is at most 2d+1. Consequently the local hit-or-pole probability is at most (2d+2g(2d+1))/(ℓ−1). Summing over the two primes bounds global coverage. Polynomial shifts use Q=1 and have no poles.

The earlier separated-prime average therefore gives O((d+1)log³ N/√N). A union bound applies even if the target selects among polynomially many such formulas after evaluating them. This argument allows a large image of the shift. Its restriction is on the degrees of P and Q; a short repeated-squaring circuit can exceed that degree range.

### Monomial shifts of arbitrary degree

For c=a y^k (a,k chosen from N alone, k≥0), the local bound can be made independent of k. Fix a prime ℓ dividing N, a sign ε=±1, and a subgroup element h≠0. Put e=2k. If r and tr are distinct nonzero solutions of a r^e+εr=h, eliminating a r^e gives ε(t−t^e)r=h(1−t^e). For t≠1 the two coefficients cannot both vanish. Thus each nonidentity ratio t occurs for at most one ordered pair of distinct solutions. A solution set of size S has S(S−1)≤ℓ−2. Summing over the g possible h and two signs gives local probability O(g/√ℓ) for uniform unit-square targets. A direct enumeration checked the inequality for all odd primes below 100, every local coefficient a, and every distinct exponent 2k modulo ℓ−1.

On the separated balanced-prime distribution the mean global norm-hit probability is O(log³ N/N^(1/4)), even when the monomial has exponential degree. A polynomial-size menu of such powers still has negligible coverage. Adding a constant b to the shift changes the right side to h−b; the proof requires that value to be nonzero.

### The exceptional value for an affine power shift

Take c=b+a y^k and put z=a y^k modulo N. At either prime ℓ, the exceptional subgroup hit h=b requires z±r=0 there. Hence ℓ divides z²−y. Compute G=gcd(z²−y,N) before the power denominator. A proper G splits N, and G=N certifies z as a root. If G=1, every remaining subgroup hit has h−b≠0, and the preceding ratio proof applies unchanged after replacing h by h−b. Consequently the *unconditional* probability of a norm hit missed by this scalar test is O(log³ N/N^(1/4)) on the same separated balanced-prime distribution with uniform unit squares. This is an additive bound, independent of k and b. It makes no smallness claim about the scalar test itself.

The scalar test has an exact local probability. Write M=(ℓ−1)/2 for the order of the local square subgroup. For a nonzero coefficient a mod ℓ it requires y^(2k−1)=a^(−2). On uniform local squares the probability is zero if the right side lies outside the image of this power map, and gcd(|2k−1|,M)/M otherwise. For a=0 mod ℓ the probability is zero. Thus a=1 maximizes the local hit probability simultaneously at both primes. Under uniform unit squares the two local events are independent by CRT. Counting either a root or a proper factor as success, the global scalar success probability is one minus the product of their miss probabilities. Finding a useful k from N alone remains open.

### Completing an exponent that covers one hidden prime

Let N=pq with distinct odd primes. Write K=|2k−1| and M=(p−1)/2. Local scalar success at least 1/B means exactly M/gcd(M,K)≤B. For an integer bound B with this property, compute E=2K lcm(1,…,B). Then p−1 divides E: the missing cofactor divides the lcm and gcd(M,K) divides K. The constructor needs only K and the bound, without knowing which prime satisfies it. If B and the bit length of K are polynomial in log N, the bit length of E and its construction cost are polynomial too.

For a fresh uniform unit a, test gcd(a^j−1,N) at j equal to the odd part of E and its successive doublings through E. If q−1 does not divide E, at most half of the units modulo q satisfy a^E=1; the endpoint alone therefore splits with probability at least one half. If q−1 divides E as well, the two local orders have independent 2-adic valuations. Each valuation distribution has maximum atom at most one half (the 2-primary part is uniform in a cyclic group of even order). Unequal valuations occur with probability at least one half and give a proper gcd at the first point where one component reaches one. Thus the index bound at one hidden prime suffices for constant splitting probability per trial using standard modular powering. Checked all unit bases for 312 prime-pair and exponent combinations, with the index bound supplied from one prime for the test.

This converts the remaining exponent-selection problem into producing a polynomial-bit integer with a polynomially bounded missing cofactor of a hidden prime's group order. A useful exponent for general inputs has yet to be constructed.

### Fixed polynomial exponents on separated balanced primes

The [earlier worst-case construction](quadratic-chirp.md#simultaneous-avoidance-of-length-indexed-polynomial-exponents) handles a broader class of exponents. For a fixed nonzero integer polynomial P there is also a mean bound on the separated balanced-prime distribution. Write K=|P(pq)| and M=(p−1)/2. Since pq≡q modulo M, gcd(M,K)=gcd(M,P(q)). For sufficiently large x, P(q)≠0 and |P(q)|≤x^C for a constant C. Its number of divisors is x^{o(1)} by the standard divisor bound. For fixed q,
the sum of gcd(M,P(q)) over eligible p is at most the sum over divisors d of |P(q)| with d≤x of d times the number of p≡1 modulo 2d. Even counting all integers in [x,2x], each summand is O(x). Dividing by the number of eligible primes and by the lower bound on M shows that the mean local hit rate is at most x^{−1+o(1)}. Interchanging the two prime ranges gives the same bound for the other prime, with d≤2x. Consequently the mean probability of either local hit on uniform unit squares is at most N^{−1/2+o(1)} for the scalar test when K is an admissible odd exponent. The same gcd bound applies without the parity restriction as a statement about group orders.

In particular the fraction of these semiprimes on which a fixed polynomial exponent has even one missing cofactor at most B is at most B times the sum of the two mean local rates. This is negligible for polylogarithmic B. The estimate concerns polynomial values as exponents; it does not apply to every polynomial-time algorithm that computes an exponent from N.

### Estimating the exponent from the square root of N

Consider E=N+1−floor(2√N)+j as an estimate for φ(N), with predetermined integer offset j. For p in [x,2x], q in [3x,4x], and M=(p−1)/2, its local coverage is gcd(M,E)/M. Modulo every divisor d of M, E equals q+1−floor(2√(pq))+j. As an integer function of q this is the ceiling of F(q)=q+1−2√(pq)+j. The derivative 1−√(p/q) is bounded above and bounded away from zero on the stated ranges. Thus each rounded value has O(1) integer preimages and the range has length O(x), uniformly in j. At most O(x/d+1) integers q give a multiple of d. Summing d times this bound over divisors of M bounds the sum of the gcds by O(x τ(M)). Restricting to primes q and dividing by their number and by M gives mean local coverage x^(−1+o(1)). Conditioning on q and varying p gives the same conclusion at the other factor: the derivative 1−√(q/p) is negative with absolute value bounded above and away from zero.

Therefore even selecting among polylogarithmically many predetermined offsets leaves mean coverage at most N^(−1/2+o(1)) on this separated-prime distribution. This applies also to small parity adjustments to obtain an odd scalar exponent. It uses the separation of the two prime ranges; when p and q nearly coincide, the derivative bound degenerates and a square-root estimate can be accurate.

### Balanced primes in the same interval

Choose distinct primes p,q uniformly without replacement from [x,2x]. For the same E and any predetermined offset j, each mean local coverage is at most x^(−1/2+o(1)). Fix p and discard |q−p|≤√x, costing at most O(log x/√x) in mean coverage even if all discarded pairs succeed. On each remaining side the absolute derivative of q+1−2√(pq)+j is |q−p|/(√q(√q+√p)), hence at least a constant times 1/√x. Each rounded value has O(√x) preimages; its range spans O(x) values independently of j. For every divisor d of M the number of remaining q giving a multiple of d is at most O(√x(x/d+1)). Sum d times this estimate over divisors of M and divide by M and the number of eligible primes. The result is x^(−1/2+o(1)); symmetry handles q. A polylogarithmic offset menu preserves the exponent.

The exponent is attained up to subpolynomial factors for j=−1. Partition [x,2x] into O(√x) intervals of length at most √x. If S is the total prime count and B the number of intervals, Cauchy–Schwarz gives at least S²/B−S ordered pairs of distinct primes in the same interval. They have probability Ω(1/√x) since S is of order x/log x. For each such pair, 0<(√q−√p)²<1, so floor(2√N)=p+q−1 and E with j=−1 equals φ(N). Both local coverages are exactly one. Thus the mean exponent N^(−1/4+o(1)) for this offset reflects genuine close-factor inputs. This lower bound is for the group-order exponent E, which may be even, rather than the odd scalar K=|2k−1|.

### Removing the exact totient hits

For the offset −1 take E=N−floor(2√N). Its error t=E−φ(N) is exactly floor((√q−√p)²), since pq is nonsquare. Exact hits are publicly recognizable: set S=N+1−E and test whether S²−4N is a nonnegative square yielding integer factors with sum S and product N.

For distinct primes from [x,2x], the mean local coverage contributed by t>0 is x^(−1+o(1)) at the level of exponents. Fix p and M=(p−1)/2, so gcd(E,M)=gcd(t,M). Split q into dyadic shells H≤|q−p|<2H on either side. On each shell t≤CH²/x, and the absolute derivative of the unrounded error is bounded below by cH/x for absolute positive constants c,C. Each integer value therefore has O(x/H) preimages. For d dividing M, the number of **positive** multiples of d in range is at most CH²/(xd). Multiplying yields O(H/d) possible integers q, and summing over the shells gives O(x/d). This removes the extra cost of the zero multiple near the stationary point. Weighting by d and summing over divisors bounds the total gcd on t>0 by O(x τ(M)). Dividing by M and the prime count yields x^(−1+o(1)) as an upper bound. The matching exponent follows from coverage at least 1/M on t>0, whose probability tends to one since t=0 requires |q−p|=O(√x). Symmetry handles the other factor. Thus the contribution beyond the exact totient test is N^(−1/2+o(1)).

### Lattice backend check

Install python-flint for larger variable-rank runs. The solver uses FLINT when available and verifies the resulting short vector's coefficient bound before recovering roots. At span floor(√N/8), the pure Python SymPy LLL backend raised an internal assertion on a 64-bit planted instance. FLINT 0.9.0 passed that input and the 40-, 48-, 80- and 96-bit cases. To reproduce, start Random(18), iterate those bit sizes in ascending order including 64, take N=nextprime(2^(b/2)−80) times nextprime(2^(b/2)−40), choose t=r.randrange(span) then l=r.randrange(N), and solve T²+lT−t(t+l)=0 mod N on [0,span]. These are local recovery checks with planted roots. At the narrower five-row cutoff the pure Python backend recovered 18 planted roots across 32 to 256 bits.
