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
