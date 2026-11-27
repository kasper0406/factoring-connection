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
