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
